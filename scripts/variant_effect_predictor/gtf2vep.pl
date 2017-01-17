#!/usr/bin/perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


use strict;

use Getopt::Long;
use FileHandle;
use File::Path qw(make_path);
use Storable qw(nstore_fd);
use Scalar::Util qw(weaken);
use Storable qw(nstore_fd fd_retrieve freeze thaw);
use MIME::Base64;

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;

our $VERSION = 87;


# set defaults
my $config = {
  cache_region_size => 1000000,
  write_cache       => 1,
  dir               => $ENV{HOME}.'/.vep/',
  compress          => 'gzip -dc',

  # db opts
  host              => 'ensembldb.ensembl.org',
  port              => 3306,
  user              => 'anonymous',
  password          => undef,
  db_version        => $VERSION,
};

my $count_args = scalar @ARGV;

GetOptions(
  $config,
  'input|i|gtf|gff|g=s',
  'fasta|f=s',
  'species|s=s',
  'db_version|d=i',
  'cache_region_size=i',
  'dir=s',
  'help',
  'verbose',
  'synonyms=s',
  'source_biotype',
  'no_transcripts',

  'host=s',
  'port=s',
  'user=s',
  'password=s',
) or die "ERROR: Failed to parse command line options\n";

if(defined($config->{help}) || !$count_args) {
  usage();
  exit(0);
}

# check for errors
die "ERROR: No species specified\n" unless defined($config->{species});
die "ERROR: No whitespace allowed in species name\n" if $config->{species} =~ /\s/;
die "ERROR: No DB version specified\n" unless defined($config->{db_version});
die "ERROR: No FASTA file/directory specified\n" unless defined($config->{fasta});


$config->{dir} .= $config->{species}.'/'.$config->{db_version};

setup_fasta($config);

# create a coord system
$config->{coord_system} = Bio::EnsEMBL::CoordSystem->new(
  -NAME => 'chromosome',
  -RANK => 1,
);

# read synonyms from file if given
read_synonyms_file($config) if $config->{synonyms};

my @fields = qw(seqname source feature start end score strand frame attributes comments);
my $line_num = 0;
$config->{dbID} = 1;
my ($prev_chr, $by_region);
my $in_file_handle = new FileHandle;

if(defined($config->{input})) {

  # check defined input file exists
  die("ERROR: Could not find input file ", $config->{input}, "\n") unless -e $config->{input};

  if($config->{input} =~ /\.gz$/){
    $in_file_handle->open($config->{compress}." ". $config->{input} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
  }
  else {
    $in_file_handle->open( $config->{input} ) or die("ERROR: Could not read from input file ", $config->{input}, "\n");
  }
}

# no file specified - try to read data off command line
else {
  $in_file_handle = 'STDIN';
  debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
}


while(<$in_file_handle>) {
  chomp;

  next if $_ =~ /^#/; #skip lines starting with comments

  my @split = split /\t/, $_;

  my $data;

  # parse split data into hash
  for my $i(0..$#split) {
    $data->{$fields[$i]} = $split[$i];
  }

  # check chr name exists
  $data->{seqname} =~ s/chr//ig if !$config->{fasta_db}->length($data->{seqname});

  # check chr synonyms
  unless(defined($config->{fasta_db}->length($data->{seqname})) && $config->{fasta_db}->length($data->{seqname}) > 0) {
    my $synonyms = get_seq_region_synonyms($config);
    $data->{seqname} = $synonyms->{$data->{seqname}} if $synonyms->{$data->{seqname}};
  }

  unless(defined($config->{fasta_db}->length($data->{seqname})) && $config->{fasta_db}->length($data->{seqname}) > 0) {
    warn("WARNING: Could not find chromosome named ".$data->{seqname}." in FASTA file\n") unless $config->{missing_chromosomes}->{$data->{seqname}} || $config->{verbose};
    $config->{missing_chromosomes}->{$data->{seqname}} = 1;
    next;
  }

  # parse attributes
  if(defined($data->{attributes})) {
    $data->{attributes} =~ s/^\s+//g;

    my %attribs;

    foreach my $pair(split /;\s*/, $data->{attributes}) {
      my ($key, $value);

      if($pair =~ /\=/) {
        ($key, $value) = split /\=/, $pair;
      }
      else {
        ($key, $value) = split /\s+/, $pair;
      }
      next unless defined($key) and defined($value);

      # remove quote marks
      $value =~ s/\"//g;

      # lowercase key to reduce chances of mess up!
      $attribs{lc($key)} = $value;
    }

    $data->{attributes} = \%attribs;
  }

  # dump if into new region or new chrom
  # this of course assumes input file is in chrom order!!!
  if(defined($prev_chr) && $prev_chr ne $data->{seqname}) {
    debug("Dumping data for chromosome ".$prev_chr) if $config->{verbose};
    export_data($config, $prev_chr);
  }

  debug("Processing chromosome ".$data->{seqname}) if !defined($prev_chr) || $data->{seqname} ne $prev_chr;

  my $ref = parse_data($config, $data);

  $prev_chr = $data->{seqname};
}

# dump remaining transcripts
debug("Dumping data for chromosome ".$prev_chr) if $config->{verbose};
export_data($config, $prev_chr);

debug("All done!");

sub parse_data {
  my ($config, $data) = @_;

  return unless defined($data);

  # check defined feature type
  if(!defined($data->{feature})) {
    debug("Feature type not described on line $.\n") if $config->{verbose};
    return;
  }

  my $method_ref;
  unless($method_ref = function_exists('parse_'.lc($data->{feature}))) {
    debug("Cannot parse feature type ".$data->{feature}."\n") if $config->{verbose};
    return;
  }


  # run data fix
  fix_data($config, $data);

  # create transcript if not already done
  # $config->{transcripts}->{$data->{attributes}->{transcript_id}} ||= create_transcript($config, $data);

  #print "running $method\n";

  return &$method_ref($config, $data);
}

sub function_exists {
  no strict 'refs';
  my $funcname = shift;
  return defined &{$funcname} ? \&{$funcname} : undef;
}

sub fix_data {
  my ($config, $data) = @_;

  # fix strand
  $data->{strand} = $data->{strand} =~ /\-/ ? -1 : 1;
}

sub parse_gene {
  my ($config, $data) = @_;

  my $id = $data->{attributes}->{id} || $data->{attributes}->{gene_id} || $data->{attributes}->{name};
  $config->{genes}->{$id} = $data;

  # pseudogenes in NCBI RefSeq GFFs dont have transcript entries, so create one here
  if(($data->{attributes}->{pseudo} || '') eq 'true') {
    return parse_transcript($config, $data, 'pseudogene');
  }
}

sub parse_transcript {
  my ($config, $data, $biotype) = @_;

  my $id = $data->{attributes}->{transcript_id} || $data->{attributes}->{name} || $data->{attributes}->{id};

  if(!exists($config->{transcripts}->{$id})) {
    my $tr = create_transcript($config, $data, $biotype);
    $config->{transcripts}->{$tr->stable_id} = $tr if $tr;

    # store by _gff_id if exists
    $config->{transcripts_by_gff_id}->{$tr->{_gff_id}} = $tr if defined($tr->{_gff_id});
  }

  return $config->{transcripts}->{$id};
}

## NCBI RefSeq GFFs, different feature types correspond to different biotypes
sub parse_mrna               { return parse_transcript(@_, 'protein_coding'); }
sub parse_rrna               { return parse_transcript(@_, 'rRNA'); }
sub parse_trna               { return parse_transcript(@_, 'tRNA'); }
sub parse_ncrna              { return parse_transcript(@_, $_[1]->{attributes}->{ncrna_class}); }
sub parse_primary_transcript { return parse_transcript(@_); }
sub parse_c_gene_segment     { return parse_gene_segment(@_, 'c'); }
sub parse_d_gene_segment     { return parse_gene_segment(@_, 'd'); }
sub parse_vd_gene_segment    { return parse_gene_segment(@_, 'd'); }
sub parse_j_gene_segment     { return parse_gene_segment(@_, 'j'); }
sub parse_v_gene_segment     { return parse_gene_segment(@_, 'v'); }

## Ensembl special types
# genes
sub parse_lincrna_gene { return parse_gene(@_); }
sub parse_mirna_gene   { return parse_gene(@_); }
sub parse_rrna_gene    { return parse_gene(@_); }
sub parse_snorna_gene  { return parse_gene(@_); }
sub parse_snrna_gene   { return parse_gene(@_); }

# transcripts
# the type field corresponds to the biotype in the attributes
# so don't technically need to pass as an argument, but can't hurt to save time
sub parse_lincrna { return parse_transcript(@_, 'lincRNA'); }
sub parse_mirna   { return parse_transcript(@_, 'miRNA'); }
sub parse_snorna  { return parse_transcript(@_, 'snoRNA'); }
sub parse_snrna   { return parse_transcript(@_, 'snRNA'); }
sub parse_nmd_transcript_variant { return parse_transcript(@_, 'nonsense_mediated_decay'); }

# these ones might vary on biotype; just trust we can find biotype in attribs
sub parse_pseudogene { return parse_transcript(@_); }
sub parse_processed_pseudogene { return parse_transcript(@_); }
sub parse_pseudogenic_transcript { return parse_transcript(@_); }
sub parse_aberrant_processed_transcript { return parse_transcript(@_); }
sub parse_nc_primary_transcript { return parse_transcript(@_); }

# this can be either a gene or a transcript
sub parse_processed_transcript { return is_gene($_[1]) ? parse_gene(@_) : parse_transcript(@_); }

# Ensembl GFF has *_gene_segment entries that are both "genes" and "transcripts"
sub is_gene {
  return ($_[0]->{attributes}->{id} || '') =~ /^gene:/i;
}

sub parse_gene_segment {
  my $type = pop @_;

  my %gene_segment_types = map {$_ => 1} qw(c d j v);
  die("ERROR: No valid type given\n") unless $type && $gene_segment_types{lc($type)};

  return is_gene($_[1]) ? parse_gene(@_) : parse_transcript(@_, 'IG_'.uc($type).'_gene');
}

# creates a new transcript object
sub create_transcript {
  my ($config, $data, $biotype) = @_;

  # get gene
  my $parent = $data->{attributes}->{parent};
  my $gene;

  if($parent) {
    $gene = $config->{genes}->{$parent};

    # Ensembl GFF has parent ID with gene version attached
    # but gene line's ID doesn't have version on it
    if(!$gene) {
      $parent =~ s/\.\d+$//;
      $gene = $config->{genes}->{$parent};
    }
  }

  # Ensembl GTF has transcript_biotype field, thanks!
  if(!$biotype) {
    $biotype = $data->{attributes}->{transcript_biotype} || $data->{attributes}->{biotype};
  }

  # NCBI RefSeq GFF
  if(!$biotype) {
    if($gene && ($gene->{attributes}->{description} || '') =~ /^microRNA/) {
      $biotype = 'miRNA';
    }

    # otherwise get biotype from gbkey
    else {
      $biotype = $data->{attributes}->{gbkey};
    }
  }

  if(!$biotype && $config->{source_biotype}) {
    $biotype = $data->{source};
  }

  # don't bother creating a transcript unless we know the biotype
  # will only create issues later
  if(!$biotype) {
    debug("No biotype identified on line $.\n") if $config->{verbose};
    return;
  }

  # usually transcript_id is an attribute
  # pseudogenes have no transcript entry so we use the name attribute
  my $transcript_id = $data->{attributes}->{transcript_id} || $data->{attributes}->{name};

  # IG_*_genes use GeneID from dbxref
  if($biotype =~ /^IG_/ && ($data->{attributes}->{dbxref} || '') =~ /GeneID:(\d+)/) {
    $transcript_id = $1;
  }

  if(!$transcript_id) {
    debug("No transcript ID found on line $.\n") if $config->{verbose};
    return;
  }

  # print "TRANSCRIPT $transcript_id $biotype \n";

  my $transcript = new Bio::EnsEMBL::Transcript(
    -STABLE_ID => $transcript_id,
    -BIOTYPE   => $biotype,
    -SLICE     => get_slice($config, $data->{seqname}),
    -STRAND    => $data->{strand},
    -VERSION   => 1,
    -dbID      => $config->{dbID}++,
  );


  # store some private attribs we use either later or in VEP
  $transcript->{_gff_id} = $data->{attributes}->{id};

  # gene ID
  $transcript->{_gene_stable_id} = $data->{attributes}->{gene_id};

  if(!$transcript->{_gene_stable_id}) {
    foreach my $pair(split(',', $data->{attributes}->{dbxref} || '')) {
      my ($k, $v) = split(':', $pair);
      if($k eq 'GeneID') {
        $transcript->{_gene_stable_id} = $v;
        last;
      }
    }
  }

  if(!$transcript->{_gene_stable_id} && $gene) {
    $transcript->{_gene_stable_id} = $gene->{attributes}->{gene_id} if $gene->{attributes}->{gene_id};
  }

  # try and get the gene symbol
  for my $key(qw(gene_name gene)) {
    if($data->{attributes}->{$key}) {
      $transcript->{_gene_symbol} = $data->{attributes}->{$key};
      last;
    }
  }

  if(!$transcript->{_gene_symbol} && $gene) {
    $transcript->{_gene_symbol} = $gene->{attributes}->{name} if $gene->{attributes}->{name};
  }

  return $transcript;
}

sub fetch_transcript {
  my ($config, $data) = @_;

  my $attribs = $data->{attributes};
  die("ERROR: Unable to look up transcript without attributes\n") unless scalar keys %$attribs;

  # have transcript_id, easy
  if(my $tr_id = $attribs->{transcript_id}) {
    return $config->{transcripts}->{$tr_id};
  }

  # otherwise look up by _gff_id
  elsif(my $gff_id = $attribs->{parent}) {
    return $config->{transcripts_by_gff_id}->{$gff_id};
  }

  return undef;
}


# creates a new exon object
sub parse_exon {
  my ($config, $data) = @_;

  my $tr = fetch_transcript($config, $data);

  # older spec Ensembl GTFs don't have explicit transcript lines
  if(!$tr && $config->{no_transcripts}) {
    $tr = parse_transcript($config, $data);
  }

  unless($tr) {
    debug("Could not fetch transcript for exon on line $.\n") if $config->{verbose};
    return;
  }

  my $exon = new Bio::EnsEMBL::Exon(
    -START  => $data->{start},
    -END    => $data->{end},
    -STRAND => $data->{strand},
    -SLICE  => get_slice($config, $data->{seqname}),
    -PHASE  => -1,    # default phase to -1
  );

  # hidden number field
  $exon->{_number} = $data->{attributes}->{exon_number};

  # get sequence
  if(defined($config->{fasta_db})) {
    my $seq;
    if($config->{fasta_db}->isa('Bio::DB::HTS::Faidx')) {
      $seq = ($config->{fasta_db}->get_sequence($data->{seqname}.':'.$data->{start}.'-'.$data->{end}))[0];
    }
    else {
      $seq = $config->{fasta_db}->seq($data->{seqname}, $data->{start} => $data->{end});
    }

    reverse_comp(\$seq) if $data->{strand} < 0;
    $exon->{_seq_cache} = $seq;
  }

  # add it to the transcript
  # sometimes this can fail if the coordinates overlap
  eval {$tr->add_Exon($exon);};
  if($@) {
    warn("WARNING: Failed to add exon to transcript ".$tr->stable_id."\n$@");

    # delete the transcript so we don't try to process it further
    delete $config->{transcripts}->{$tr->stable_id};
    return;
  }

  return $exon;
}

# modifies a transcript with coding start/end info
sub parse_cds {
  my ($config, $data) = @_;

  # update the coding_region_start/end
  my $tr = fetch_transcript($config, $data);

  # older spec Ensembl GTFs don't have explicit transcript lines
  if(!$tr && $config->{no_transcripts}) {
    $tr = parse_transcript($config, $data);
  }

  unless($tr) {
    debug("Could not fetch transcript for CDS on line $.\n") if $config->{verbose};
    return;
  }

  # check transcript has exons
  return unless defined($tr->{_trans_exon_array});

  # update/create the translation
  if(!exists($tr->{translation})) {
    $tr->{translation} ||= new Bio::EnsEMBL::Translation(
      -TRANSCRIPT => $tr,
      -VERSION    => 1,
    );

    weaken($tr->{translation}->{transcript});
  }

  my $translation = $tr->{translation};

  # return unless $tr->biotype eq 'protein_coding';

  # get overlapping exon
  my ($matched_exon) = grep {overlap($_->start, $_->end, $data->{start}, $data->{end})} @{$tr->get_all_Exons};

  unless($matched_exon) {
    debug("Could not fetch exon matching CDS on line $.\n") if $config->{verbose};
    return;
  }

  $translation->start_Exon($matched_exon) unless defined($translation->start_Exon);
  $translation->end_Exon($matched_exon);

  my ($start_exon, $end_exon) = ($translation->start_Exon, $translation->end_Exon);

  # overlap start exon?
  if(overlap($start_exon->start, $start_exon->end, $data->{start}, $data->{end})) {
    my $offset;

    if($data->{strand} > 0) {
      $offset = ($data->{start} - $start_exon->start) + 1;
    }
    else {
      $offset = ($start_exon->end - $data->{end}) + 1;
    }

    $translation->start($offset);

    # update phase
    $start_exon->phase(frame_to_phase($data->{frame}));
  }

  # overlap end exon?
  if(overlap($end_exon->start, $end_exon->end, $data->{start}, $data->{end})) {
    my $offset;

    if($data->{strand} > 0) {
      $offset = ($data->{end} - $end_exon->start) + 1;
    }
    else {
      $offset = ($end_exon->end - $data->{start}) + 1;
    }

    $translation->end($offset);

    # update phase
    $end_exon->phase(frame_to_phase($data->{frame}));
  }
}

# converts GTF/GFF frame to phase
# basically 0 => 0, 1 => 2, 2 => 1; why? Because.
sub frame_to_phase {
  my $phase = shift;

  if($phase == 1) {
    $phase = 2;
  }
  elsif($phase == 2) {
    $phase = 1;
  }

  return $phase;
}

sub get_slice {
  my $config = shift;
  my $chr = shift;

  if(!defined($config->{slice_cache}->{$chr})) {
    $config->{slice_cache}->{$chr} = Bio::EnsEMBL::Slice->new(
      -COORD_SYSTEM      => $config->{coord_system},
      -START             => 1,
      -END               => $config->{fasta_db}->length($chr),
      -SEQ_REGION_NAME   => $chr,
      -SEQ_REGION_LENGTH => $config->{fasta_db}->length($chr)
    );
  }

  return $config->{slice_cache}->{$chr};
}

sub get_regions {
  my $config = shift;
  my $obj = shift;

  my ($obj_start, $obj_end) = ($obj->start, $obj->end);

  my @regions;
  my ($s, $e);

  # initial start, end
  $s = $config->{cache_region_size} * int($obj->{start} / $config->{cache_region_size});
  $e = $s + $config->{cache_region_size};
  $s += 1;
  push @regions, $s.'-'.$e;

  while($e < $obj_end) {
    $s += $config->{cache_region_size};
    $e += $config->{cache_region_size};
    push @regions, $s.'-'.$e;
  }

  return \@regions;
}

sub get_dump_file_name {
  my $config = shift;
  my $chr    = shift;
  my $region = shift;
  my $type   = shift;

  $type ||= 'transcript';

  if($type eq 'transcript') {
    $type = '';
  }
  else {
    $type = '_'.$type;
  }

  my $dir = $config->{dir}.'/'.$chr;
  my $dump_file = $dir.'/'.$region.$type.(defined($config->{tabix}) ? '_tabix' : '').'.gz';

  # make directory if it doesn't exist
  if(!(-e $dir) && defined($config->{write_cache})) {
    make_path($dir);
  }

  return $dump_file;
}

sub export_data {
  my $config = shift;
  my $chr = shift;

  my $hash;
  foreach my $tr(values %{$config->{transcripts}}) {
    next unless $tr->seq_region_name eq $chr;

    fix_transcript($tr);
    next unless check_transcript($config, $tr);

    foreach my $region(@{get_regions($config, $tr)}) {
      push @{$hash->{$region}}, $tr;
    }
  }

  foreach my $region(keys %{$hash}) {

    # unique and sort list
    my @array =
      sort {$a->start <=> $b->start || $a->end <=> $b->end}
      values %{{
        map {$_->stable_id => $_}
        grep {defined($_)}
        @{$hash->{$region}}
      }};

    dump_transcript_cache($config, {$chr => \@array}, $chr, $region);

    # remove all the dumped transcripts
    delete $hash->{$region};
  }

  delete $config->{transcripts};
  delete $config->{transcripts_by_gff_id};
  delete $config->{genes};
}

sub fix_transcript {
  my $tr = shift;

  # no CDS processed, can't be protein coding
  if(!defined($tr->{translation}->{start_exon})) {
    delete $tr->{translation};
    $tr->{biotype} ||= 'pseudogene';
  }
}

sub check_transcript {
  my $config = shift;
  my $tr = shift;

  my @errors;

  push @errors, "Object is not a transcript" unless $tr->isa('Bio::EnsEMBL::Transcript');
  if(!(defined($tr->{_trans_exon_array}) && scalar @{$tr->get_all_Exons})) {
    push @errors, "No exons found";
  }
  else {
    push @errors, "Exon missing phase" if grep {not defined $_->phase} @{$tr->get_all_Exons};
  }

  if($tr->biotype eq 'protein_coding') {
    push @errors, "No start_exon defined on translation" unless defined $tr->{translation}->{start_exon};
    push @errors, "No end_exon defined on translation" unless defined $tr->{translation}->{end_exon};
  }

  if(scalar @errors) {
    warn "WARNING: Transcript ".$tr->stable_id." fails checks: ".(join ", ", @errors)."\n";
    return 0;
  }

  return 1;
}

# dumps out transcript cache to file
sub dump_transcript_cache {
  my $config = shift;
  my $tr_cache = shift;
  my $chr = shift;
  my $region = shift;

  #debug("Dumping cached transcript data") unless defined($config->{quiet});

  my $dump_file = get_dump_file_name($config, $chr, $region, 'transcript');

  #debug("Writing to $dump_file") unless defined($config->{quiet});
  #$DB::single = 1;

  # storable
  open my $fh, "| gzip -9 -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
  nstore_fd($tr_cache, $fh);
  close $fh;
}

sub get_seq_region_synonyms {
  my $config = shift;

  if(!exists($config->{seq_region_synonyms})) {
    my $reg = 'Bio::EnsEMBL::Registry';

    $reg->load_registry_from_db(
      -host       => $config->{host},
      -user       => $config->{user},
      -pass       => $config->{password},
      -port       => $config->{port},
      -db_version => $config->{db_version},
      -species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
    );

    my $sa;

    eval {$sa = $reg->get_adaptor($config->{species}, 'core', 'slice')};

    my %synonyms = ();

    if(!$@ && $sa) {
      my $slices = $sa->fetch_all('toplevel');

      foreach my $slice(@$slices) {
        my $slice_name = $slice->seq_region_name;
        $synonyms{$_->name} = $slice_name for @{$slice->get_all_synonyms};
      }
    }

    $config->{seq_region_synonyms} = \%synonyms;
  }

  return $config->{seq_region_synonyms};
}

sub read_synonyms_file {
  my $config = shift;
  my $file = shift || $config->{synonyms};

  die("ERROR: Could not read from synonyms file $file\n") unless $file && -e $file && open IN, $file;

  while(<IN>) {
    chomp;

    my ($name, $synonym) = split(/\s+/, $_);

    # store both ways around
    $config->{seq_region_synonyms}->{$name} = $synonym;
    $config->{seq_region_synonyms}->{$synonym} = $name;
  }

  close IN;
}

# setup FASTA file
sub setup_fasta() {
  my $config = shift;

  die "ERROR: Specified FASTA file/directory not found" unless -e $config->{fasta};

  my $index_type = 'faidx';

  eval q{ use Bio::DB::HTS::Faidx; };

  if($@) {

    # if FASTA file is gzipped, we can't index it without Faidx
    if($config->{fasta} =~ /\.gz$/) {

      # but first check that the unpacked file doesn't exist
      my $unpacked_fa = $config->{fasta};
      $unpacked_fa =~ s/\.gz$//;

      # if it does, we can use it instead
      if(-e $unpacked_fa) {
        $config->{fasta} = $unpacked_fa;
      }
      else {
        die("ERROR: Cannot index gzipped FASTA file without Faidx\n");
      }
    }


    debug("Unable to use Faidx, falling back to Bio::DB::Fasta\n") unless defined($config->{quiet});
    $index_type = 'bioperl';

    # try and fall back to
    eval q{ use Bio::DB::Fasta; };

    if($@) {
      die("ERROR: Could not load required Faidx or BioPerl module\n");
    }
  }

  if($index_type eq 'faidx') {

    # try to overwrite sequence method in Slice
    eval q{
      package Bio::EnsEMBL::Slice;

      # define a global variable so that we can pull in config hash
      our $config;

      {
        # don't want a redefine warning spat out, thanks
        no warnings 'redefine';

        # overwrite seq method to read from FASTA DB
        sub seq {
          my $self = shift;

          # special case for in-between (insert) coordinates
          return '' if($self->start() == $self->end() + 1);

          my $seq ;
          my $length = 0 ;
          if(defined($config->{fasta_db})) {
            my $location_string = $self->seq_region_name.":".$self->start."-".$self->end ;
            ($seq, $length) = $config->{fasta_db}->get_sequence($location_string) ;
            reverse_comp(\$seq) if $self->strand < 0;
          }

          else {
            return $self->{'seq'} if($self->{'seq'});

            if($self->adaptor()) {
              my $seqAdaptor = $self->adaptor()->db()->get_SequenceAdaptor();
              return ${$seqAdaptor->fetch_by_Slice_start_end_strand($self,1,undef,1)};
            }
          }

          # default to a string of Ns if we couldn't get sequence
          $seq ||= 'N' x $self->length();

          return $seq;
        }
      }

      1;
    };
  }
  else {
    eval q{
      package Bio::EnsEMBL::Slice;

      # define a global variable so that we can pull in config hash
      our $config;

      {
        # don't want a redefine warning spat out, thanks
        no warnings 'redefine';

        # overwrite seq method to read from FASTA DB
        sub seq {
          my $self = shift;

          # special case for in-between (insert) coordinates
          return '' if($self->start() == $self->end() + 1);

          my $seq;

          if(defined($config->{fasta_db})) {
            $seq = $config->{fasta_db}->seq($self->seq_region_name, $self->start => $self->end);
            reverse_comp(\$seq) if $self->strand < 0;
          }

          else {
            return $self->{'seq'} if($self->{'seq'});

            if($self->adaptor()) {
              my $seqAdaptor = $self->adaptor()->db()->get_SequenceAdaptor();
              return ${$seqAdaptor->fetch_by_Slice_start_end_strand($self,1,undef,1)};
            }
          }

          # default to a string of Ns if we couldn't get sequence
          $seq ||= 'N' x $self->length();

          return $seq;
        }
      }

      1;
    };
  }

  if($@) {
    debug($@) unless defined($config->{quiet});
    die("ERROR: Could not redefine sequence method\n");
  }

  # copy to Slice for offline sequence fetching
  {
    no warnings 'once';
    $Bio::EnsEMBL::Slice::config = $config;
  }

  # spoof a coordinate system
  $config->{coord_system} = Bio::EnsEMBL::CoordSystem->new(
    -NAME => 'chromosome',
    -RANK => 1,
  );

  # check lock file
  my $lock_file = $config->{fasta};
  $lock_file .= -d $config->{fasta} ? '/.vep.lock' : '.vep.lock';

  # lock file exists, indexing failed
  if(-e $lock_file) {
    for(qw(.fai .index .gzi /directory.index /directory.fai .vep.lock)) {
      unlink($config->{fasta}.$_) if -e $config->{fasta}.$_;
    }
  }

  my $index_exists = 0;

  for my $fn(map {$config->{fasta}.$_} qw(.fai .index .gzi /directory.index /directory.fai)) {
    if(-e $fn) {
      $index_exists = 1;
      last;
    }
  }

  # create lock file
  unless($index_exists) {
    debug("Creating FASTA index") unless defined($config->{quiet});

    open LOCK, ">$lock_file" or die("ERROR: Could not write to FASTA lock file $lock_file\n");
    print LOCK "1\n";
    close LOCK;
  }

  # run indexing
  $config->{fasta_db} = $index_type eq 'faidx' ? Bio::DB::HTS::Faidx->new($config->{fasta}) : Bio::DB::Fasta->new($config->{fasta});

  # remove lock file
  unlink($lock_file) unless $index_exists;
}


# DEBUG AND STATUS METHODS
##########################

# gets time
sub get_time() {
    my @time = localtime(time());

    # increment the month (Jan = 0)
    $time[4]++;

    # add leading zeroes as required
    for my $i(0..4) {
        $time[$i] = "0".$time[$i] if $time[$i] < 10;
    }

    # put the components together in a string
    my $time =
        ($time[5] + 1900)."-".
        $time[4]."-".
        $time[3]." ".
        $time[2].":".
        $time[1].":".
        $time[0];

    return $time;
}

# prints debug output with time
sub debug {
    my $text = (@_ ? (join "", @_) : "No message");
    my $time = get_time;

    print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}



sub usage {
    my $usage =<<END;
#--------------------------------#
# GTF/GFF to VEP cache converter #
#--------------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

This script creates a VEP cache from a GTF or GFF file containing transcript/exon
definitions and a FASTA file containing the reference sequence for the same
species.

Usage: perl gtf2vep.pl [arguments]

Options
=======

-h | --help               Display this message and quit
-i | --input [file]       GTF/GFF file (may be gzipped)
-f | --fasta [file]       FASTA file or directory containing FASTA files
-s | --species [species]  Species name
-d | --db_version [n]     Database version - must match version of API in use
--dir [dir]               Root directory for cache (default = '\$HOME/.vep/')

--host [host]             Database connection settings; used to look up synonyms
--port [port]             for chromosome names. Connects by default to
--user [user]             anonymous\@ensembldb.ensembl.org:3306
--password [pass]

--verbose                 Show more warnings and status output
--no_transcripts          Enable this flag for GTF files with no explicit
                          transcript or gene lines
--source_biotype          Older Ensembl GTFs have the transcript biotype in the
                          GTF source field - enable this flag for such files

END

    print $usage;
}
