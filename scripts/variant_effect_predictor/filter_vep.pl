#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

filter_vep.pl - a script to filter results from the Variant Effect Predictor

http://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;

# set output autoflush for progress bars
$| = 1;

# configure from command line opts
my $config = configure(scalar @ARGV);

# run the main sub routine
main($config);

sub configure {
  my $args = shift;
  
  my $config = {};
  
  GetOptions(
    $config,
    'help|h',                  # displays help message
    
    'test=i',                  # test run on n lines
    'count|c',                 # only print a count
    'list|l',                  # list available columns
    
    'input_file|i=s',          # input file
    'output_file|o=s',         # output file
    'force_overwrite',         # force overwrite output
    'format=s',                # input format
    'gz',                      # force read as gzipped
    'only_matched',            # rewrite CSQ field in VCF with only matched "blobs"
    
    'ontology|y',              # use ontology for matching consequence terms
    'host=s',                  # DB options
    'user=s',
    'pass=s',
    'port=i',
    'version=i',
    'registry=s',
    
    'start|s=i',               # skip first N results
    'limit|l=i',               # return max N results
    
    'filter|f=s@',             # filter
  ) or die "ERROR: Failed to parse command-line flags\n";
  
  # print usage message if requested or no args supplied
  if(defined($config->{help}) || !$args) {
    &usage;
    exit(0);
  }
  
  # set defaults
  $config->{start} ||= 0;
  $config->{limit} ||= 1e12;
  
  # default to stdout
  $config->{output_file} ||= 'stdout';
  
  die("ERROR: No valid filters given\n") if !defined($config->{filter}) && !defined($config->{list});
  
  return $config;
}

sub main {
  my $config = shift;
  
  # input
  my $in_fh = new FileHandle;
    
  if(defined($config->{input_file})) {
    # check defined input file exists
    die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
    
    if($config->{input_file} =~ /\.gz$/ || defined($config->{gz})){
      $in_fh->open("gzip -dc ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
    }
    else {
      $in_fh->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
    }
  }
  else {
    $in_fh = 'STDIN';
  }
  
  # output
  my $out_fh = new FileHandle;
  
  if(-e $config->{output_file} && !defined($config->{force_overwrite})) {
    die("ERROR: Output file ", $config->{output_file}, " already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
  }
  elsif($config->{output_file} =~ /stdout/i) {
    $out_fh = *STDOUT;
  }
  else {
    $out_fh->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ".$config->{output_file}."\n");
  }
  
  my (@raw_headers, @headers, $count, $line_number);
  
  while(<$in_fh>) {
    chomp;
    
    $config->{line_number}++;
    
    # header line?
    if(/^\#/) {
      push @raw_headers, $_;
    }
    else {
      $line_number++;
      last if defined($config->{test}) && $line_number > $config->{test};
      
      # parse headers before processing input
      if(!(scalar @headers)) {
        die("ERROR: No headers found in input file\n") unless scalar @raw_headers;
        
        # write headers to output
        unless(defined($config->{count}) || defined($config->{list})) {
          print $out_fh join("\n", @raw_headers);
          print $out_fh "\n";
        }
        
        # parse into data structures
        parse_headers($config, \@raw_headers);
        @headers = @{$config->{headers} || $config->{col_headers}};
        $config->{allowed_fields}->{$_} = 1 for (@{$config->{headers}}, @{$config->{col_headers}});
        
        if(defined($config->{list})) {
          print "Available fields:\n\n";
          print "$_\n" for sort keys %{$config->{allowed_fields}};
          exit(0);
        }
        
        # parse filters
        $config->{filters} = parse_filters($config, $config->{filter});
      }
      
      my $line = $_;
      
      # get format
      $config->{format} ||= detect_format($line);
      die("ERROR: Could not detect input file format - perhaps you need to specify it with --format?\n") unless defined($config->{format});
      die("ERROR: --only_matched is compatible only with VCF files\n") if $config->{format} ne 'vcf' && defined($config->{only_matched});
      
      my (@data, @chunks);
      
      if($config->{format} eq 'vep') {
        push @data, parse_line($line, \@headers, "\t");
        push @chunks, $line;
      }
      elsif($config->{format} eq 'vcf') {
        
        # get main data
        my $main_data = parse_line($line, $config->{col_headers}, "\t");
        
        # get info fields
        foreach my $info_field(split /\;/, (split /\s+/, $line)[7]) {
          my ($field, $value) = split /\=/, $info_field;
          $main_data->{$field} = $value;
        }
        
        # get CSQ stuff
        if($line =~ m/CSQ\=(.+?)(\;|$|\s)/) {
          @chunks = split('\,', $1);
          push @data,
            map {merge_hashes($_, $main_data)}
            map {parse_line($_, \@headers, '\|')}
            @chunks;
        }
      }
      else {
        die("ERROR: Unable to parse data in format ".$config->{format}."\n");
      }
      
      my ($line_pass, @new_chunks);
      
      # test each chunk
      foreach my $i(0..$#data) {
        my $chunk_pass;
        my $parsed_chunk = $data[$i];
        my $raw_chunk = $chunks[$i];
        
        $chunk_pass = run_filters($config, $parsed_chunk, $config->{filters});
        $line_pass += $chunk_pass;
        
        push @new_chunks, $raw_chunk if $chunk_pass;
      }
      
      # update CSQ if using only_matched
      if(defined($config->{only_matched}) && scalar @new_chunks != scalar @chunks) {
        my $new_csq = join(",", @new_chunks);
        $line =~ s/CSQ\=(.+?)(\;|\s|$)/CSQ\=$new_csq$2/;
      }
      
      $count++ if $line_pass;
      
      next unless $count >= $config->{start};
      
      print $out_fh "$line\n" if $line_pass && !defined($config->{count});
      
      last if $count >= $config->{limit} + $config->{start} - 1;
    }
  }
  
  print $out_fh "$count\n" if defined($config->{count});
}

sub parse_headers {
  my $config = shift;
  my $raw_headers = shift;
  
  foreach my $raw_header(@$raw_headers) {
    
    # remove and count hash characters
    my $hash_count = $raw_header =~ s/\#//g;
    
    # field definition (VCF)
    if($hash_count == 2) {
      if($raw_header =~ /INFO\=\<ID\=CSQ/) {
        $raw_header =~ m/Format\: (.+?)\"/;
        $config->{headers} = [split '\|', $1];
      }
      elsif($raw_header =~ /INFO\=\<ID\=(.+?)\,/) {
        $config->{allowed_fields}->{$1} = 1;
      }
      elsif($raw_header =~ m/ (.+?) \:/) {
        $config->{allowed_fields}->{$1} = 1;
      }
    }
    
    # column headers
    else {
      $config->{col_headers} = [split "\t", $raw_header];
    }
  }
}

sub parse_line {
  my $line = shift;
  my $headers = shift;
  my $delimiter = shift;
  
  chomp $line;
  my @split = split($delimiter, $line);
  
  my %data = map {$headers->[$_] => $split[$_]} 0..$#split;
  
  if(defined($data{Extra})) {
    foreach my $chunk(split /\;/, $data{Extra}) {
      my ($key, $value) = split /\=/, $chunk;
      $data{$key} = $value;
    }
  }
  
  map {$data{$_} = undef if $data{$_} eq ''} keys %data;
  
  return \%data;
}

sub parse_filters {
  my $config = shift;
  my $filters = shift;
  
  my @return;
  
  my %filter_synonyms = (
    '>' => 'gt',
    '>=' => 'gte',
    '<' => 'lt',
    '<=' => 'lte',
    
    'is' => 'eq',
    '=' => 'eq',
    
    '!=' => 'ne',
    
    'exists' => 'ex',
    'defined' => 'ex',
    
    'match' => 're',
    'matches' => 're',
    'regex' => 're',
  );
  
  foreach my $filter_list(@$filters) {
    
    while($filter_list =~ m/(.+?)(\sand\s|\sor\s|$)/g) {
      my $filter_string = $1;
      my @words = split /\s+/, $filter_string;
      
      # logic
      my $logic = lc($2 || 'and');
      $logic =~ s/\s+//g;
      
      # invert?
      my $invert;
      if($words[0] =~ /^not$/i) {
        $invert = 1;
        shift @words;
      }
      else {
        $invert = 0;
      }
      
      # field
      my $field = shift @words;
      
      # operator?
      my $operator = scalar @words ? shift @words : 'ex';
      
      # value?
      my $value = scalar @words ? shift @words : undef;
      
      if(!defined($value)) {
        $operator = 'nex' if $operator eq 'ne';
        $operator = 'ex' if $operator eq 'is';
      }
      
      # sub
      my $sub_name = 'filter_'.$operator;
      
      if(!defined(&$sub_name)) {
        if(defined($filter_synonyms{$operator})) {
          $operator = $filter_synonyms{$operator};
          $sub_name = 'filter_'.$operator;
        }
        else {
          die("ERROR: No such operator \"$operator\"\n") unless defined(&$sub_name);
        }
      }
      
      # match field
      if(!defined($config->{allowed_fields}->{$field})) {
        my @matched_fields = grep {$_ =~ /^$field/i} keys %{$config->{allowed_fields}};
        
        if(scalar @matched_fields == 1) {
          $field = $matched_fields[0];
        }
        #else {
        #  die("ERROR: No field matching $field found\n");
        #}
      }
      
      if(defined($config->{ontology}) && $field eq 'Consequence' && $operator eq 'eq') {
        $operator = 'is_child';
        $sub_name = 'filter_is_child';
      }
      
      push @return, {
        predicate => \&$sub_name,
        field     => $field,
        value     => $value,
        invert    => $invert,
        logic     => $logic,
      };
    }
  }
  
  return \@return;
}

sub run_filters {
  my $config = shift;
  my $data = shift;
  my $filters = shift;
  
  my $pass;
  my $logic = 'and';
  
  foreach my $filter(@$filters) {
    
    my $predicate = $filter->{predicate};
    
    # process input
    my $input = $data->{$filter->{field}};
    my $value = $filter->{value};
    
    if(defined($input) && $input =~ /([\w\.]+)?\:?\(?([\-\d\.]*)\)?/ && $filter->{field} ne 'CELL_TYPE') {
      my ($text, $num) = ($1, $2);
      
      if($value =~ /^[\-\d\.]+$/) {
        $input = $text =~ /^[\-\d\.]+$/ ? $text : $num;
      }
      else {
        $input = $text;
      }
    }
    
    # run filter
    my $this_pass;
    
    if(defined($value) && !defined($input)) {
      $this_pass = 0;
    }
    else {
      $this_pass = &$predicate($input, $value, $config);
    }
    
    # invert
    $this_pass = 1 - $this_pass if $filter->{invert};
    
    if(defined($pass)) {
      if($logic eq 'and') {
        $pass *= $this_pass;
      }
      elsif($logic eq 'or') {
        $pass += $this_pass;
      }
    }
    else {
      $pass = $this_pass;
    }
    
    $logic = $filter->{logic};
  }
  
  return $pass;
}

sub merge_hashes {
  my ($x, $y, $add) = @_;

  foreach my $k (keys %$y) {
    if (!defined($x->{$k})) {
      $x->{$k} = $y->{$k};
    } else {
      if(ref($x->{$k}) eq 'ARRAY') {
        $x->{$k} = merge_arrays($x->{$k}, $y->{$k});
      }
      elsif(ref($x->{$k}) eq 'HASH') {
        $x->{$k} = merge_hashes($x->{$k}, $y->{$k}, $add);
      }
      else {
        $x->{$k} = ($add ? $x->{$k} + $y->{$k} : $y->{$k});
      }
    }
  }
  return $x;
}

sub merge_arrays {
  my ($x, $y) = @_;
  
  my %tmp = map {$_ => 1} (@$x, @$y);
  
  return [keys %tmp];
}


# sub-routine to detect format of input
sub detect_format {
  my $line = shift;
  my @data = split /\s+/, $line;
  
  # VCF: 20  14370  rs6054257  G  A  29  0  NS=58;DP=258;AF=0.786;DB;H2  GT:GQ:DP:HQ
  if (
    $data[0] =~ /(chr)?\w+/ &&
    $data[1] =~ /^\d+$/ &&
    $data[3] =~ /^[ACGTN-]+$/i &&
    $data[4] =~ /^([\.ACGTN-]+\,?)+$/i
  ) {
    return 'vcf';
  }
  
  # vep output: ID  1:142849179   -     -     -     -     INTERGENIC
  elsif (
    $data[0] =~ /\w+/ &&
    $data[1] =~ /^\w+?\:\d+(\-\d+)*$/ &&
    scalar @data == 14
  ) {
    return 'vep';
  }
  
  else {
    die("ERROR: Could not detect input file format\n");
  }
}

# basic filters
sub filter_eq  { return $_[0] eq $_[1] }
sub filter_ne  { return $_[0] ne $_[1] }
sub filter_gt  { return $_[0] >  $_[1] }
sub filter_lt  { return $_[0] <  $_[1] }
sub filter_gte { return $_[0] >= $_[1] }
sub filter_lte { return $_[0] <= $_[1] }
sub filter_ex  { return defined($_[0]) }
sub filter_nex { return !defined($_[0]) }

# string
sub filter_re  { return $_[0] =~ /$_[1]/i }
sub filter_nre { return $_[0] !~ /$_[1]/i }

# checks if a value exists in a file or list
sub filter_in {
  my ($data, $list, $config) = @_;
  
  my %compare;
  
  if($list =~ /.+\,.+/) {
    %compare = map {$_ => 1} split(',', $list);
  }
  
  # file?
  elsif(-e $list) {
    if(defined($config->{files}->{$list})) {
      %compare = %{$config->{files}->{$list}};
    }
    else {
      open IN, $list or die("ERROR: Could not read from file $list\n");
      while(<IN>) {
        chomp;
        $compare{$_} = 1;
      }
      close IN;
      $config->{files}->{$list} = \%compare;
    }
  }
  
  else {
    die("ERROR: Could not find/parse list $list\n");
  }
  
  return defined($compare{$data});
}

# uses ontology to see if term is a child term of given parent
sub filter_is_child {
  my ($child, $parent, $config) = @_;
  
  # exact match, don't need to use ontology
  return 1 if filter_re($child, $parent);
  
  # get parent term and descendants and cache it on $config
  if(!defined($config->{descendants}) || !defined($config->{descendants}->{$parent})) {
    
    # connect to DBs here
    if(!defined($config->{ontology_adaptor})) {
      eval {
        use Bio::EnsEMBL::Registry;
      };
      
      if($@) {
        die("ERROR: Could not load Ensembl API modules\n");
      }
      
      $config->{reg} = 'Bio::EnsEMBL::Registry';
      
      # registry file for local DBs?
      if(defined($config->{registry})) {
        $config->{reg}->load_all($config->{registry});
      }
      
      # otherwise manually connect to DB server
      else {
        $config->{host} ||= 'ensembldb.ensembl.org';
        $config->{port} ||= '3306';
        $config->{user} ||= 'anonymous';
        
        $config->{reg}->load_registry_from_db(
          -host       => $config->{host},
          -user       => $config->{user},
          -pass       => $config->{password},
          -port       => $config->{port},
          -db_version => $config->{version},
        );
      }
      
      # get ontology adaptor
      my $oa = $config->{reg}->get_adaptor('Multi','Ontology','OntologyTerm');
      die("ERROR: Could not fetch OntologyTerm adaptor\n") unless defined($oa);
      $config->{ontology_adaptor} = $oa;
    }
    
    my $terms = $config->{ontology_adaptor}->fetch_all_by_name($parent, 'SO');
    die("ERROR: No matching SO terms found for $parent\n") unless $terms && scalar @$terms;
    die("ERROR: Found more than one SO term matching $parent: ".join(", ", map {$_->name} @$terms)."\n") if scalar @$terms > 1;
    my $parent_term = $terms->[0];
    $config->{descendants}->{$parent} = $parent_term->descendants;
  }  
  
  return grep {$_->name =~ /^$child$/i} @{$config->{descendants}->{$parent}};
}

sub usage {
  print qq{#---------------#
# filter_vep.pl #
#---------------#

By Will McLaren (wm2\@ebi.ac.uk)

http://www.ensembl.org/info/docs/variation/vep/vep_script.html#filter

Usage:
perl filter_vep.pl [arguments]
  
--help               -h   Print usage message and exit

--input_file [file]  -i   Specify the input file (i.e. the VEP results file).
                          If no input file is specified, the script will
                          attempt to read from STDIN. Input may be gzipped - to
                          force the script to read a file as gzipped, use --gz
--format [format]         Specify input file format (vep or vcf)

--output_file [file] -o   Specify the output file to write to. If no output file
                          is specified, the script will write to STDOUT
--force_overwrite         Force the script to overwrite the output file if it
                          already exists

--filter [filters]   -f   Add filter. Multiple --filter flags may be used, and
                          are treated as logical ANDs, i.e. all filters must
                          pass for a line to be printed

--list               -l   List allowed fields from the input file
--count              -c   Print only a count of matched lines

--only_matched            In VCF files, the CSQ field that contains the
                          consequence data will often contain more than one
                          "block" of consequence data, where each block
                          corresponds to a variant/feature overlap. Using
                          --only_matched will remove blocks that do not pass the
                          filters. By default, the script prints out the entire
                          VCF line if any of the blocks pass the filters.
                          
--ontology           -y   Use Sequence Ontology to match consequence terms. Use
                          with operator "is" to match against all child terms of
                          your value.
                          e.g. "Consequence is coding_sequence_variant" will
                          match missense_variant, synonymous_variant etc.
                          Requires database connection; defaults to connecting
                          to ensembldb.ensembl.org. Use --host, --port, --user,
                          --password, --version as per
                          variant_effect_predictor.pl to change connection
                          parameters.
};
}
