#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Variant Effect Predictor - a script to predict the consequences of genomic variants

http://www.ensembl.org/info/docs/variation/vep/vep_script.html

Version 71

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use CGI qw/:standard/;
use FindBin qw($Bin);
use lib $Bin;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Variation::Utils::VEP qw(
    parse_line
    vf_to_consequences
    validate_vf
    convert_to_vcf
    load_dumped_adaptor_cache
    dump_adaptor_cache
    get_all_consequences
    get_slice
    build_full_cache
    read_cache_info
    get_time
    debug
    @OUTPUT_COLS
    @REG_FEAT_TYPES
    %FILTER_SHORTCUTS
);

# global vars
my $VERSION = '71';

 
# define headers that would normally go in the extra field
# keyed on the config parameter used to turn it on
my %extra_headers = (
    protein         => ['ENSP'],
    canonical       => ['CANONICAL'],
    ccds            => ['CCDS'],
    hgvs            => ['HGVSc','HGVSp'],
    hgnc            => ['HGNC'],
    sift            => ['SIFT'],
    polyphen        => ['PolyPhen'],
    numbers         => ['EXON','INTRON'],
    domains         => ['DOMAINS'],
    regulatory      => ['MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE'],
    cell_type       => ['CELL_TYPE'],
    individual      => ['IND','ZYG'],
    xref_refseq     => ['RefSeq'],
    check_svs       => ['SV'],
    check_frequency => ['FREQS'],
    gmaf            => ['GMAF'],
    maf_1kg         => ['AFR_MAF','AMR_MAF','ASN_MAF','EUR_MAF'],
    user            => ['DISTANCE'],
    check_existing  => ['CLIN_SIG'],
    biotype         => ['BIOTYPE'],
);

my %extra_descs = (
    'CANONICAL'    => 'Indicates if transcript is canonical for this gene',
    'CCDS'         => 'Indicates if transcript is a CCDS transcript',
    'HGNC'         => 'HGNC gene identifier',
    'ENSP'         => 'Ensembl protein identifer',
    'HGVSc'        => 'HGVS coding sequence name',
    'HGVSp'        => 'HGVS protein sequence name',
    'SIFT'         => 'SIFT prediction',
    'PolyPhen'     => 'PolyPhen prediction',
    'EXON'         => 'Exon number(s) / total',
    'INTRON'       => 'Intron number(s) / total',
    'DOMAINS'      => 'The source and identifer of any overlapping protein domains',
    'MOTIF_NAME'   => 'The source and identifier of a transcription factor binding profile (TFBP) aligned at this position',
    'MOTIF_POS'    => 'The relative position of the variation in the aligned TFBP',
    'HIGH_INF_POS' => 'A flag indicating if the variant falls in a high information position of the TFBP',
    'MOTIF_SCORE_CHANGE' => 'The difference in motif score of the reference and variant sequences for the TFBP',
    'CELL_TYPE'    => 'List of cell types and classifications for regulatory feature',
    'IND'          => 'Individual name',
    'ZYG'          => 'Zygosity of individual genotype at this locus',
    'SV'           => 'IDs of overlapping structural variants',
    'FREQS'        => 'Frequencies of overlapping variants used in filtering',
    'GMAF'         => 'Minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined population',
    'AFR_MAF'      => 'Minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined African population',
    'AMR_MAF'      => 'Minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined American population',
    'ASN_MAF'      => 'Minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined Asian population',
    'EUR_MAF'      => 'Minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined European population',
    'DISTANCE'     => 'Shortest distance from variant to transcript',
    'CLIN_SIG'     => 'Clinical significance of variant from dbSNP',
    'BIOTYPE'      => 'Biotype of transcript',
);

my %ts_tv = (
  'A/G' => 'Ts',
  'G/A' => 'Ts',
  'C/T' => 'Ts',
  'T/C' => 'Ts',
  'A/C' => 'Tv',
  'C/A' => 'Tv',
  'G/T' => 'Tv',
  'T/G' => 'Tv',
  'C/G' => 'Tv',
  'G/C' => 'Tv',
  'A/T' => 'Tv',
  'T/A' => 'Tv',
);

my %colour_keys = (
  'polyphen' => {
    'unknown' => 'blue',
    'benign' => 'green',
    'possibly damaging' => 'orange',
    'probably damaging' => 'red',
  },
  'sift' => {
    'tolerated' => 'green',
    'deleterious' => 'red',
  },
  
  # copied from COLOUR.ini in web code via browser to check colours
  'consequences' => {
    'intergenic_variant'                => 'gray',
    'intron_variant'                    => '#02599c',
    'upstream_gene_variant'             => '#a2b5cd',
    'downstream_gene_variant'           => '#a2b5cd',
    '5_prime_utr_variant'               => '#7ac5cd',
    '3_prime_utr_variant'               => '#7ac5cd',
    'splice_region_variant'             => '#ff7f50',
    'splice_donor_variant'              => '#ff7f50',
    'splice_acceptor_variant'           => '#ff7f50',
    'frameshift_variant'                => '#ff69b4',
    'transcript_ablation'               => '#ff0000',
    'transcript_amplification'          => '#ff69b4',
    'inframe_insertion'                 => '#ff69b4',
    'inframe_deletion'                  => '#ff69b4',
    'synonymous_variant'                => '#76ee00',
    'stop_retained_variant'             => '#76ee00',
    'missense_variant'                  => '#ffd700',
    'initiator_codon_variant'           => '#ffd700',
    'stop_gained'                       => '#ff0000',
    'stop_lost'                         => '#ff0000',
    'mature_mirna_variant'              => '#458b00',
    'non_coding_exon_variant'           => '#32cd32',
    'nc_transcript_variant'             => '#32cd32',
    'incomplete_terminal_codon_variant' => '#ff00ff',
    'nmd_transcript_variant'            => '#ff4500',
    'coding_sequence_variant'           => '#458b00',
    'tfbs_ablation'                     => 'brown',
    'tfbs_amplification'                => 'brown',
    'tf_binding_site_variant'           => 'brown',
    'regulatory_region_variant'         => 'brown',
    'regulatory_region_ablation'        => 'brown',
    'regulatory_region_amplification'   => 'brown',
  },
);

# set output autoflush for progress bars
$| = 1;

# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine
&main($config);

# this is the main sub-routine - it needs the configured $config hash
sub main {
    my $config = shift;
    
    debug("Starting...") unless defined $config->{quiet};
    
    # this is for counting seconds
    $config->{start_time} = time();
    $config->{last_time} = time();
    
    # this is for stats
    $config->{stats}->{start_time} = get_time();
    
    my $tr_cache = {};
    my $rf_cache = {};
    
    # create a hash to hold slices so we don't get the same one twice
    my %slice_cache = ();
    
    my @vfs;    
    my ($vf_count, $total_vf_count);
    my $in_file_handle = $config->{in_file_handle};
    
    # initialize line number in config
    $config->{line_number} = 0;
    
    # read the file
    while(<$in_file_handle>) {
        chomp;
        
        $config->{line_number}++;
        
        # header line?
        if(/^\#/) {
            
            # retain header lines if we are outputting VCF
            if(defined($config->{vcf})) {
                push @{$config->{headers}}, $_;
            }
            
            # line with sample labels in VCF
            if(defined($config->{individual}) && /^#CHROM/) {
                my @split = split /\s+/;
                
                # no individuals
                die("ERROR: No individual data found in VCF\n") if scalar @split <= 9;
                
                # get individual column indices
                my %ind_cols = map {$split[$_] => $_} (9..$#split);
                
                # all?
                if(scalar @{$config->{individual}} == 1 && $config->{individual}->[0] =~ /^all$/i) {
                    $config->{ind_cols} = \%ind_cols;
                }
                else {
                    my %new_ind_cols;
                    
                    # check we have specified individual(s)
                    foreach my $ind(@{$config->{individual}}) {
                        die("ERROR: Individual named \"$ind\" not found in VCF\n") unless defined $ind_cols{$ind};
                        $new_ind_cols{$ind} = $ind_cols{$ind};
                    }
                    
                    $config->{ind_cols} = \%new_ind_cols;
                }
            }
            
            next;
        }
        
        # configure output file
        $config->{out_file_handle} ||= &get_out_file_handle($config);
        
        # some lines (pileup) may actually parse out into more than one variant
        foreach my $vf(@{&parse_line($config, $_)}) {
            
            $vf->{_line} = $_ ;#if defined($config->{vcf}) || defined($config->{original});
            
            # now get the slice
            if(!defined($vf->{slice})) {
                my $slice;
                
                # don't get slices if we're using cache
                # we can steal them from transcript objects later
                if((!defined($config->{cache}) && !defined($config->{whole_genome})) || defined($config->{check_ref}) || defined($config->{convert})) {
                    
                    # check if we have fetched this slice already
                    if(defined $slice_cache{$vf->{chr}}) {
                        $slice = $slice_cache{$vf->{chr}};
                    }
                    
                    # if not create a new one
                    else {
                        
                        $slice = &get_slice($config, $vf->{chr});
                        
                        # if failed, warn and skip this line
                        if(!defined($slice)) {
                            warn("WARNING: Could not fetch slice named ".$vf->{chr}." on line ".$config->{line_number}."\n") unless defined $config->{quiet};
                            next;
                        }    
                        
                        # store the hash
                        $slice_cache{$vf->{chr}} = $slice;
                    }
                }
                
                $vf->{slice} = $slice;
            }
            
            # validate the VF
            next unless validate_vf($config, $vf);
            
            # make a name if one doesn't exist
            $vf->{variation_name} ||= $vf->{chr}.'_'.$vf->{start}.'_'.($vf->{allele_string} || $vf->{class_SO_term});
            
            # jump out to convert here
            if(defined($config->{convert})) {
                &convert_vf($config, $vf);
                next;
            }
            
            if(defined $config->{whole_genome}) {
                push @vfs, $vf;
                $vf_count++;
                $total_vf_count++;
                
                if($vf_count == $config->{buffer_size}) {
                    debug("Read $vf_count variants into buffer") unless defined($config->{quiet});
                    
                    print_line($config, $_) foreach @{get_all_consequences($config, \@vfs)};
                    
                    # calculate stats
                    my $total_rate = sprintf("%.0f vars/sec", $total_vf_count / ((time() - $config->{start_time}) || 1));
                    my $rate = sprintf("%.0f vars/sec", $vf_count / ((time() - $config->{last_time}) || 1));
                    $config->{last_time} = time();
                    
                    debug("Processed $total_vf_count total variants ($rate, $total_rate total)") unless defined($config->{quiet});
                    
                    @vfs = ();
                    $vf_count = 0;
                }
            }
            else {
                print_line($config, $_) foreach @{vf_to_consequences($config, $vf)};
                $vf_count++;
                $total_vf_count++;
                debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
            }
        }
    }
    
    # if in whole-genome mode, finish off the rest of the buffer
    if(defined $config->{whole_genome} && scalar @vfs) {
        debug("Read $vf_count variants into buffer") unless defined($config->{quiet});
        
        print_line($config, $_) foreach @{get_all_consequences($config, \@vfs)};
        
        # calculate stats
        my $total_rate = sprintf("%.0f vars/sec", $total_vf_count / ((time() - $config->{start_time}) || 1));
        my $rate = sprintf("%.0f vars/sec", $vf_count / ((time() - $config->{last_time}) || 1));
        $config->{last_time} = time();
        
        debug("Processed $total_vf_count total variants ($rate, $total_rate total)") unless defined($config->{quiet});
    }
    
    debug($config->{stats}->{filter_count}, "/$total_vf_count variants remain after filtering") if (defined($config->{filter}) || defined($config->{check_frequency})) && !defined($config->{quiet});
    
    debug("Executed ", defined($Bio::EnsEMBL::DBSQL::StatementHandle::count_queries) ? $Bio::EnsEMBL::DBSQL::StatementHandle::count_queries : 'unknown number of', " SQL statements") if defined($config->{count_queries}) && !defined($config->{quiet});
    
    # finalise run-time stats
    $config->{stats}->{var_count} = $total_vf_count;
    $config->{stats}->{end_time} = get_time();
    $config->{stats}->{run_time} = time() - $config->{start_time};
    
    # write stats
    unless(defined($config->{no_stats})) {
      summarise_stats($config);
      debug("Wrote stats summary to ".$config->{stats_file}) unless defined($config->{quiet});
    }
    
    # close HTML output
    if(defined($config->{html}) && defined($config->{html_file_handle})) {
      my $fh = $config->{html_file_handle};
      print $fh "</tbody><tfoot><tr>".$config->{_th}."</tr></tfoot></table><p>&nbsp;</p></div></html></body>\n</html>\n";
      $fh->close;
    }
    
    debug("Finished!") unless defined $config->{quiet};
}

# sets up configuration hash that is used throughout the script
sub configure {
    my $args = shift;
    
    my $config = {};
    
    my @ARGV_copy = @ARGV;
    $config->{stats}->{options} = \@ARGV_copy;
    
    GetOptions(
        $config,
        'help',                    # displays help message
        
        # input options,
        'config=s',                # config file name
        'input_file|i=s',          # input file name
        'format=s',                # input file format
        
        # DB options
        'species=s',               # species e.g. human, homo_sapiens
        'registry=s',              # registry file
        'host=s',                  # database host
        'port=s',                  # database port
        'user=s',                  # database user name
        'password=s',              # database password
        'db_version=i',            # Ensembl database version to use e.g. 62
        'genomes',                 # automatically sets DB params for e!Genomes
        'refseq',                  # use otherfeatures RefSeq DB instead of Ensembl
        
        # runtime options
        'most_severe',             # only return most severe consequence
        'summary',                 # only return one line per variation with all consquence types
        'per_gene',                # only return most severe per gene
        'buffer_size=i',           # number of variations to read in before analysis
        'chunk_size=s',            # size in bases of "chunks" used in internal hash structure
        'failed=i',                # include failed variations when finding existing
        'no_whole_genome',         # disables now default whole-genome mode
        'whole_genome',            # proxy for whole genome mode - now just warns user
        'gp',                      # read coords from GP part of INFO column in VCF (probably only relevant to 1KG)
        'chr=s',                   # analyse only these chromosomes, e.g. 1-5,10,MT
        'check_ref',               # check supplied reference allele against DB
        'check_existing',          # find existing co-located variations
        'check_svs',               # find overlapping structural variations
        'check_alleles',           # only attribute co-located if alleles are the same
        'check_frequency',         # enable frequency checking
        'gmaf',                    # add global MAF of existing var
        'maf_1kg',                 # add 1KG MAFs of existing vars
        'freq_filter=s',           # exclude or include
        'freq_freq=f',             # frequency to filter on
        'freq_gt_lt=s',            # gt or lt (greater than or less than)
        'freq_pop=s',              # population to filter on
        'filter_common',           # shortcut to MAF filtering
        'allow_non_variant',       # allow non-variant VCF lines through
        'individual=s',            # give results by genotype for individuals
        'phased',                  # force VCF genotypes to be interpreted as phased
        'fork=i',                  # fork into N processes
        
        # verbosity options
        'verbose|v',               # print out a bit more info while running
        'quiet',                   # print nothing to STDOUT (unless using -o stdout)
        'no_progress',             # don't display progress bars
        
        # output options
        'everything|e',            # switch on EVERYTHING :-)
        'output_file|o=s',         # output file name
        'html',                    # generate an HTML version of output
        'stats_file|sf=s',         # stats file name
        'no_stats',                # don't write stats file
        'force_overwrite',         # force overwrite of output file if already exists
        'terms|t=s',               # consequence terms to use e.g. NCBI, SO
        'coding_only',             # only return results for consequences in coding regions
        'canonical',               # indicates if transcript is canonical
        'ccds',                    # output CCDS identifer
        'xref_refseq',             # output refseq mrna xref
        'protein',                 # add e! protein ID to extra column
        'biotype',                 # add biotype of transcript to output
        'hgnc',                    # add HGNC gene ID to extra column
        'hgvs',                    # add HGVS names to extra column
        'sift=s',                  # SIFT predictions
        'polyphen=s',              # PolyPhen predictions
        'condel=s',                # Condel predictions
        'regulatory',              # enable regulatory stuff
        'cell_type=s' => ($config->{cell_type} ||= []),             # filter cell types for regfeats
        'convert=s',               # convert input to another format (doesn't run VEP)
        'filter=s',                # run in filtering mode
        'no_intergenic',           # don't print out INTERGENIC consequences
        'gvf',                     # produce gvf output
        'vcf',                     # produce vcf output
        'original',                # produce output in input format
        'no_consequences',         # don't calculate consequences
        'lrg',                     # enable LRG-based features
        'fields=s',                # define your own output fields
        'domains',                 # output overlapping protein features
        'numbers',                 # include exon and intron numbers
        
        # cache stuff
        'database',                # must specify this to use DB now
        'cache',                   # use cache
        'write_cache',             # enables writing to the cache
        'build=s',                 # builds cache from DB from scratch; arg is either all (all top-level seqs) or a list of chrs
        'no_adaptor_cache',        # don't write adaptor cache
        'prefetch',                # prefetch exons, translation, introns, codon table etc for each transcript
        'strip',                   # strips adaptors etc from objects before caching them
        'rebuild=s',               # rebuilds cache by reading in existing then redumping - probably don't need to use this any more
        'dir=s',                   # dir where cache is found (defaults to $HOME/.vep/)
        'dir_cache=s',             # specific directory for cache
        'dir_plugins=s',           # specific directory for plugins
        'cache_region_size=i',     # size of region in bases for each cache file
        'no_slice_cache',          # tell API not to cache features on slice
        'standalone',              # standalone renamed offline
        'offline',                 # offline mode uses minimal set of modules installed in same dir, no DB connection
        'skip_db_check',           # don't compare DB parameters with cached
        'compress=s',              # by default we use zcat to decompress; user may want to specify gzcat or "gzip -dc"
        'custom=s' => ($config->{custom} ||= []), # specify custom tabixed bgzipped file with annotation
        'tmpdir=s',                # tmp dir used for BigWig retrieval
        'plugin=s' => ($config->{plugin} ||= []), # specify a method in a module in the plugins directory
        'fasta=s',                 # file or dir containing FASTA files with reference sequence
        'freq_file=s',             # file containing freqs to add to cache build
        
        # debug
        'cluck',                   # these two need some mods to Bio::EnsEMBL::DBSQL::StatementHandle to work. Clucks callback trace and SQL
        'count_queries',           # counts SQL queries executed
        'admin',                   # allows me to build off public hosts
        'debug',                   # print out debug info
        'tabix',                   # experimental use tabix cache files
    ) or die "ERROR: Failed to parse command-line flags\n";
    
    # print usage message if requested or no args supplied
    if(defined($config->{help}) || !$args) {
        &usage;
        exit(0);
    }
    
    # config file?
    if(defined $config->{config}) {
        read_config_from_file($config, $config->{config});
    }
    
    # dir is where the cache and plugins live
    my $default_dir = join '/', ($ENV{'HOME'}, '.vep');
    $config->{dir_plugins} ||= $config->{dir} || $default_dir;
    $config->{dir} ||= $config->{dir_cache} || $default_dir;

    # ini file?
    my $ini_file = $config->{dir}.'/vep.ini';
    
    if(-e $ini_file) {
        read_config_from_file($config, $ini_file);
    }

    # can't be both quiet and verbose
    die "ERROR: Can't be both quiet and verbose!\n" if defined($config->{quiet}) && defined($config->{verbose});
    
    # check forking
    if(defined($config->{fork})) {
        die "ERROR: Fork number must be greater than 1\n" if $config->{fork} <= 1;
        
        # check we can use MIME::Base64
        eval q{ use MIME::Base64; };
        
        if($@) {
            debug("WARNING: Unable to load MIME::Base64, forking disabled") unless defined($config->{quiet});
            delete $config->{fork};
        }
        else {
            
            # try a practice fork
            my $pid = fork;
            
            if(!defined($pid)) {
                debug("WARNING: Fork test failed, forking disabled") unless defined($config->{quiet});
                delete $config->{fork};
            }
            elsif($pid) {
                waitpid($pid, 0);
            }
            elsif($pid == 0) {
                exit(0);
            }
        }
    }
    
    # check file format
    if(defined $config->{format}) {
        die "ERROR: Unrecognised input format specified \"".$config->{format}."\"\n" unless $config->{format} =~ /^(pileup|vcf|guess|hgvs|ensembl|id|vep)$/i;
    }
    
    # check convert format
    if(defined $config->{convert}) {
        die "ERROR: Unrecognised output format for conversion specified \"".$config->{convert}."\"\n" unless $config->{convert} =~ /vcf|ensembl|pileup|hgvs/i;
        
        # disable stats
        $config->{no_stats} = 1;
    }
    
    
    if(defined($config->{fasta})) {
        die "ERROR: Specified FASTA file/directory not found" unless -e $config->{fasta};
        
        eval q{ use Bio::DB::Fasta; };
        
        if($@) {
            die("ERROR: Could not load required BioPerl module\n");
        }
        
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
        
        if($@) {
            die("ERROR: Could not redefine sequence method\n");
        }
        
        # copy to Slice for offline sequence fetching
        $Bio::EnsEMBL::Slice::config = $config;
        
        # spoof a coordinate system
        $config->{coord_system} = Bio::EnsEMBL::CoordSystem->new(
            -NAME => 'chromosome',
            -RANK => 1,
        );
        
        debug("Checking/creating FASTA index") unless defined($config->{quiet});
        $config->{fasta_db} = Bio::DB::Fasta->new($config->{fasta});
    }
    
    # check if user still using --standalone
    if(defined $config->{standalone}) {
        die "ERROR: --standalone replaced by --offline\n";
    }
    
    # connection settings for Ensembl Genomes
    if($config->{genomes}) {
        $config->{host} ||= 'mysql.ebi.ac.uk';
        $config->{port} ||= 4157;
    }
    
    # connection settings for main Ensembl
    else {
        $config->{species} ||= "homo_sapiens";
        $config->{host}    ||= 'ensembldb.ensembl.org';
        $config->{port}    ||= 5306;
    }
    
    # refseq or core?
    if(defined($config->{refseq})) {
        $config->{core_type} = 'otherfeatures';
    }
    else {
        $config->{core_type} = 'core';
    }
    
    # check one of database/cache/offline/build
    if(!grep {defined($config->{$_})} qw(database cache offline build convert)) {
      die qq{
IMPORTANT INFORMATION:

The VEP can read gene data from either a local cache or local/remote databases.

Using a cache is the fastest and most efficient way to use the VEP. The
included INSTALL.pl script can be used to fetch and set up cache files from the
Ensembl FTP server. Simply run "perl INSTALL.pl" and follow the instructions, or
see the documentation pages listed below.

If you have already set up a cache, use "--cache" or "--offline" to use it.

It is possible to use the public databases hosted at ensembldb.ensembl.org, but
this is slower than using the cache and concurrent and/or long running VEP jobs
can put strain on the Ensembl servers, limiting availability to other users.

To enable using databases, add the flag "--database".

Documentation
Installer: http://www.ensembl.org/info/docs/variation/vep/vep_script.html#installer
Cache: http://www.ensembl.org/info/docs/variation/vep/vep_script.html#cache

      }
    };
    
    # output term
    if(defined $config->{terms}) {
        die "ERROR: Unrecognised consequence term type specified \"".$config->{terms}."\" - must be one of ensembl, so, ncbi\n" unless $config->{terms} =~ /ensembl|display|so|ncbi/i;
        if($config->{terms} =~ /ensembl|display/i) {
            $config->{terms} = 'display';
        }
        else {
            $config->{terms} = uc($config->{terms});
        }
    }
    
    # everything?
    if(defined($config->{everything})) {
        my %everything = (
            sift       => 'b',
            polyphen   => 'b',
            ccds       => 1,
            hgvs       => 1,
            hgnc       => 1,
            numbers    => 1,
            domains    => 1,
            regulatory => 1,
            canonical  => 1,
            protein    => 1,
            biotype    => 1,
            gmaf       => 1,
        );
        
        $config->{$_} = $everything{$_} for keys %everything;
        
        # these ones won't work with offline
        delete $config->{hgvs} if defined($config->{offline}) && !defined($config->{fasta_db});
    }
    
    # check nsSNP tools
    foreach my $tool(grep {defined $config->{lc($_)}} qw(SIFT PolyPhen Condel)) {
        die "ERROR: Unrecognised option for $tool \"", $config->{lc($tool)}, "\" - must be one of p (prediction), s (score) or b (both)\n" unless $config->{lc($tool)} =~ /^(s|p|b)/;
        
        #die "ERROR: $tool not available for this species\n" unless $config->{species} =~ /human|homo/i;
        
        die "ERROR: $tool functionality is now available as a VEP Plugin - see http://www.ensembl.org/info/docs/variation/vep/vep_script.html#plugins\n" if $tool eq 'Condel';
    }
    
    # force quiet if outputting to STDOUT
    if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
        delete $config->{verbose} if defined($config->{verbose});
        $config->{quiet} = 1;
    }
    
    # individual(s) specified?
    if(defined($config->{individual})) {
        $config->{individual} = [split /\,/, $config->{individual}];
        
        # force allow_non_variant
        $config->{allow_non_variant} = 1;
    }
    
    # regulatory has to be on for cell_type
    if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
        $config->{regulatory} = 1;
        $config->{cell_type} = [map {split /\,/, $_} @{$config->{cell_type}}];
    }
    
    # summarise options if verbose
    if(defined $config->{verbose}) {
        my $header =<<INTRO;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

Configuration options:

INTRO
        print $header;
        
        my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
        
        foreach my $key(sort keys %$config) {
            next if ref($config->{$key}) eq 'ARRAY' && scalar @{$config->{$key}} == 0;
            print $key.(' ' x (($max_length - length($key)) + 4)).(ref($config->{$key}) eq 'ARRAY' ? join "\t", @{$config->{$key}} : $config->{$key})."\n";
        }
        
        print "\n".("-" x 20)."\n\n";
    }
    
    # check custom annotations
    for my $i(0..$#{$config->{custom}}) {
        my $custom = $config->{custom}->[$i];
        
        my ($filepath, $shortname, $format, $type, $coords) = split /\,/, $custom;
        $type ||= 'exact';
        $format ||= 'bed';
        $coords ||= 0;
        
        # check type
        die "ERROR: Type $type for custom annotation file $filepath is not allowed (must be one of \"exact\", \"overlap\")\n" unless $type =~ /exact|overlap/;
        
        # check format
        die "ERROR: Format $format for custom annotation file $filepath is not allowed (must be one of \"bed\", \"vcf\", \"gtf\", \"gff\", \"bigwig\")\n" unless $format =~ /bed|vcf|gff|gtf|bigwig/;
        
        # bigwig format
        if($format eq 'bigwig') {
            # check for bigWigToWig
            die "ERROR: bigWigToWig does not seem to be in your path - this is required to use bigwig format custom annotations\n" unless `which bigWigToWig 2>&1` =~ /bigWigToWig$/;
        }
        
        else {
            # check for tabix
            die "ERROR: tabix does not seem to be in your path - this is required to use custom annotations\n" unless `which tabix 2>&1` =~ /tabix$/;
            
            # remote files?
            if($filepath =~ /tp\:\/\//) {
                my $remote_test = `tabix $filepath 1:1-1 2>&1`;
                if($remote_test =~ /fail/) {
                    die "$remote_test\nERROR: Could not find file or index file for remote annotation file $filepath\n";
                }
                elsif($remote_test =~ /get_local_version/) {
                    debug("Downloaded tabix index file for remote annotation file $filepath") unless defined($config->{quiet});
                }
            }
        
            # check files exist
            else {
                die "ERROR: Custom annotation file $filepath not found\n" unless -e $filepath;
                die "ERROR: Tabix index file $filepath\.tbi not found - perhaps you need to create it first?\n" unless -e $filepath.'.tbi';
            }
        }
        
        $config->{custom}->[$i] = {
            'file'   => $filepath,
            'name'   => $shortname || 'CUSTOM'.($i + 1),
            'type'   => $type,
            'format' => $format,
            'coords' => $coords,
        };
    }
    
    # check if using filter and original
    die "ERROR: You must also provide output filters using --filter to use --original\n" if defined($config->{original}) && !defined($config->{filter});
    
    # filter by consequence?
    if(defined($config->{filter})) {
        
        my %filters = map {$_ => 1} split /\,/, $config->{filter};
        
        # add in shortcuts
        foreach my $filter(keys %filters) {
            my $value = 1;
            if($filter =~ /^no_/) {
                delete $filters{$filter};
                $filter =~ s/^no_//g;
                $value = 0;
                $filters{$filter} = $value;
            }
            
            if(defined($FILTER_SHORTCUTS{$filter})) {
                delete $filters{$filter};
                $filters{$_} = $value for keys %{$FILTER_SHORTCUTS{$filter}};
            }
        }
        
        $config->{filter} = \%filters;
        
        $config->{stats}->{filter_count} = 0;
    }
    
    # set defaults
    $config->{user}              ||= 'anonymous';
    $config->{buffer_size}       ||= 5000;
    $config->{chunk_size}        ||= '50kb';
    $config->{output_file}       ||= "variant_effect_output.txt";
    $config->{stats_file}        ||= $config->{output_file}."_summary.html";
    $config->{tmpdir}            ||= '/tmp';
    $config->{format}            ||= 'guess';
    $config->{terms}             ||= 'SO';
    $config->{cache_region_size} ||= 1000000;
    $config->{compress}          ||= 'zcat';
    
    # can't use a whole bunch of options with most_severe
    if(defined($config->{most_severe})) {
        foreach my $flag(qw(no_intergenic protein hgnc sift polyphen coding_only ccds canonical xref_refseq numbers domains summary)) {
            die "ERROR: --most_severe is not compatible with --$flag\n" if defined($config->{$flag});
        }
    }
    
    # can't use a whole bunch of options with summary
    if(defined($config->{summary})) {
        foreach my $flag(qw(no_intergenic protein hgnc sift polyphen coding_only ccds canonical xref_refseq numbers domains most_severe)) {
            die "ERROR: --summary is not compatible with --$flag\n" if defined($config->{$flag});
        }
    }
    
    # frequency filtering
    if(defined($config->{filter_common})) {
        $config->{check_frequency} = 1;
        
        # set defaults
        $config->{freq_freq}   ||= 0.01;
        $config->{freq_filter} ||= 'exclude';
        $config->{freq_pop}    ||= '1KG_ALL';
        $config->{freq_gt_lt}  ||= 'gt';
    }
    
    if(defined($config->{check_frequency})) {
        foreach my $flag(qw(freq_freq freq_filter freq_pop freq_gt_lt)) {
            die "ERROR: To use --check_frequency you must also specify flag --$flag\n" unless defined $config->{$flag};
        }
        
        # need to set check_existing
        $config->{check_existing} = 1;
    }
    
    $config->{check_existing} = 1 if defined $config->{check_alleles} || defined $config->{gmaf} || defined $config->{maf_1kg};
    
    # warn users still using whole_genome flag
    if(defined($config->{whole_genome})) {
        debug("INFO: Whole-genome mode is now the default run-mode for the script. To disable it, use --no_whole_genome") unless defined($config->{quiet});
    }
    
    $config->{whole_genome}      = 1 unless defined $config->{no_whole_genome};
    $config->{failed}            = 0 unless defined $config->{failed};
    $config->{chunk_size}        =~ s/mb?/000000/i;
    $config->{chunk_size}        =~ s/kb?/000/i;
    $config->{cache_region_size} =~ s/mb?/000000/i;
    $config->{cache_region_size} =~ s/kb?/000/i;
    
    # cluck and display executed SQL?
    $Bio::EnsEMBL::DBSQL::StatementHandle::cluck = 1 if defined($config->{cluck});
    
    # offline needs cache, can't use HGVS
    if(defined($config->{offline})) {
        $config->{cache} = 1;
        
        die("ERROR: Cannot generate HGVS coordinates in offline mode without a FASTA file (see --fasta)\n") if defined($config->{hgvs}) && !defined($config->{fasta_db});
        die("ERROR: Cannot use HGVS as input in offline mode\n") if $config->{format} eq 'hgvs';
        die("ERROR: Cannot use variant identifiers as input in offline mode\n") if $config->{format} eq 'id';
        die("ERROR: Cannot do frequency filtering in offline mode\n") if defined($config->{check_frequency}) && $config->{freq_pop} !~ /1kg.*(all|afr|amr|asn|eur)/i;
        die("ERROR: Cannot retrieve overlapping structural variants in offline mode\n") if defined($config->{check_sv});
        die("ERROR: Cannot check reference sequences without a FASTA file (see --fasta)\n") if defined($config->{check_ref}) && !defined($config->{fasta_db});
    }
    
    # write_cache needs cache
    $config->{cache} = 1 if defined $config->{write_cache};
    
    # no_slice_cache, prefetch and whole_genome have to be on to use cache
    if(defined($config->{cache})) {
        $config->{prefetch} = 1;
        $config->{no_slice_cache} = 1;
        $config->{whole_genome} = 1;
        $config->{strip} = 1;
    }
    
    $config->{build} = $config->{rebuild} if defined($config->{rebuild});
    
    # force options for full build
    if(defined($config->{build})) {
        $config->{prefetch} = 1;
        $config->{hgnc} = 1;
        $config->{no_slice_cache} = 1;
        $config->{cache} = 1;
        $config->{strip} = 1;
        $config->{write_cache} = 1;
        $config->{cell_type} = [1] if defined($config->{regulatory});
    }
    
    # connect to databases
    $config->{reg} = &connect_to_dbs($config);
    
    # complete dir with species name and db_version
    $config->{dir} .= '/'.(
        join '/', (
            defined($config->{offline}) ? $config->{species} : ($config->{reg}->get_alias($config->{species}) || $config->{species}),
            $config->{db_version} || $config->{reg}->software_version
        )
    );
    
    # warn user cache directory doesn't exist
    if(!-e $config->{dir}) {
        
        # if using write_cache
        if(defined($config->{write_cache})) {
            debug("INFO: Cache directory ", $config->{dir}, " not found - it will be created") unless defined($config->{quiet});
        }
        
        # want to read cache, not found
        elsif(defined($config->{cache})) {
            die("ERROR: Cache directory ", $config->{dir}, " not found");
        }
    }
    
    if(defined($config->{cache})) {
        # read cache info
        if(read_cache_info($config)) {
            debug("Read existing cache info") unless defined $config->{quiet};
        }
    }
   
    # we configure plugins here because they can sometimes switch on the 
    # regulatory config option
    configure_plugins($config);
    
    # include regulatory modules if requested
    if(defined($config->{regulatory})) {
        # do the use statements here so that users don't have to have the
        # funcgen API installed to use the rest of the script
        eval q{
            use Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor;
            use Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor;
            use Bio::EnsEMBL::Funcgen::MotifFeature;
            use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
            use Bio::EnsEMBL::Funcgen::BindingMatrix;
        };
        
        if($@) {
            die("ERROR: Ensembl Funcgen API must be installed to use --regulatory or plugins that deal with regulatory features\n");
        }
    }
    
    # user defined custom output fields
    if(defined($config->{fields})) {
        $config->{fields} = [split ',', $config->{fields}];
        debug("Output fields redefined (".scalar @{$config->{fields}}." defined)") unless defined($config->{quiet});
        $config->{fields_redefined} = 1;
    }
    $config->{fields} ||= \@OUTPUT_COLS;
    
    # suppress warnings that the FeatureAdpators spit if using no_slice_cache
    Bio::EnsEMBL::Utils::Exception::verbose(1999) if defined($config->{no_slice_cache});
    
    # get adaptors (don't get them in offline mode)
    unless(defined($config->{offline})) {
        
        if(defined($config->{cache}) && !defined($config->{write_cache})) {
            
            # try and load adaptors from cache
            if(!&load_dumped_adaptor_cache($config)) {
                &get_adaptors($config);
                &dump_adaptor_cache($config) if defined($config->{write_cache}) && !defined($config->{no_adaptor_cache});
            }
            
            # check cached adaptors match DB params
            else {
                my $dbc = $config->{sa}->{dbc};
            
                my $ok = 1;
                
                if($dbc->{_host} ne $config->{host}) {
                    
                    # ens-livemirror, useastdb and ensembldb should all have identical DBs
                    unless(
                        (
                            $dbc->{_host} eq 'ens-livemirror'
                            || $dbc->{_host} eq 'ensembldb.ensembl.org'
                            || $dbc->{_host} eq 'useastdb.ensembl.org'
                        ) && (
                            $config->{host} eq 'ens-livemirror'
                            || $config->{host} eq 'ensembldb.ensembl.org'
                            || $config->{host} eq 'useastdb.ensembl.org'
                        )
                    ) {
                        $ok = 0;
                    }
                    
                    unless(defined($config->{skip_db_check})) {
                        # but we still need to reconnect
                        debug("INFO: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}, " - reconnecting to host") unless defined($config->{quiet});
                        
                        &get_adaptors($config);
                    }
                }
                
                if(!$ok) {
                    if(defined($config->{skip_db_check})) {
                        debug("INFO: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}) unless defined($config->{quiet});
                    }
                    else {
                        die "ERROR: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}, ". If you are sure this is OK, rerun with -skip_db_check flag set";
                    }
                }
            }
        }
        else {
            &get_adaptors($config);
            &dump_adaptor_cache($config) if defined($config->{write_cache}) && !defined($config->{no_adaptor_cache});
        }
        
        # reg adaptors (only fetches if not retrieved from cache already)
        &get_reg_adaptors($config) if defined($config->{regulatory});
    }
    
    # check cell types
    if(defined($config->{cell_type}) && scalar @{$config->{cell_type}} && !defined($config->{build})) {
        my $cls = '';
        
        if(defined($config->{cache})) {
            $cls = $config->{cache_cell_types};
        }
        else {
            my $cta = $config->{RegulatoryFeature_adaptor}->db->get_CellTypeAdaptor();
            $cls = join ",", map {$_->name} @{$cta->fetch_all};
        }
        
        foreach my $cl(@{$config->{cell_type}}) {
            die "ERROR: cell type $cl not recognised; available cell types are:\n$cls\n" unless $cls =~ /(^|,)$cl(,|$)/;
        }
    }
    
    # get terminal width for progress bars
    unless(defined($config->{quiet})) {
        my $width;
        
        # module may not be installed
        eval q{
            use Term::ReadKey;
        };
        
        if(!$@) {
            my ($w, $h);
            
            # module may be installed, but e.g.
            eval {
                ($w, $h) = GetTerminalSize();
            };
            
            $width = $w if defined $w;
        }
        
        $width ||= 60;
        $width -= 12;
        $config->{terminal_width} = $width;
    }
    
    # jump out to build cache if requested
    if(defined($config->{build})) {
        
        if($config->{host} =~ /^(ensembl|useast)db\.ensembl\.org$/ && !defined($config->{admin})) {
            die("ERROR: Cannot build cache using public database server ", $config->{host}, "\n");
        }
        
        # get 1KG freqs
        if(defined($config->{'freq_file'})) {
            my ($freq_file, @file_pops) = split /\,/, $config->{'freq_file'};
            debug("Loading extra frequencies from $freq_file") unless defined($config->{quiet});
            
            open IN, $freq_file or die "ERROR: Could not open frequencies file $freq_file\n";
            while(<IN>) {
                chomp;
                my @data = split /\t/;
                my $id = shift @data;
                
                # sanity check
                die("ERROR: column count in frequency file $freq_file does not match specified populations ".(join ",", @file_pops)."\n") unless scalar @data == scalar @file_pops;
                
                $config->{'freqs'}->{$id} = join(" ", @data);
            }
            close IN;
            
            # add pops to $config
            $config->{'freq_file_pops'} = \@file_pops;
        }
        
        # build the cache
        debug("Building cache for ".$config->{species}) unless defined($config->{quiet});
        build_full_cache($config);
        
        # exit script
        debug("Finished building cache") unless defined($config->{quiet});
        exit(0);
    }
    
    
    # warn user DB will be used for SIFT/PolyPhen/HGVS/frequency/LRG
    if(defined($config->{cache})) {
        
        # these two def depend on DB
        foreach my $param(grep {defined $config->{$_}} qw(hgvs lrg check_sv check_ref)) {
            debug("INFO: Database will be accessed when using --$param") unless defined($config->{quiet}) or ($param =~ /hgvs|check_ref/ and defined($config->{fasta_db}));
        }
        
        debug("INFO: Database will be accessed when using --check_frequency with population ".$config->{freq_pop}) if !defined($config->{quiet}) and defined($config->{check_frequency}) && $config->{freq_pop} !~ /1kg.*(all|afr|amr|asn|eur)/i;
        
        # as does using HGVS or IDs as input
        debug("INFO: Database will be accessed when using --format ", $config->{format}) if ($config->{format} eq 'id' || $config->{format} eq 'hgvs') && !defined($config->{quiet});
        
        # the rest may be in the cache
        foreach my $param(grep {defined $config->{$_}} qw(sift polyphen regulatory)) {
            next if defined($config->{'cache_'.$param});
            debug("INFO: Database will be accessed when using --$param; consider using the complete cache containing $param data (see documentation for details)") unless defined($config->{quiet});
        }
    }
    
    # get list of chrs if supplied
    if(defined($config->{chr})) {
        my %chrs;
        
        foreach my $val(split /\,/, $config->{chr}) {
            my @nnn = split /\-/, $val;
            
            foreach my $chr($nnn[0]..$nnn[-1]) {
                $chrs{$chr} = 1;
            }
        }
        
        $config->{chr} = \%chrs;
    }
    
    # get input file handle
    $config->{in_file_handle} = &get_in_file_handle($config);
    
    return $config;
}

# reads config from a file
sub read_config_from_file {
    my $config = shift;
    my $file = shift;
    
    open CONFIG, $file or die "ERROR: Could not open config file \"$file\"\n";
    
    while(<CONFIG>) {
        next if /^\#/;
        my @split = split /\s+|\=/;
        my $key = shift @split;
        $key =~ s/^\-//g;
        
        if(defined($config->{$key}) && ref($config->{$key}) eq 'ARRAY') {
            push @{$config->{$key}}, @split;
        }
        else {
            $config->{$key} ||= $split[0];
        }
    }
    
    close CONFIG;
    
    # force quiet if outputting to STDOUT
    if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
        delete $config->{verbose} if defined($config->{verbose});
        $config->{quiet} = 1;
    }
    
    debug("Read configuration from $file") unless defined($config->{quiet});
}

# configures custom VEP plugins
sub configure_plugins {

    my $config = shift;
    
    $config->{plugins} = [];
    
    if (my @plugins = @{ $config->{plugin} }) {

        # add the Plugins directory onto @INC

        unshift @INC, $config->{dir_plugins}."/Plugins";

        for my $plugin (@plugins) {

            # parse out the module name and parameters

            my ($module, @params) = split /,/, $plugin;

            # check we can use the module
            
            eval qq{
                use $module;
            };
            if ($@) {
                debug("Failed to compile plugin $module: $@") unless defined($config->{quiet});
                next;
            }
            
            # now check we can instantiate it, passing any parameters to the constructor
            
            my $instance;
            
            eval {
                $instance = $module->new($config, @params);
            };
            if ($@) {
                debug("Failed to instantiate plugin $module: $@") unless defined($config->{quiet});
                next;
            }

            # check that the versions match
            
            my $plugin_version;
            
            if ($instance->can('version')) {
                $plugin_version = $instance->version;
            }
            
            my $version_ok = 1;

            if ($plugin_version) {
                my ($plugin_major, $plugin_minor, $plugin_maintenance) = split /\./, $plugin_version;
                my ($major, $minor, $maintenance) = split /\./, $VERSION;
    
                if ($plugin_major != $major) {
                    debug("Warning: plugin $plugin version ($plugin_version) does not match the current VEP version ($VERSION)") unless defined($config->{quiet});
                    $version_ok = 0;
                }
            }
            else {
                debug("Warning: plugin $plugin does not define a version number") unless defined($config->{quiet});
                $version_ok = 0;
            }

            debug("You may experience unexpected behaviour with this plugin") unless defined($config->{quiet}) || $version_ok;

            # check that it implements all necessary methods
            
            for my $required(qw(run get_header_info check_feature_type check_variant_feature_type)) {
                unless ($instance->can($required)) {
                    debug("Plugin $module doesn't implement a required method '$required', does it inherit from BaseVepPlugin?") unless defined($config->{quiet});
                    next;
                }
            }
           
            # all's good, so save the instance in our list of plugins
            
            push @{ $config->{plugins} }, $instance;
            
            debug("Loaded plugin: $module") unless defined($config->{quiet}); 

            # for convenience, check if the plugin wants regulatory stuff and turn on the config option if so
            
            if (grep { $_ =~ /motif|regulatory/i } @{ $instance->feature_types }) {
                debug("Fetching regulatory features for plugin: $module") unless defined($config->{quiet});
                $config->{regulatory} = 1;
            }
        }
    }
} 

# connects to DBs (not done in offline mode)
sub connect_to_dbs {
    my $config = shift;
    
    # get registry
    my $reg = 'Bio::EnsEMBL::Registry';
    
    unless(defined($config->{offline})) {
        # load DB options from registry file if given
        if(defined($config->{registry})) {
            debug("Loading DB config from registry file ", $config->{registry}) unless defined($config->{quiet});
            $reg->load_all(
                $config->{registry},
                $config->{verbose},
                undef,
                $config->{no_slice_cache}
            );
        }
        
        # otherwise manually connect to DB server
        else {
            $reg->load_registry_from_db(
                -host       => $config->{host},
                -user       => $config->{user},
                -pass       => $config->{password},
                -port       => $config->{port},
                -db_version => $config->{db_version},
                -species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
                -verbose    => $config->{verbose},
                -no_cache   => $config->{no_slice_cache},
            );
        }
        
        eval { $reg->set_reconnect_when_lost() };
        
        if(defined($config->{verbose})) {
            # get a meta container adaptors to check version
            my $core_mca = $reg->get_adaptor($config->{species}, 'core', 'metacontainer');
            my $var_mca = $reg->get_adaptor($config->{species}, 'variation', 'metacontainer');
            
            if($core_mca && $var_mca) {
                debug(
                    "Connected to core version ", $core_mca->get_schema_version, " database ",
                    "and variation version ", $var_mca->get_schema_version, " database"
                );
            }
        }
    }
    
    return $reg;
}

# get adaptors from DB
sub get_adaptors {
    my $config = shift;
    
    die "ERROR: No registry" unless defined $config->{reg};
    
    # try fetching a slice adaptor
    eval {
        $config->{sa}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'slice');
    };
    
    if($@) {
        if($@ =~ /not find internal name for species/) {
            my %register = %Bio::EnsEMBL::Registry::registry_register;
            my @species_list = sort keys %{$register{_SPECIES}};
            
            debug("ERROR: Could not find database for species ".$config->{species});
            
            if(scalar @species_list) {
                debug("List of valid species for this server:\n\n".(join "\n", @species_list));
            }
            
            die("\nExiting\n");
        }
        
        else {
            die $@;
        }
    }
    
    # get the remaining core adaptors
    $config->{ga}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'gene');
    $config->{ta}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'transcript');
    $config->{mca} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'metacontainer');
    $config->{csa} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'coordsystem');
    $config->{tra} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'translation');
    
    # get variation adaptors
    $config->{vfa}   = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variationfeature');
    $config->{svfa}  = $config->{reg}->get_adaptor($config->{species}, 'variation', 'structuralvariationfeature');
    $config->{tva}   = $config->{reg}->get_adaptor($config->{species}, 'variation', 'transcriptvariation');
    $config->{pfpma} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'proteinfunctionpredictionmatrix');
    $config->{va}    = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variation');
    
    # get fake ones for species with no var DB
    if(!defined($config->{vfa})) {
        $config->{vfa}  = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($config->{species});
        $config->{svfa} = Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor->new_fake($config->{species});
        $config->{tva}  = Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor->new_fake($config->{species});
    }
    
    # cache schema version
    $config->{mca}->get_schema_version if defined $config->{mca};
    
    # check we got slice adaptor - can't continue without a core DB
    die("ERROR: Could not connect to core database\n") unless defined $config->{sa};
}

# gets regulatory adaptors
sub get_reg_adaptors {
    my $config = shift;

    foreach my $type(@REG_FEAT_TYPES) {
        next if defined($config->{$type.'_adaptor'});
        
        my $adaptor = $config->{reg}->get_adaptor($config->{species}, 'funcgen', $type);
        if(defined($adaptor)) {
            $config->{$type.'_adaptor'} = $adaptor;
        }
        else {
            delete $config->{regulatory};
            last;
        }
    }
}

# gets file handle for input
sub get_in_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;
    
    if(defined($config->{input_file})) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
        
        if($config->{input_file} =~ /\.gz$/){
            $in_file_handle->open($config->{compress}." ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
        else {
            $in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
    }
    
    # no file specified - try to read data off command line
    else {
        $in_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
    }
    
    return $in_file_handle;
}

# gets file handle for output and adds header
sub get_out_file_handle {
    my $config = shift;
    
    # define filehandle to write to
    my $out_file_handle = new FileHandle;
    
    # check if file exists
    if(-e $config->{output_file} && !defined($config->{force_overwrite})) {
        die("ERROR: Output file ", $config->{output_file}, " already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
    }
    
    # do same for stats file
    if(-e $config->{stats_file} && !defined($config->{force_overwrite})) {
        die("ERROR: Stats file ", $config->{stats_file}, " already exists. Specify a different output file with --stats_file or overwrite existing file with --force_overwrite\n");
    }
    
    if($config->{output_file} =~ /stdout/i) {
        $out_file_handle = *STDOUT;
    }
    else {
        $out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
    }
    
    # get stats file handle
    die("ERROR: Stats file name ", $config->{stats_file}, " doesn't end in \".htm\" or \".html\" - some browsers may not be able to open this file\n") unless $config->{stats_file} =~ /htm(l)?$/;
    my $stats_file_handle = new FileHandle;
    $stats_file_handle->open(">".$config->{stats_file}) or die("ERROR: Could not write to stats file ", $config->{stats_file}, "\n");
    $config->{stats_file_handle} = $stats_file_handle;
    
    # HTML output?
    my $html_file_handle;
    
    if(defined($config->{html})) {
      if(-e $config->{output_file}.'.html' && !defined($config->{force_overwrite})) {
          die("ERROR: Stats file ", $config->{stats_file}, " already exists. Specify a different output file with --stats_file or overwrite existing file with --force_overwrite\n");
      }
      
      $html_file_handle = new FileHandle;
      $html_file_handle->open(">".$config->{output_file}.'.html') or die("ERROR: Could not write to HTML file ", $config->{output_file}, ".html\n");
      $config->{html_file_handle} = $html_file_handle;
      
      print $html_file_handle html_head();
    }
    
    # define headers for a VCF file
    my @vcf_headers = (
        '#CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO'
    );
    
    # file conversion, don't want to add normal headers
    if(defined($config->{convert})) {
        # header for VCF
        if($config->{convert} =~ /vcf/i) {
            print $out_file_handle "##fileformat=VCFv4.0\n";
            print $out_file_handle join "\t", @vcf_headers;
            print $out_file_handle "\n";
        }
        
        return $out_file_handle;
    }
    
    # GVF output, no header
    elsif(defined($config->{gvf}) || defined($config->{original})) {
        if(defined($config->{headers}) && defined($config->{original})) {
            print $out_file_handle join "\n", @{$config->{headers}};
            print $html_file_handle join("\n", @{$config->{headers}})."\n</pre>" if defined($config->{html});
        }
        return $out_file_handle;
    }
    
    elsif(defined($config->{vcf})) {
        
        # create an info string for the VCF header        
        my @new_headers;
        
        # if the user has defined the fields themselves, we don't need to worry
        if(defined $config->{fields_redefined}) {
            @new_headers = @{$config->{fields}};
        }
        else {
            @new_headers = (
                
                # get default headers, minus variation name and location (already encoded in VCF)
                grep {
                    $_ ne 'Uploaded_variation' and
                    $_ ne 'Location' and
                    $_ ne 'Extra'
                } @{$config->{fields}},
                
                # get extra headers
                map {@{$extra_headers{$_}}}
                grep {defined $config->{$_}}
                keys %extra_headers
            );
            
            # plugin headers
            foreach my $plugin_header(split /\n/, get_plugin_headers($config)) {
                $plugin_header =~ /\#\# (.+?)\t\:.+/;
                push @new_headers, $1;
            }
            
            # redefine the main headers list in config
            $config->{fields} = \@new_headers;
        }
        
        # add the newly defined headers as a header to the VCF
        my $string = join '|', @{$config->{fields}};
        my @vcf_info_strings = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: '.$string.'">');
        
        # add custom headers
        foreach my $custom(@{$config->{custom}}) {
            push @vcf_info_strings, '##INFO=<ID='.$custom->{name}.',Number=.,Type=String,Description="'.$custom->{file}.' ('.$custom->{type}.')">';
        }
        
        # if this is already a VCF file, we need to add our new headers in the right place
        if(defined($config->{headers})) {
            
            for my $i(0..$#{$config->{headers}}) {
                if($config->{headers}->[$i] =~ /^\#CHROM\s+POS\s+ID/) {
                    splice(@{$config->{headers}}, $i, 0, @vcf_info_strings);
                }
            }
            
            print $out_file_handle join "\n", @{$config->{headers}};
            print $out_file_handle "\n";
            
            if(defined($config->{html})) {
                my @tmp = @{$config->{headers}};
                my @cols = split /\s+/, pop @tmp;
                print $html_file_handle join "\n", @tmp;
                print $html_file_handle "\n";
                print $html_file_handle html_table_headers($config, \@cols);
            }
        }
        
        else {
            print $out_file_handle "##fileformat=VCFv4.0\n";
            print $out_file_handle join "\n", @vcf_info_strings;
            print $out_file_handle "\n";
            print $out_file_handle join "\t", @vcf_headers;
            print $out_file_handle "\n";
            
            if(defined($config->{html})) {
                print $html_file_handle "##fileformat=VCFv4.0\n";
                print $html_file_handle join "\n", @vcf_info_strings;
                print $html_file_handle "\n";
                print $html_file_handle html_table_headers($config, \@vcf_headers);
            }
        }
        
        return $out_file_handle;
    }
    
    # make header
    my $time = &get_time;
    my $db_string = $config->{mca}->dbc->dbname." on ".$config->{mca}->dbc->host if defined $config->{mca};
    $db_string .= "\n## Using cache in ".$config->{dir} if defined($config->{cache});
    my $version_string =
        "Using API version ".$config->{reg}->software_version.
        ", DB version ".(defined $config->{mca} && $config->{mca}->get_schema_version ? $config->{mca}->get_schema_version : '?');
        
    # sift/polyphen versions
    foreach my $tool(qw(sift polyphen)) {
        if(defined($config->{$tool})) {
            my $string = 'config_'.$tool.'_version';
            
            if(!defined($config->{$string}) && !defined($config->{offline})) {
                my $var_mca = $config->{reg}->get_adaptor($config->{species}, 'variation', 'metacontainer');
                my $values = $var_mca->list_value_by_key($tool.'_version') if defined($var_mca);
                $config->{$string} = $values->[0] if scalar @$values;
            }
            $version_string .= "\n## $tool version ".$config->{$string} if defined($config->{$string});
        }
    }
    
    # add key for extra column headers based on config
    my $extra_column_keys = join "\n",
        map {'## '.$_.' : '.$extra_descs{$_}}
        sort map {@{$extra_headers{$_}}}
        grep {defined $config->{$_}}
        keys %extra_headers;
    
    my $header =<<HEAD;
## ENSEMBL VARIANT EFFECT PREDICTOR v$VERSION
## Output produced at $time
## Connected to $db_string
## $version_string
## Extra column keys:
$extra_column_keys
HEAD
   
    $header .= get_plugin_headers($config);
    
    # add headers
    print $out_file_handle $header;
    print $html_file_handle $header if defined($config->{html});
    
    # add custom data defs
    if(defined($config->{custom})) {
        foreach my $custom(@{$config->{custom}}) {
            print $out_file_handle '## '.$custom->{name}."\t: ".$custom->{file}.' ('.$custom->{type}.")\n";
            print $html_file_handle '## '.$custom->{name}."\t: ".$custom->{file}.' ('.$custom->{type}.")\n" if defined($config->{html});
        }
    }
    
    # add column headers
    print $out_file_handle '#', (join "\t", @{$config->{fields}});
    print $out_file_handle "\n";
    
    if(defined($config->{html})) {
        print $html_file_handle html_table_headers($config, $config->{fields});
    }
    
    return $out_file_handle;
}

sub get_plugin_headers {

    my $config = shift;

    my $header = "";

    for my $plugin (@{ $config->{plugins} }) {
        if (my $hdr = $plugin->get_header_info) {
            for my $key (keys %$hdr) {
                my $val = $hdr->{$key};
                
                $header .= "## $key\t: $val\n";
            }
        }
    }

    return $header;
}

# convert a variation feature to a line of output
sub convert_vf {
    my $config = shift;
    my $vf = shift;
    
    my $convert_method = 'convert_to_'.lc($config->{convert});
    my $method_ref   = \&$convert_method; 
    
    my $line = &$method_ref($config, $vf);
    my $handle = $config->{out_file_handle};
    
    if(scalar @$line) {
        print $handle join "\t", @$line;
        print $handle "\n";
    }
}

# converts to Ensembl format
sub convert_to_ensembl {
    my $config = shift;
    my $vf = shift;
    
    return [
        $vf->{chr} || $vf->seq_region_name,
        $vf->start,
        $vf->end,
        $vf->allele_string,
        $vf->strand,
        $vf->variation_name
    ];
}


# converts to pileup format
sub convert_to_pileup {
    my $config = shift;
    my $vf = shift;
    
    # look for imbalance in the allele string
    my %allele_lengths;
    my @alleles = split /\//, $vf->allele_string;
    
    foreach my $allele(@alleles) {
        $allele =~ s/\-//g;
        $allele_lengths{length($allele)} = 1;
    }
    
    # in/del
    if(scalar keys %allele_lengths > 1) {
        
        if($vf->allele_string =~ /\-/) {
            
            # insertion?
            if($alleles[0] eq '-') {
                shift @alleles;
            
                for my $i(0..$#alleles) {
                    $alleles[$i] =~ s/\-//g;
                    $alleles[$i] = '+'.$alleles[$i];
                }
            }
            
            else {
                @alleles = grep {$_ ne '-'} @alleles;
                
                for my $i(0..$#alleles) {
                    $alleles[$i] =~ s/\-//g;
                    $alleles[$i] = '-'.$alleles[$i];
                }
            }
            
            @alleles = grep {$_ ne '-' && $_ ne '+'} @alleles;
            
            return [
                $vf->{chr} || $vf->seq_region_name,
                $vf->start - 1,
                '*',
                (join "/", @alleles),
            ];
        }
        
        else {
            warn "WARNING: Unable to convert variant to pileup format on line number ", $config->{line_number} unless defined($config->{quiet});
            return [];
        }
        
    }
    
    # balanced sub
    else {
        return [
            $vf->{chr} || $vf->seq_region_name,
            $vf->start,
            shift @alleles,
            (join "/", @alleles),
        ];
    }
}

# converts to HGVS (hackily returns many lines)
sub convert_to_hgvs {
    my $config = shift;
    my $vf = shift;
    
    # ensure we have a slice
    $vf->{slice} ||= get_slice($config, $vf->{chr});
    
    my $tvs = $vf->get_all_TranscriptVariations;
    
    my @return = values %{$vf->get_all_hgvs_notations()};
    
    if(defined($tvs)) {
        push @return, map {values %{$vf->get_all_hgvs_notations($_->transcript, 'c')}} @$tvs;
        push @return, map {values %{$vf->get_all_hgvs_notations($_->transcript, 'p')}} @$tvs;
    }
    
    return [join "\n", @return];
}

# prints a line of output from the hash
sub print_line {
    my $config = shift;
    my $line = shift;
    return unless defined($line);
    
    my $output;
    my $html_fh = $config->{html_file_handle};
    
    # normal
    if(ref($line) eq 'HASH') {
        my %extra = %{$line->{Extra}};
        
        $line->{Extra} = join ';', map { $_.'='.$line->{Extra}->{$_} } keys %{ $line->{Extra} || {} };
        
        # if the fields have been redefined we need to search through in case
        # any of the defined fields are actually part of the Extra hash
        $output = join "\t", map {
            (defined $line->{$_} ? $line->{$_} : (defined $extra{$_} ? $extra{$_} : '-'))
        } @{$config->{fields}};
        
        if(defined($config->{html})) {
          print $html_fh Tr(
            map {td($_)}
            map {linkify($config, $_)}
            map {
              (defined $line->{$_} ? $line->{$_} : (defined $extra{$_} ? $extra{$_} : '-'))
            }
            @{$config->{fields}}
          );
        }
    }
    
    # gvf/vcf
    else {
        $output = $$line;
    }
    
    my $fh = $config->{out_file_handle};
    print $fh "$output\n";
}

sub summarise_stats {
    my $config = shift;
    
    # convert gene and transcript hashes to counts
    for my $type(qw(genes transcripts regfeats)) {
      $config->{stats}->{$type} = scalar keys %{$config->{stats}->{$type}} if defined $config->{stats}->{$type};
    }
    
    # tot up chromosome counts
    foreach my $chr(keys %{$config->{stats}->{chr}}) {
      $config->{stats}->{chr_totals}->{$chr} += $config->{stats}->{chr}->{$chr}->{$_} for keys %{$config->{stats}->{chr}->{$chr}};
      
      my $start = 0;
      my %tmp;
      
      while($start <= $config->{chr_lengths}->{$chr}) {
        $tmp{$start / 1e6} = $config->{stats}->{chr}->{$chr}->{$start} || 0;
        $start += 1e6;
      }
      
      $config->{stats}->{chr}->{$chr} = \%tmp;
    }
    
    # convert allele changes to Ts/Tv
    map {$config->{stats}->{ts_tv}->{$ts_tv{$_}} += $config->{stats}->{allele_changes}->{$_}} keys %{$config->{stats}->{allele_changes}} if defined($config->{stats}->{allele_changes});
    
    # flesh out protein_pos
    if(defined($config->{stats}->{protein_pos})) {
      if(defined($config->{stats}->{protein_pos}->{10})) {
        $config->{stats}->{protein_pos}->{9} += $config->{stats}->{protein_pos}->{10};
        delete $config->{stats}->{protein_pos}->{10};
      }
      $config->{stats}->{protein_pos}->{$_} ||= 0 for (0..9);
      
      my %tmp = map {$_.'0-'.($_+1).'0%' => $config->{stats}->{protein_pos}->{$_}} keys %{$config->{stats}->{protein_pos}};
      $config->{stats}->{protein_pos} = \%tmp;
    }
    
    # coding cons
    foreach my $con(qw(missense_variant synonymous_variant coding_sequence_variant stop_lost stop_gained frameshift_variant inframe_insertion inframe_deletion)) {
      $config->{stats}->{coding}->{$con} = $config->{stats}->{consequences}->{$con} if defined($config->{stats}->{consequences}->{$con});
    }
    
    # get ranks to sort
    my %cons_ranks = map { $_->{SO_term} => $_->{rank} } values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;
    
    # create pie chart hashes
    my @charts = (
      {
        id => 'var_class',
        title => 'Variant classes',
        header => ['Variant class', 'Count'],
        data => $config->{stats}->{classes},
        type => 'pie',
        sort => 'value',
        height => 200,
      },
      {
        id => 'consequences',
        title => 'Variant consequences',
        header => ['Consequence type', 'Count'],
        data => $config->{stats}->{consequences},
        type => 'pie',
        sort => \%cons_ranks,
        colours => $colour_keys{consequences},
      },
      {
        id => 'coding',
        title => 'Coding consequences',
        header => ['Consequence type', 'Count'],
        data => $config->{stats}->{coding},
        type => 'pie',
        sort => \%cons_ranks,
        colours => $colour_keys{consequences},
      }
    );
    
    foreach my $tool(qw(SIFT PolyPhen)) {
      my $lc_tool = lc($tool);
      
      push @charts, {
        id => $lc_tool,
        title => $tool.' summary',
        header => ['Prediction', 'Count'],
        data => $config->{stats}->{$tool},
        type => 'pie',
        height => 200,
        sort => 'value',
        colours => $colour_keys{$lc_tool},
      } if defined($config->{$lc_tool});
    }
    
    push @charts, {
      id => 'chr',
      title => 'Variants by chromosome',
      header => ['Chromosome','Count'],
      data => $config->{stats}->{chr_totals},
      sort => 'chr',
      type => 'bar',
      options => '{legend: {position: "none"}}',
    };
    
    foreach my $chr(sort {($a !~ /^\d+$/ || $b !~ /^\d+/) ? $a cmp $b : $a <=> $b} keys %{$config->{stats}->{chr}}) {
      push @charts, {
        id => 'chr_'.$chr,
        title => 'Distribution of variants on chromosome '.$chr,
        header => ['Position (mb)', 'Count'],
        data => $config->{stats}->{chr}->{$chr},
        sort => 'chr',
        type => 'line',
        options => '{hAxis: {title: "Position (mb)"}, legend: {position: "none"}}',
        no_table => 1,
        no_link => 1,
      };
    }
    
    push @charts, {
      id => 'protein',
      title => 'Position in protein',
      header => ['Position in protein (percentile)','Count'],
      data => $config->{stats}->{protein_pos},
      sort => 'chr',
      type => 'bar',
      no_table => 1,
      options => '{hAxis: {title: "Position in protein (percentile)", textStyle: {fontSize: 10}}, legend: {position: "none"}}',
    };
    
    my $fh = $config->{stats_file_handle};
    print $fh stats_html_head($config, \@charts);
    
    # create menu
    print $fh div(
      {class => 'sidemenu'},
      div(
        {class => 'sidemenu_head'},
        "Links"
      ),
      div(
        {class => 'sidemenu_body'},
        ul(
          li([
            a({href => '#masthead'}, "Top of page"),
            a({href => '#run_stats'}, "VEP run statistics"),
            a({href => '#gen_stats'}, "General statistics"),
            map {
              a({href => '#'.$_->{id}}, $_->{title})
            } grep { !$_->{no_link} } @charts,
          ])
        ),
      )
    );
    
    print $fh "<div class='main_content'>";
    
    # run stats
    print $fh h3({id => 'run_stats'}, "VEP run statistics");
    
    my @rows = (
      ['VEP version (API)', $VERSION.' ('.$config->{reg}->software_version.')'],
      ['Cache/Database', ($config->{cache} ? $config->{dir} : ($config->{mca} ? $config->{mca}->dbc->dbname." on ".$config->{mca}->dbc->host : '?'))],
      ['Species', $config->{species}],
      ['Command line options', pre(join(" ", @{$config->{stats}->{options}}))],
      ['Start time', $config->{stats}->{start_time}],
      ['End time', $config->{stats}->{end_time}],
      ['Run time', $config->{stats}->{run_time}." seconds"],
      ['Input file (format)', $config->{input_file}.' ('.uc($config->{format}).')'],
      [
        'Output file',
        $config->{output_file}.
        (defined($config->{html}) ? ' '.a({href => $config->{output_file}.'.html'}, '[HTML]') : '').
        ' '.a({href => $config->{output_file}}, '[text]')
      ],
    );
    print $fh table({class => 'stats_table'}, Tr([map {td($_)} @rows]));
    
    # vars in/out stats
    print $fh h3({id => 'gen_stats'}, "General statistics");
    
    @rows = (
      ['Lines of input read', $config->{line_number}],
      ['Variants processed', $config->{stats}->{var_count}],
      ['Variants remaining after filtering', $config->{stats}->{filter_count}],
      [
        'Novel / known variants',
        defined($config->{stats}->{existing}) ?
        sprintf("%s (%.1f\%) / %s (%.1f\%)",
          $config->{stats}->{var_count} - $config->{stats}->{existing},
          100 * (($config->{stats}->{var_count} - $config->{stats}->{existing}) / $config->{stats}->{var_count}),
          $config->{stats}->{existing},
          100 * ($config->{stats}->{existing} / $config->{stats}->{var_count}),
        )
        : '-'
      ],
      ['Overlapped genes', $config->{stats}->{genes}],
      ['Overlapped transcripts', $config->{stats}->{transcripts}],
      ['Overlapped regulatory features', $config->{stats}->{regfeats} || '-'],
    );
    print $fh table({class => 'stats_table'}, Tr([map {td($_)} @rows]));
    
    foreach my $chart(@charts) {
      my $height = $chart->{height} || ($chart->{type} eq 'pie' ? '400' : '200');
      
      print $fh hr();
      print $fh h3({id => $chart->{id}}, $chart->{title});
      print $fh div({id => $chart->{id}."_".$chart->{type}, style => 'width: 800px; height: '.$height.'px'}, '&nbsp;');
      print $fh div({id => $chart->{id}."_table", style => 'width: 800px; height: 200px'}, '&nbsp;') unless $chart->{no_table};
    }
    
    print $fh '</div>';
    print $fh stats_html_tail();
    $config->{stats_file_handle}->close;
}

sub stats_html_head {
    my $config = shift;
    my $charts = shift;
    
    my ($js);
    foreach my $chart(@$charts) {
      my @keys;
      
      # sort data
      if(defined($chart->{sort})) {
        if($chart->{sort} eq 'chr') {
          @keys = sort {($a !~ /^\d+$/ || $b !~ /^\d+/) ? $a cmp $b : $a <=> $b} keys %{$chart->{data}};
        }
        elsif($chart->{sort} eq 'value') {
          @keys = sort {$chart->{data}->{$a} <=> $chart->{data}->{$b}} keys %{$chart->{data}};
        }
        elsif(ref($chart->{sort}) eq 'HASH') {
          @keys = sort {$chart->{sort}->{$a} <=> $chart->{sort}->{$b}} keys %{$chart->{data}};
        }
      }
      else {
        @keys = keys %{$chart->{data}};
      }
      
      my $type = ucfirst($chart->{type});
      
      # add colour
      if(defined($chart->{colours})) {
        my $co = 'slices: ['.join(", ", map { $chart->{colours}->{$_} ? '{color: "'.$chart->{colours}->{$_}.'"}' : '{}' } @keys).']';
        
        if(defined($chart->{options})) {
          $chart->{options} =~ s/}$/, $co}/;
        }
        else {
          $chart->{options} = "{$co}";
        }
      }
      
      # code to draw chart
      $js .= sprintf(
        "var %s = draw$type('%s', '%s', google.visualization.arrayToDataTable([['%s','%s'],%s]), %s);\n",
        $chart->{id}.'_'.$chart->{type},
        $chart->{id}.'_'.$chart->{type},
        $chart->{title},
        $chart->{header}->[0], $chart->{header}->[1],
        join(",", map {"['".$_."',".$chart->{data}->{$_}."]"} @keys),
        $chart->{options} || 'null',
      );
      
      unless($chart->{no_table}) {
        
        # code to draw table
        $js .= sprintf(
          "var %s = drawTable('%s', '%s', google.visualization.arrayToDataTable([['%s','%s'],%s]));\n",
          $chart->{id}.'_table',
          $chart->{id}.'_table',
          $chart->{title},
          $chart->{header}->[0], $chart->{header}->[1],
          join(",", map {"['".$_."',".$chart->{data}->{$_}."]"} @keys)
        );
        
        # interaction between table/chart
        $js .= sprintf(
          qq{
            google.visualization.events.addListener(%s, 'select', function() {
              %s.setSelection(%s.getSelection());
            });
            google.visualization.events.addListener(%s, 'select', function() {
              %s.setSelection(%s.getSelection());
            });
          },
          $chart->{id}.'_'.$chart->{type},
          $chart->{id}.'_table',
          $chart->{id}.'_'.$chart->{type},
          $chart->{id}.'_table',
          $chart->{id}.'_'.$chart->{type},
          $chart->{id}.'_table',
        );
      }
    }
    
    my $html =<<SHTML;
<html>
<head>
  <title>VEP summary</title>
  <script type="text/javascript" src="http://www.google.com/jsapi"></script>
  <script type="text/javascript">
    google.load('visualization', '1', {packages: ['corechart','table']});
  </script>
  <script type="text/javascript">
    
    function init() {
      // charts
      $js
    }
    
    function drawPie(id, title, data, options) {    
      var pie = new google.visualization.PieChart(document.getElementById(id));
      pie.draw(data, options);
      return pie;
    }
    function drawBar(id, title, data, options) {
      var bar = new google.visualization.ColumnChart(document.getElementById(id));
      bar.draw(data, options);
      return bar;
    }
    function drawTable(id, title, data) {
      var table = new google.visualization.Table(document.getElementById(id));
      table.draw(data, null);
      return table;
    }
    function drawLine(id, title, data, options) {
      var line = new google.visualization.LineChart(document.getElementById(id));
      line.draw(data, options);
      return line;
    }
    google.setOnLoadCallback(init);
  </script>
  
  
  <style type="text/css">
    body {
      font-family: arial, sans-serif;
      margin: 0px;
      padding: 0px;
    }
    
    a {color: #36b;}
    a.visited {color: #006;}
    
    .stats_table {
      margin: 5px;
    }
    
    tr:nth-child(odd) {
      background-color: rgb(238, 238, 238);
    }
    
    td {
      padding: 5px;
    }
    
    td:nth-child(odd) {
      font-weight: bold;
    }
    
    h3 {
      color: #666;
    }
    
    .masthead {
      background-color: rgb(51, 51, 102);
      color: rgb(204, 221, 255);
      height: 80px;
      width: 100%;
      padding: 0px;
    }
    
    .main {
      padding: 10px;
    }
    
    .gradient {
      background: #333366; /* Old browsers */
      background: -moz-linear-gradient(left,  #333366 0%, #ffffff 100%); /* FF3.6+ */
      background: -webkit-gradient(linear, left top, right top, color-stop(0%,#333366), color-stop(100%,#ffffff)); /* Chrome,Safari4+ */
      background: -webkit-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Chrome10+,Safari5.1+ */
      background: -o-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Opera 11.10+ */
      background: -ms-linear-gradient(left,  #333366 0%,#ffffff 100%); /* IE10+ */
      background: linear-gradient(to right,  #333366 0%,#ffffff 100%); /* W3C */
      filter: progid:DXImageTransform.Microsoft.gradient( startColorstr='#333366', endColorstr='#ffffff',GradientType=1 ); /* IE6-9 */
      
      padding: 0px;
      height: 80px;
      width: 500px;
      float: right;
      display: inline;
    }
    
    .main_content {
      margin-left: 300px;
    }
    
    .sidemenu {
      width: 260px;
      position: fixed;
      border-style: solid;
      border-width: 2px;
      border-color: rgb(51, 51, 102);
    }
    
    .sidemenu_head {
      width: 250px;
      background-color: rgb(51, 51, 102);
      color: rgb(204, 221, 255);
      padding: 5px;
    }
    
    .sidemenu_body {
      width: 250px;
      padding: 5px;
    }
  </style>
</head>
<body>
<div id="masthead" class="masthead">
  <div style="float: left; display: inline; padding: 10px; height: 80px;">
    <a href="http://www.ensembl.org/"><img src="http://static.ensembl.org/i/e-ensembl.png"></a>
  </div>
  
  <div style="float: right; display: inline; height: 80px; background: white; padding: 10px;">
    <a href="http://www.ensembl.org/info/docs/variation/vep/vep_script.html"><img src="http://www.ensembl.org/img/vep_logo.png"></a>
  </div>
  <div class="gradient">
  </div>
</div>
<div class="main">
SHTML

    return $html;
}

sub stats_html_tail {
  return "\n</div></body>\n</html>\n";
}

sub html_head {
    my $txt_file = $config->{output_file};
    my $stats_file = $config->{stats_file};
    my $html =<<HTML;
<html>
<head>
  <title>VEP output</title>
  <script type="text/javascript" language="javascript" src="http://www.datatables.net/media/javascript/complete.min.js"></script>
  <script class="jsbin" src="http://datatables.net/download/build/jquery.dataTables.nightly.js"></script>
  <script type="text/javascript" language="javascript">
    \$(document).ready(function() {
      \$('#data').dataTable({
        "sPaginationType": "full_numbers"
      });
    });
    
    function fnShowHide( iCol ) {
      /* Get the DataTables object again - this is not a recreation, just a get of the object */
      var oTable = \$('#data').dataTable();
       
      var bVis = oTable.fnSettings().aoColumns[iCol].bVisible;
      oTable.fnSetColumnVis( iCol, bVis ? false : true );
    }
    
    function showAllCols() {
      var oTable = \$('#data').dataTable();
      for (var i=0;i<oTable.fnSettings().aoColumns.length;i++) { 
        oTable.fnSetColumnVis(i, true);
      }
    }
    
    function showHide(lyr) {
      var lyrobj = document.getElementById(lyr);
      
      if(lyrobj.style.height == "0px") {
        lyrobj.style.height = "";
        lyrobj.style.display = "";
      }
      
      else {
        lyrobj.style.height = "0px";
        lyrobj.style.display = "none";
      }
    }
  </script>
  <style type="text/css">
    \@import "http://www.datatables.net/release-datatables/media/css/demo_table.css";
    body {
      font-family: arial, sans-serif;
      margin: 0px;
      padding: 0px;
    }
    
    a {color: #36b;}
    a.visited {color: #006;}
    
    th {
      font-size: 11px;
    }
    td {
      font-size: 11px;
    }
    
    .masthead {
      background-color: rgb(51, 51, 102);
      color: rgb(204, 221, 255);
      height: 80px;
      width: 100%;
      padding: 0px;
    }
    
    .main {
      padding: 10px;
    }
    
    .gradient {
      background: #333366; /* Old browsers */
      background: -moz-linear-gradient(left,  #333366 0%, #ffffff 100%); /* FF3.6+ */
      background: -webkit-gradient(linear, left top, right top, color-stop(0%,#333366), color-stop(100%,#ffffff)); /* Chrome,Safari4+ */
      background: -webkit-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Chrome10+,Safari5.1+ */
      background: -o-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Opera 11.10+ */
      background: -ms-linear-gradient(left,  #333366 0%,#ffffff 100%); /* IE10+ */
      background: linear-gradient(to right,  #333366 0%,#ffffff 100%); /* W3C */
      filter: progid:DXImageTransform.Microsoft.gradient( startColorstr='#333366', endColorstr='#ffffff',GradientType=1 ); /* IE6-9 */
      
      padding: 0px;
      height: 80px;
      width: 500px;
      float: right;
      display: inline;
    }
  </style>
  </head>
  <body>
  <div id="masthead" class="masthead">
    <div style="float: left; display: inline; padding: 10px; height: 80px;">
      <a href="http://www.ensembl.org/"><img src="http://static.ensembl.org/i/e-ensembl.png"></a>
    </div>
    
    <div style="float: right; display: inline; height: 80px; background: white; padding: 10px;">
      <a href="http://www.ensembl.org/info/docs/variation/vep/vep_script.html"><img src="http://www.ensembl.org/img/vep_logo.png"></a>
    </div>
    <div class="gradient">
    </div>
  </div>
  <div class="main">
  <p>
    View: <a href="$stats_file">Summary statistics</a> |
    <a href="$txt_file">as text</a> |
    <a href="javascript:void();" onclick="showHide('header')">Show/hide header</a> |
    <a href="javascript:void();" onclick="showAllCols()">Restore columns</a>
  </p>
  <hr/>
  <pre id="header" style="height:0px; display: none;">
HTML
  return $html;
}

sub html_table_headers {
  my $config = shift;
  my $cols = shift;
  
  my @cols_copy = @$cols;
  
  my $html = qq{</pre><table id="data" class="display"><thead><tr>};
  
  $config->{_th} = join("", map {
    $cols_copy[$_] =~ s/\_/ /g;
    '<th>'.$cols_copy[$_].' '.a({href => 'javascript:void();', onclick => "fnShowHide($_);"}, img({src => 'http://www.ensembl.org/i/16/cross.png', height => 6, width => 6, style => 'border: 1px solid gray; padding: 1px;'})).'</th>'
  } (0..$#cols_copy));
  
  $html .= $config->{_th};
  $html .= qq{</thead></tr><tbody>};
  
  return $html;
}

sub linkify {
  my $config = shift;
  my $string = shift;
  
  my $species = ucfirst($config->{species});
  
  # Ensembl genes
  $string =~ s/(ENS.{0,3}G\d+|CCDS\d+\.?\d+?|N[MP]_\d+\.?\d+?)/a({href => "http:\/\/www.ensembl.org\/$species\/Gene\/Summary\?g=$1", target => "_blank"}, $1)/ge;
  
  # Ensembl transcripts
  $string =~ s/(ENS.{0,3}T\d+)/a({href => "http:\/\/www.ensembl.org\/$species\/Transcript\/Summary\?t=$1", target => "_blank"}, $1)/ge;
  
  # Ensembl regfeats
  $string =~ s/(ENS.{0,3}R\d+)/a({href => "http:\/\/www.ensembl.org\/$species\/Regulation\/Summary\?rf=$1", target => "_blank"}, $1)/ge;
  
  # variant identifiers
  $string =~ s/(rs\d+|COSM\d+|C[DMIX]\d+)/a({href => "http:\/\/www.ensembl.org\/$species\/Variation\/Summary\?v=$1", target => "_blank"}, $1)/gie;
  
  # locations
  while($string =~ m/(^[A-Z\_\d]+?:[1-9]\d+)(\-\d+)?/g) {
    my $loc = $1.($2 ? $2 : '');
    my ($chr, $start, $end) = split /\-|\:/, $loc;
    $end ||= $start;
    
    # adjust +/- 1kb
    $start -= 1000;
    $end   += 1000;
    
    my $link = a({href => "http://www.ensembl.org/$species/Location/View?r=$chr:$start\-$end", target => "_blank"}, $string);
    $string =~ s/$loc/$link/;
  }
  
  # split strings
  $string =~ s/([,;])/$1 /g;
  
  return $string;
}

# outputs usage message
sub usage {
    my $usage =<<END;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

http://www.ensembl.org/info/docs/variation/vep/vep_script.html

Usage:
perl variant_effect_predictor.pl [arguments]

Options
=======

--help                 Display this message and quit
--verbose              Display verbose output as the script runs [default: off]
--quiet                Suppress status and warning messages [default: off]
--no_progress          Suppress progress bars [default: off]

--config               Load configuration from file. Any command line options
                       specified overwrite those in the file [default: off]
                       
--everything           Shortcut switch to turn on commonly used options. See web
                       documentation for details [default: off]
                       
--fork [num_forks]     Use forking to improve script runtime [default: off]

-i | --input_file      Input file - if not specified, reads from STDIN. Files
                       may be gzip compressed.
--format               Specify input file format - one of "ensembl", "pileup",
                       "vcf", "hgvs", "id" or "guess" to try and work out format.
-o | --output_file     Output file. Write to STDOUT by specifying -o STDOUT - this
                       will force --quiet [default: "variant_effect_output.txt"]
--force_overwrite      Force overwriting of output file [default: quit if file
                       exists]
--original             Writes output as it was in input - must be used with --filter
                       since no consequence data is added [default: off]
--vcf                  Write output as VCF [default: off]
--gvf                  Write output as GVF [default: off]
--html                 Write output also as HTML (filename: [output_file].html)
--stats_file           Specify stats summary file [default: [output_file]_summary.html]
--no_stats             Don't write stats summary file
--fields [field list]  Define a custom output format by specifying a comma-separated
                       list of field names. Field names normally present in the
                       "Extra" field may also be specified, including those added by
                       plugin modules. Can also be used to configure VCF output
                       columns [default: off]
                       
--species [species]    Species to use [default: "human"]

-t | --terms           Type of consequence terms to output - one of "SO", "ensembl"
                       [default: SO]
 
--sift=[p|s|b]         Add SIFT [p]rediction, [s]core or [b]oth [default: off]
--polyphen=[p|s|b]     Add PolyPhen [p]rediction, [s]core or [b]oth [default: off]

NB: SIFT predictions are only available for some species, PolyPhen for human only
NB: Condel support has been moved to a VEP plugin module - see documentation

--regulatory           Look for overlaps with regulatory regions. The script can
                       also call if a variant falls in a high information position
                       within a transcription factor binding site. Output lines have
                       a Feature type of RegulatoryFeature or MotifFeature
                       [default: off]
--cell_type [types]    Report only regulatory regions that are found in the given cell
                       type(s). Can be a single cell type or a comma-separated list.
                       The functional type in each cell type is reported under
                       CELL_TYPE in the output. To retrieve a list of cell types, use
                       "--cell_type list" [default: off]
                       
NB: Regulatory consequences are currently available for human and mouse only

--custom [file list]   Add custom annotations from tabix-indexed files. See
                       documentation for full details [default: off]
--plugin [plugin_name] Use named plugin module [default: off]
--hgnc                 Add HGNC gene identifiers to output [default: off]
--hgvs                 Output HGVS identifiers (coding and protein). Requires database
                       connection or --fasta [default: off]
--ccds                 Output CCDS transcript identifiers [default: off]
--xref_refseq          Output aligned RefSeq mRNA identifier for transcript. NB: the
                       RefSeq and Ensembl transcripts aligned in this way MAY NOT, AND
                       FREQUENTLY WILL NOT, match exactly in sequence, exon structure
                       and protein product [default: off]
--protein              Output Ensembl protein identifer [default: off]
--biotype              Output transcript biotype [default: off]
--canonical            Indicate if the transcript for this consequence is the canonical
                       transcript for this gene [default: off]
--domains              Include details of any overlapping protein domains [default: off]
--numbers              Include exon & intron numbers [default: off]

--no_intergenic        Excludes intergenic consequences from the output [default: off]
--coding_only          Only return consequences that fall in the coding region of
                       transcripts [default: off]
--most_severe          Ouptut only the most severe consequence per variation.
                       Transcript-specific columns will be left blank. [default: off]
--summary              Output only a comma-separated list of all consequences per
                       variation. Transcript-specific columns will be left blank.
                       [default: off]
--per_gene             Output only the most severe consequence per gene. Where more
                       than one transcript has the same consequence, the transcript
                       chosen is arbitrary. [default: off]


--check_ref            If specified, checks supplied reference allele against stored
                       entry in Ensembl Core database [default: off]
--check_existing       If specified, checks for existing co-located variations in the
                       Ensembl Variation database [default: off]
--failed [0|1]         Include (1) or exclude (0) variants that have been flagged as
                       failed by Ensembl when checking for existing variants.
                       [default: exclude]
--check_alleles        If specified, the alleles of existing co-located variations
                       are compared to the input; an existing variation will only
                       be reported if no novel allele is in the input (strand is
                       accounted for) [default: off]
--check_svs            Report overlapping structural variants [default: off]

--filter [filters]     Filter output by consequence type. Use this to output only
                       variants that have at least one consequence type matching the
                       filter. Multiple filters can be used separated by ",". By
                       combining this with --original it is possible to run the VEP
                       iteratively to progressively filter a set of variants. See
                       documentation for full details [default: off]

--filter_common        Shortcut flag for the filters below - this will exclude
                       variants that have a co-located existing variant with global
                       MAF > 0.01 (1%). May be modified using any of the following
                       freq_* filters. For human, this can be used in offline mode
                       for the following populations: 1KG_ALL, 1KG_AFR, 1KG_AMR,
                       1KG_ASN, 1KG_EUR
--check_frequency      Turns on frequency filtering. Use this to include or exclude
                       variants based on the frequency of co-located existing
                       variants in the Ensembl Variation database. You must also
                       specify all of the following --freq flags [default: off]
--freq_pop [pop]       Name of the population to use e.g. hapmap_ceu for CEU HapMap,
                       1kg_yri for YRI 1000 genomes. See documentation for more
                       details
--freq_freq [freq]     Frequency to use in filter. Must be a number between 0 and 0.5
--freq_gt_lt [gt|lt]   Specify whether the frequency should be greater than (gt) or
                       less than (lt) --freq_freq
--freq_filter          Specify whether variants that pass the above should be included
  [exclude|include]    or excluded from analysis
--gmaf                 Include global MAF of existing variant from 1000 Genomes
                       Phase 1 in output
--maf_1kg              Include MAF from continental populations (AFR,AMR,ASN,EUR) of
                       1000 Genomes Phase 1 in output
  
--individual [id]      Consider only alternate alleles present in the genotypes of the
                       specified individual(s). May be a single individual, a comma-
                       separated list or "all" to assess all individuals separately.
                       Each individual and variant combination is given on a separate
                       line of output. Only works with VCF files containing individual
                       genotype data; individual IDs are taken from column headers.
--allow_non_variant    Prints out non-variant lines when using VCF input
--phased               Force VCF individual genotypes to be interpreted as phased.
                       For use with plugins that depend on phased state.
                       
--chr [list]           Select a subset of chromosomes to analyse from your file. Any
                       data not on this chromosome in the input will be skipped. The
                       list can be comma separated, with "-" characters representing
                       a range e.g. 1-5,8,15,X [default: off]
--gp                   If specified, tries to read GRCh37 position from GP field in the
                       INFO column of a VCF file. Only applies when VCF is the input
                       format and human is the species [default: off]
                       
--convert              Convert the input file to the output format specified.
  [ensembl|vcf|pileup] Converted output is written to the file specified in
                       --output_file. No consequence calculation is carried out when
                       doing file conversion. [default: off]
                       
--cache                Enables read-only use of cache [default: off]
--dir [directory]      Specify the base cache directory to use [default: "\$HOME/.vep/"]
--dir_cache [dir]      Specify cache directory (if different from --dir)
--dir_plugins [dir]    Specify plugins directory (if different from --dir)
--fasta [file|dir]     Specify a FASTA file or a directory containing FASTA files
                       to use to look up reference sequence. The first time you
                       run the script with this parameter an index will be built
                       which can take a few minutes. This is required if
                       fetching HGVS annotations (--hgvs) or checking reference
                       sequences (--check_ref) in offline mode (--offline), and
                       optional with some performance increase in cache mode
                       (--cache). See documentation for more details

--refseq               Use the otherfeatures database to retrieve transcripts - this
                       database contains RefSeq transcripts (as well as CCDS and
                       Ensembl EST alignments) [default: off]
--database             Enable using databases [default: off]
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]
--genomes              Sets DB connection params for Ensembl Genomes [default: off]
--registry             Registry file to use defines DB connections [default: off]
                       Defining a registry file overrides above connection settings.
--db_version=[number]  Force script to load DBs from a specific Ensembl version. Not
                       advised due to likely incompatibilities between API and DB
                       
--write_cache          Enable writing to cache [default: off]
--build [all|list]     Build a complete cache for the selected species. Build for all
                       chromosomes with --build all, or a list of chromosomes (see
                       --chr). DO NOT USE WHEN CONNECTED TO PUBLIC DB SERVERS AS THIS
                       VIOLATES OUR FAIR USAGE POLICY [default: off]
                       
--compress             Specify utility to decompress cache files - may be "gzcat" or
                       "gzip -dc" Only use if default does not work [default: zcat]
                       
--skip_db_check        ADVANCED! Force the script to use a cache built from a different
                       database than specified with --host. Only use this if you are
                       sure the hosts are compatible (e.g. ensembldb.ensembl.org and
                       useastdb.ensembl.org) [default: off]
--cache_region_size    ADVANCED! The size in base-pairs of the region covered by one
                       file in the cache. [default: 1MB]
                       
--buffer_size          Sets the number of variants sent in each batch [default: 5000]
                       Increasing buffer size can retrieve results more quickly
                       but requires more memory. Only applies to whole genome mode.
--no_whole_genome      Run in old-style, non-whole genome mode [default: off]
END

    print $usage;
}
