#!/usr/bin/perl

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Variant Effect Predictor - a script to predict the consequences of genomic variants

http://www.ensembl.org/info/docs/variation/vep/vep_script.html

Version 2.1

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Storable qw(nstore_fd fd_retrieve);

# we need to manually include all these modules for caching to work
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::TranslationAdaptor;
use Bio::EnsEMBL::DBSQL::TranscriptAdaptor;
use Bio::EnsEMBL::DBSQL::MetaContainer;
use Bio::EnsEMBL::DBSQL::CoordSystemAdaptor;


# debug
#use Time::HiRes qw(tv_interval gettimeofday);

# output columns
my @OUTPUT_COLS = qw(
    Uploaded_variation
    Location
    Allele
    Gene
    Feature
    Feature_type
    Consequence
    cDNA_position
    CDS_position
    Protein_position
    Amino_acids
    Codons
    Existing_variation
    Extra
);

# global vars
my $VERSION = '2.1';

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
    
    my ($include_regions, $transcript_cache);
    
    # scan file if requested
    $include_regions = &scan_file($config) if defined($config->{scan});
    
    # build transcript cache upfront if requested
    $transcript_cache = &cache_transcripts($config, $include_regions) if defined($config->{upfront});
    
    # create a hash to hold slices so we don't get the same one twice
    my %slice_cache = ();
    
    # load slices from the transcript cache if we have it
    # saves us fetching them again
    %slice_cache = %{&build_slice_cache($config, $transcript_cache)} if defined($transcript_cache);
    
    my @new_vfs;
    my %vf_hash;
    
    my $line_number = 0;
    my ($vf_count, $total_vf_count);
    my $in_file_handle = $config->{in_file_handle};
    
    # read the file
    while(<$in_file_handle>) {
      chomp;
      
      $line_number++;
      
      # header line?
      next if /^\#/;
      
      # some lines (pileup) may actually parse out into more than one variant
      foreach my $sub_line(@{&parse_line($config, $_)}) {
        
        # get the sub-line into named variables
        my ($chr, $start, $end, $allele_string, $strand, $var_name) = @{$sub_line};
        
        next if defined($config->{chr}) && !$config->{chr}->{$chr};
        
        # non-variant line from VCF
        next if $chr eq 'non-variant';
        
        # fix inputs
        $chr =~ s/chr//ig unless $chr =~ /^chromosome$/i;
        $chr = 'MT' if $chr eq 'M';
        $strand = ($strand =~ /\-/ ? "-1" : "1");
        $allele_string =~ tr/acgt/ACGT/;
        
        # sanity checks
        unless($start =~ /^\d+$/ && $end =~ /^\d+$/) {
          warn("WARNING: Start $start or end $end coordinate invalid on line $line_number\n") unless defined $config->{quiet};
          next;
        }
        
        unless($allele_string =~ /([ACGT-]+\/*)+/) {
          warn("WARNING: Invalid allele string $allele_string on line $line_number\n") unless defined $config->{quiet};
          next;
        }
        
        # now get the slice
        my $slice;
        
        # don't get slices if we're using cache
        # we can steal them from transcript objects later
        if((!defined($config->{cache}) && !defined($config->{whole_genome})) || defined($config->{check_ref})) {
            
            # check if we have fetched this slice already
            if(defined $slice_cache{$chr}) {
                $slice = $slice_cache{$chr};
            }
            
            # if not create a new one
            else {
                
                $slice = &get_slice($config, $chr);
                
                # if failed, warn and skip this line
                if(!defined($slice)) {
                    warn("WARNING: Could not fetch slice named $chr on line $line_number\n") unless defined $config->{quiet};
                    next;
                }    
                
                # store the hash
                $slice_cache{$chr} = $slice;
            }
        }
        
        # check reference allele if requested
        if(defined $config->{check_ref}) {
            my $ref_allele = (split /\//, $allele_string)[0];
            
            my $ok = 0;
            my $slice_ref_allele;
            
            # insertion, therefore no ref allele to check
            if($ref_allele eq '-') {
                $ok = 1;
            }
            else {
                my $slice_ref = $slice->sub_Slice($start, $end, $strand);
                
                if(!defined($slice_ref)) {
                    warn "WARNING: Could not fetch sub-slice from $start\-$end\($strand\) on line $line_number" unless defined $config->{quiet};
                }
                
                else {
                    $slice_ref_allele = $slice_ref->seq;
                    $ok = ($slice_ref_allele eq $ref_allele ? 1 : 0);
                }
            }
            
            if(!$ok) {
                warn
                    "WARNING: Specified reference allele $ref_allele ",
                    "does not match Ensembl reference allele",
                    ($slice_ref_allele ? " $slice_ref_allele" : ""),
                    " on line $line_number" unless defined $config->{quiet};
                next;
            }
        }
       
        # create a new VariationFeature object
        my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
          -start => $start,
          -end => $end,
          -slice => $slice,           # the variation must be attached to a slice
          -allele_string => $allele_string,
          -strand => $strand,
          -map_weight => 1,
          -adaptor => $config->{vfa},           # we must attach a variation feature adaptor
          -variation_name => (defined $var_name ? $var_name : $chr.'_'.$start.'_'.$allele_string),
        );
        
        if(defined $config->{whole_genome}) {
            push @{$vf_hash{$chr}{int($start / $config->{chunk_size})}{$start}}, $new_vf;
            $vf_count++;
            $total_vf_count++;
            
            if($vf_count == $config->{buffer_size}) {
                debug("Read $vf_count variants into buffer") unless defined($config->{quiet});
                
                $include_regions ||= &regions_from_hash($config, \%vf_hash);
                
                &check_existing_hash($config, \%vf_hash) if defined($config->{check_existing});
                &whole_genome_fetch($config, \%vf_hash, $transcript_cache, $include_regions);
                
                debug("Processed $total_vf_count total variants") unless defined($config->{quiet});
                
                undef $include_regions unless defined($config->{scan});
                %vf_hash = ();
                $vf_count = 0;
            }
        }
        else {
            &print_consequences($config, [$new_vf]);
            $vf_count++;
            $total_vf_count++;
            debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
        }
      }
    }
    
    # if in whole-genome mode, finish off the rest of the buffer
    if(defined $config->{whole_genome} && %vf_hash) {
        debug("Read $vf_count variants into buffer") unless defined($config->{quiet});
        $include_regions ||= &regions_from_hash($config, \%vf_hash);
        &check_existing_hash($config, \%vf_hash) if defined($config->{check_existing});
        &whole_genome_fetch($config, \%vf_hash, $transcript_cache, $include_regions);
    }
    
    debug("Executed ", defined($Bio::EnsEMBL::DBSQL::StatementHandle::count_queries) ? $Bio::EnsEMBL::DBSQL::StatementHandle::count_queries : 'unknown number of', " SQL statements") if defined($config->{count_queries}) && !defined($config->{quiet});
    
    debug("Finished!") unless defined $config->{quiet};
}

# takes a listref of variation features and prints out consequence information
sub print_consequences {
    my $config = shift;
    my $vfs = shift;
    
    my $out_file_handle = $config->{out_file_handle};
    
    # method name for consequence terms
    my $term_method = $config->{terms}.'_term';
    
    my ($vf_count, $vf_counter);
    $vf_count = scalar @$vfs;
    
    foreach my $new_vf(@$vfs) {
        
        &progress($config, $vf_counter++, $vf_count) unless $vf_count == 1;
        
        # find any co-located existing VFs
        my $existing_vf = $new_vf->{existing};
        $existing_vf ||= &find_existing($config, $new_vf) if defined $config->{check_existing};
        
        # initiate line hash for this variation
        my $line = {
            Uploaded_variation  => $new_vf->variation_name,
            Location            => $new_vf->seq_region_name.':'.&format_coords($new_vf->start, $new_vf->end),
            Existing_variation  => $existing_vf,
            Extra               => {},
        };
        
        # force empty hash into object's transcript_variations if undefined from whole_genome_fetch
        # this will stop the API trying to go off and fill it again
        $new_vf->{transcript_variations} ||= {} if defined $config->{whole_genome};
        
        # regulatory stuff
        if(!defined $config->{coding_only} && defined $config->{regulatory}) {
            
            for my $rfv (@{ $new_vf->get_all_RegulatoryFeatureVariations }) {
                
                my $rf = $rfv->regulatory_feature;
                
                $line->{Feature_type}   = 'RegulatoryFeature';
                $line->{Feature}        = $rf->stable_id;
                
                # this currently always returns 'RegulatoryFeature', so we ignore it for now
                #$line->{Extra}->{REG_FEAT_TYPE} = $rf->feature_type->name;
                
                for my $rfva (@{ $rfv->get_all_alternate_RegulatoryFeatureVariationAlleles }) {
                    
                    $line->{Allele}         = $rfva->variation_feature_seq;
                    $line->{Consequence}    = join ',', 
                        map { $_->$term_method || $_->display_term } 
                            @{ $rfva->get_all_OverlapConsequences };
                            
                    print_line($line);
                }
            }
            
            for my $mfv (@{ $new_vf->get_all_MotifFeatureVariations }) {
                
                my $mf = $mfv->motif_feature;
                
                $line->{Feature_type}   = 'MotifFeature';
                $line->{Feature}        = $mf->binding_matrix->name;
               
                for my $mfva (@{ $mfv->get_all_alternate_MotifFeatureVariationAlleles }) {
                   
                    $line->{Extra}->{MATRIX} = $mf->binding_matrix->description.'_'.$mf->display_label,
                    $line->{Extra}->{MATRIX} =~ s/\s+/\_/g;

                    my $high_inf_pos = $mfva->in_informative_position;

                    if (defined $high_inf_pos) {
                        $line->{Extra}->{HIGH_INF_POS}  = ($high_inf_pos ? 'Y' : 'N');
                    }
                    
                    $line->{Allele}         = $mfva->variation_feature_seq;
                    $line->{Consequence}    = join ',', 
                        map { $_->$term_method || $_->display_term } 
                            @{ $mfva->get_all_OverlapConsequences };
                            
                    print_line($line);
                }
            }
        }
        
        
        # get TVs
        my $tvs = $new_vf->get_all_TranscriptVariations;
        
        # no TVs (intergenic) or only most severe
        if(!@$tvs || defined($config->{most_severe}) || defined($config->{summary})) {
            if(defined($config->{summary})) {
                $line->{Consequence} = join ",", @{$new_vf->consequence_type($config->{terms}) || $new_vf->consequence_type};
            }
            else {
                $line->{Consequence} = $new_vf->display_consequence($config->{terms}) || $new_vf->display_consequence;
            }
            
            &print_line($line);
        }
        
        else {
            foreach my $tv(@$tvs) {
                
                next if(defined $config->{coding_only} && !($tv->affects_transcript));
                
                my $t = $tv->transcript;
                
                $line->{Feature_type}       = 'Transcript';
                $line->{Feature}            = $t->stable_id if defined $t;
                $line->{cDNA_position}      = &format_coords($tv->cdna_start, $tv->cdna_end);
                $line->{CDS_position}       = &format_coords($tv->cds_start, $tv->cds_end);
                $line->{Protein_position}   = &format_coords($tv->translation_start, $tv->translation_end);
                
                # get gene
                my $gene;
                
                if(defined($config->{gene})) {
                    $line->{Gene} = $tv->transcript->{_gene_stable_id};
                    
                    if(!defined($line->{Gene})) {
                        $gene = $config->{ga}->fetch_by_transcript_stable_id($t->stable_id);
                        $line->{Gene}= $gene->stable_id;
                    }
                }
                
                foreach my $tva(@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
                    
                    # basic stuff
                    $line->{Allele}         = $tva->variation_feature_seq;
                    $line->{Amino_acids}    = $tva->pep_allele_string;
                    $line->{Codons}         = $tva->display_codon_allele_string;
                    $line->{Consequence}    = join ",", map {$_->$term_method || $_->display_term} @{$tva->get_all_OverlapConsequences};
                    
                    # HGNC
                    if(defined $config->{hgnc}) {
                        my $hgnc;
                        $hgnc = $tv->transcript->{_gene_hgnc};
                        
                        if(!defined($hgnc)) {
                            if(!defined($gene)) {
                                $gene = $config->{ga}->fetch_by_transcript_stable_id($tv->transcript->stable_id);
                            }
                            
                            my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};
                            if(scalar @entries) {
                                $hgnc = $entries[0]->display_id;
                            }
                        }
                        
                        $hgnc = undef if $hgnc eq '-';
                        
                        $line->{Extra}->{HGNC} = $hgnc if defined($hgnc);
                    }
                    
                    # protein ID
                    if(defined $config->{protein} && $t->translation) {
                        $line->{Extra}->{ENSP} = $t->translation->stable_id;
                    }
                    
                    # HGVS
                    if(defined $config->{hgvs}) {
                        $line->{Extra}->{HGVSc} = $tva->hgvs_coding if defined($tva->hgvs_coding);
                        $line->{Extra}->{HGVSp} = $tva->hgvs_protein if defined($tva->hgvs_protein);
                    }
                    
                    foreach my $tool (qw(SIFT PolyPhen Condel)) {
                        my $lc_tool = lc($tool);
                        
                        if (my $opt = $config->{$lc_tool}) {
                            my $want_pred   = $opt =~ /^p/i;
                            my $want_score  = $opt =~ /^s/i;
                            my $want_both   = $opt =~ /^b/i;
                            
                            if ($want_both) {
                                $want_pred  = 1;
                                $want_score = 1;
                            }
                            
                            next unless $want_pred || $want_score;
                            
                            my $pred_meth   = $lc_tool.'_prediction';
                            my $score_meth  = $lc_tool.'_score';
                            
                            my $pred = $tva->$pred_meth;
                            
                            if($pred) {
                                
                                if ($want_pred) {
                                    $pred =~ s/\s+/\_/;
                                    $line->{Extra}->{$tool} = $pred;
                                }
                                    
                                if ($want_score) {
                                    my $score = $tva->$score_meth;
                                    
                                    if(defined $score) {
                                        if($want_pred) {
                                            $line->{Extra}->{$tool} .= "($score)";
                                        }
                                        else {
                                            $line->{Extra}->{$tool} = $score;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    &print_line($line);
                }
            }
        }
    }
    
    &end_progress($config) unless $vf_count == 1;
}

# prints a line from the hash
sub print_line {
    my $line = shift;

    $line->{Extra} = join ';', map { $_.'='.$line->{Extra}->{$_} } keys %{ $line->{Extra} || {} };

    my $output = join "\t", map { $line->{$_} || '-' } @OUTPUT_COLS;

    my $fh = $config->{out_file_handle};

    print $fh "$output\n";

    # clear out the Extra column for the next line
    $line->{Extra} = {};
}

# sets up configuration hash that is used throughout the script
sub configure {
    my $args = shift;
    
    my $config = {};
    
    GetOptions(
        $config,
        'help',                    # displays help message
        
        # input options,
        'config=s',                # config file name
        'input_file=s',            # input file name
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
        #'no_disconnect',           # disables disconnect_when_inactive
        
        # runtime options
        'most_severe',             # only return most severe consequence
        'summary',                 # only return one line per variation with all consquence types
        'buffer_size=i',           # number of variations to read in before analysis
        'chunk_size=s',            # size in bases of "chunks" used in internal hash structure
        'check_ref',               # check supplied reference allele against DB
        'check_existing',          # find existing co-located variations
        'check_alleles',           # only attribute co-located if alleles are the same
        'failed=i',                # include failed variations when finding existing
        'no_whole_genome',         # disables now default whole-genome mode
        'whole_genome',            # proxy for whole genome mode - now just warns user
        'gp',                      # read coords from GP part of INFO column in VCF (probably only relevant to 1KG)
        'chr=s',                   # analyse only these chromosomes, e.g. 1-5,10,MT
        
        # verbosity options
        'verbose',                 # print out a bit more info while running
        'quiet',                   # print nothing to STDOUT (unless using -o stdout)
        'no_progress',             # don't display progress bars
        
        # output options
        'output_file=s',           # output file name
        'force_overwrite',         # force overwrite of output file if already exists
        'terms=s',                 # consequence terms to use e.g. NCBI, SO
        'coding_only',             # only return results for consequences in coding regions
        'protein',                 # add e! protein ID to extra column
        'hgnc',                    # add HGNC gene ID to extra column
        'hgvs',                    # add HGVS names to extra column
        'sift=s',                  # SIFT predictions
        'polyphen=s',              # PolyPhen predictions
        'condel=s',                # Condel predictions
        'gene',                    # force gene column to be populated (disabled by default, enabled when using cache)
        'regulatory',              # enable regulatory stuff
        
        # cache stuff
        'cache',                   # use cache
        'write_cache',             # enables writing to the cache
        'build=s',                 # builds cache from DB from scratch; arg is either all (all top-level seqs) or a list of chrs
        'scan',                    # scan the whole input file at the beginning to get regions
        'upfront',                 # fetch transcripts and prefetch upfront before analysis starts (requires scan)
        'prefetch',                # prefetch exons, translation, introns, codon table etc for each transcript
        'strip',                   # strips adaptors etc from objects before caching them
        'rebuild=s',               # rebuilds cache by reading in existing then redumping - probably don't need to use this any more
        'dir=s',                   # dir where cache is found (defaults to $HOME/.vep/)
        'cache_region_size=i',     # size of region in bases for each cache file
        'no_slice_cache',          # tell API not to cache features on slice
        'standalone',              # standalone mode uses minimal set of modules installed in same dir, no DB connection
        'skip_db_check',           # don't compare DB parameters with cached
        'compress=s',              # by default we use zcat to decompress; user may want to specify gzcat or "gzip -dc"
        
        # debug
        'cluck',                   # these two need some mods to Bio::EnsEMBL::DBSQL::StatementHandle to work. Clucks callback trace and SQL
        'count_queries',           # counts SQL queries executed
    );
    
    # print usage message if requested or no args supplied
    if(defined($config->{help}) || !$args) {
        &usage;
        exit(0);
    }
    
    # config file?
    if(defined $config->{config}) {
        
        open CONFIG, $config->{config} or die "ERROR: Could not open config file \"".$config->{config}."\"\n";
        
        while(<CONFIG>) {
            next if /^\#/;
            my ($key, $value) = split /\s+|\=/;
            $key =~ s/^\-//g;
            $config->{$key} = $value unless defined $config->{$key};
        }
        
        close CONFIG;
    }

    # can't be both quiet and verbose
    die "ERROR: Can't be both quiet and verbose!" if defined($config->{quiet}) && defined($config->{verbose});
    
    # check file format
    if(defined $config->{format}) {
        die "ERROR: Unrecognised input format specified \"".$config->{format}."\"\n" unless $config->{format} =~ /pileup|vcf|guess/i;
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
    
    # check nsSNP tools
    foreach my $tool(grep {defined $config->{lc($_)}} qw(SIFT PolyPhen Condel)) {
        die "ERROR: Unrecognised option for $tool \"", $config->{lc($tool)}, "\" - must be one of p (prediction), s (score) or b (both)\n" unless $config->{lc($tool)} =~ /^(s|p|b)/;
        
        die "ERROR: $tool not available for this species\n" unless $config->{species} =~ /human|homo/i;
        
        die "ERROR: $tool not available in standalone mode\n" if defined($config->{standalone});
    }
    
    # force quiet if outputting to STDOUT
    if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
        delete $config->{verbose} if defined($config->{verbose});
        $config->{quiet} = 1;
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
            print $key.(' ' x (($max_length - length($key)) + 4)).$config->{$key}."\n";
        }
        
        print "\n".("-" x 20)."\n\n";
    }
    
    # set defaults
    $config->{user}              ||= 'anonymous';
    $config->{buffer_size}       ||= 5000;
    $config->{chunk_size}        ||= '50kb';
    $config->{output_file}       ||= "variant_effect_output.txt";
    $config->{tmpdir}            ||= '/tmp';
    $config->{format}            ||= 'guess';
    $config->{terms}             ||= 'display';
    $config->{gene}              ||= 1 unless defined($config->{whole_genome});
    $config->{cache_region_size} ||= 1000000;
    $config->{dir}          ||= join '/', ($ENV{'HOME'}, '.vep');
    $config->{compress}          ||= 'zcat';
    
    # warn users still using whole_genome flag
    if(defined($config->{whole_genome})) {
        debug("INFO: Whole-genome mode is now the default run-mode for the script. To disable it, use --no_whole_genome") unless defined($config->{quiet});
    }
    
    $config->{whole_genome}      = 1 unless defined $config->{no_whole_genome};
    $config->{include_failed}    = 1 unless defined $config->{include_failed};
    $config->{chunk_size}        =~ s/mb?/000000/i;
    $config->{chunk_size}        =~ s/kb?/000/i;
    $config->{cache_region_size} =~ s/mb?/000000/i;
    $config->{cache_region_size} =~ s/kb?/000/i;
    
    # cluck and display executed SQL?
    $Bio::EnsEMBL::DBSQL::StatementHandle::cluck = 1 if defined($config->{cluck});
    
    # standalone needs cache, can't use HGVS
    if(defined($config->{standalone})) {
        $config->{cache} = 1;
        
        die("ERROR: Cannot generate HGVS coordinates in standalone mode") if defined($config->{hgvs});
        
        die("ERROR: Cannot analyse regulatory features in standalone mode") if defined($config->{regulatory});
    }
    
    # no_slice_cache, prefetch and whole_genome have to be on to use cache or upfront
    if(defined($config->{cache}) || defined($config->{upfront})) {
        $config->{prefetch} = 1;
        $config->{no_slice_cache} = 1;
        $config->{whole_genome} = 1;
        $config->{strip} = 1;
        
        # scan should also be on for upfront
        $config->{scan} = 1 if defined($config->{upfront});
    }
    
    $config->{build} = $config->{rebuild} if defined($config->{rebuild});
    
    # force options for full build
    if(defined($config->{build})) {
        $config->{prefetch} = 1;
        $config->{gene} = 1;
        $config->{hgnc} = 1;
        $config->{no_slice_cache} = 1;
        $config->{cache} = 1;
        $config->{strip} = 1;
        $config->{write_cache} = 1;
    }
    
    # connect to databases
    $config->{reg} = &connect_to_dbs($config);
    
    # complete dir with species name and db_version
    $config->{dir} .= '/'.(
        join '/', (
            defined($config->{standalone}) ? $config->{species} : ($config->{reg}->get_alias($config->{species}) || $config->{species}),
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
    
    # suppress warnings that the FeatureAdpators spit if using no_slice_cache
    Bio::EnsEMBL::Utils::Exception::verbose(1999) if defined($config->{no_slice_cache});
    
    # get adaptors
    if(defined($config->{cache})) {
        
        # try and load adaptors from cache
        if(!&load_dumped_adaptor_cache($config)) {
            &get_adaptors($config);
            &dump_adaptor_cache($config) if defined($config->{write_cache});
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
                
                # but we still need to reconnect
                debug("INFO: Defined host ", $config->{host}, " is different from cached ", $dbc->{_host}, " - reconnecting to host") unless defined($config->{quiet});
                
                &get_adaptors($config);
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
        &dump_adaptor_cache($config) if defined($config->{write_cache})
    }
    
    # get terminal width for progress bars
    unless(defined($config->{quiet})) {
        my $width;
        
        # module may not be installed
        eval {
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
        
        # build the cache
        debug("Building cache for ".$config->{species}) unless defined($config->{quiet});
        &build_full_cache($config, $config->{rebuild});
        
        # exit script
        debug("Finished building cache") unless defined($config->{quiet});
        exit(0);
    }
    
    # warn user DB will be used for SIFT/PolyPhen/Condel
    if(defined($config->{cache}) && (defined($config->{sift}) || defined($config->{polyphen}) || defined($config->{condel}) || defined($config->{hgvs}) || defined($config->{regulatory}))) {
        debug("INFO: Database will be accessed for SIFT/PolyPhen/Condel, HGVS and regulatory features") unless defined($config->{quiet});
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
    
    # configure output file
    $config->{out_file_handle} = &get_out_file_handle($config);
    
    return $config;
}

# connects to DBs; in standalone mode this just loads registry module
sub connect_to_dbs {
    my $config = shift;
    
    # get registry
    my $reg = 'Bio::EnsEMBL::Registry';
    
    unless(defined($config->{standalone})) {
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
    
    $config->{vfa} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variationfeature');
    $config->{tva} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'transcriptvariation');
    
    # get fake ones for species with no var DB
    if(!defined($config->{vfa})) {
        $config->{vfa} = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($config->{species});
    }
    else {
        $config->{vfa}->db->include_failed_variations($config->{include_failed}) if defined($config->{vfa}->db) && $config->{vfa}->db->can('include_failed_variations');
    }
    
    $config->{sa}  = $config->{reg}->get_adaptor($config->{species}, 'core', 'slice');
    $config->{ga}  = $config->{reg}->get_adaptor($config->{species}, 'core', 'gene');
    $config->{ta}  = $config->{reg}->get_adaptor($config->{species}, 'core', 'transcript');
    $config->{mca} = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
    $config->{csa} = $config->{reg}->get_adaptor($config->{species}, 'core', 'coordsystem');
    
    # cache schema version
    $config->{mca}->get_schema_version if defined $config->{mca};
    
    # check we got slice adaptor - can't continue without a core DB
    die("ERROR: Could not connect to core database\n") unless defined $config->{sa};
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
            $in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{in_file}, "\n");
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
    
    if($config->{output_file} =~ /stdout/i) {
        $out_file_handle = *STDOUT;
    }
    else {
        $out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
    }
    
    # make header
    my $time = &get_time;
    my $db_string = $config->{mca}->dbc->dbname." on ".$config->{mca}->dbc->host if defined $config->{mca};
    $db_string .= "\n## Using cache in ".$config->{dir} if defined($config->{cache});
    my $version_string =
        "Using API version ".$config->{reg}->software_version.
        ", DB version ".(defined $config->{mca} && $config->{mca}->get_schema_version ? $config->{mca}->get_schema_version : '?');
    
    my $header =<<HEAD;
## ENSEMBL VARIANT EFFECT PREDICTOR v$VERSION
## Output produced at $time
## Connected to $db_string
## $version_string
## Extra column keys:
## HGNC         : HGNC gene identifier
## ENSP         : Ensembl protein identifer
## HGVSc        : HGVS coding sequence name
## HGVSp        : HGVS protein sequence name
## SIFT         : SIFT prediction
## PolyPhen     : PolyPhen prediction
## Condel       : Condel SIFT/PolyPhen consensus prediction
## MATRIX       : The source and identifier of a transcription factor binding profile aligned at this position
## HIGH_INF_POS : A flag indicating if the variant falls in a high information position of a transcription factor binding profile
HEAD
    
    # add headers
    print $out_file_handle $header;
    
    # add column headers
    print $out_file_handle '#', (join "\t", @OUTPUT_COLS);
    print $out_file_handle "\n";
    
    return $out_file_handle;
}

# parses a line of input
sub parse_line {
    my $config = shift;
    my $line   = shift;
    
    my @data = (split /\s+/, $_);
    
    # pileup: chr1 60 T A
    if(
       ($config->{format} =~ /pileup/i) ||
       (
            $data[0] =~ /(chr)?\w+/ &&
            $data[1] =~ /\d+/ &&
            $data[2] =~ /^[ACGTN-]+$/ &&
            $data[3] =~ /^[ACGTNRYSWKM*+\/-]+$/
        )
    ) {
        my @return = ();
        
        if($data[2] ne "*"){
            my $var;
            
            if($data[3] =~ /^[A|C|G|T]$/) {
                $var = $data[3];
            }
            else {
                ($var = unambiguity_code($data[3])) =~ s/$data[2]//ig;
            }
            if(length($var)==1){
                push @return, [$data[0], $data[1], $data[1], $data[2]."/".$var, 1, undef];
            }
            else{
                for my $nt(split //,$var){
                    push @return, [$data[0], $data[1], $data[1], $data[2]."/".$nt, 1, undef];
                }
            }
        }
        else{ #indel
            my @genotype=split /\//,$data[3];
            foreach my $allele(@genotype){
                if(substr($allele,0,1) eq "+") { #ins
                    push @return, [$data[0], $data[1]+1, $data[1], "-/".substr($allele,1), 1, undef];
                }
                elsif(substr($allele,0,1) eq "-"){ #del
                    push @return, [$data[0], $data[1], $data[1]+length($data[3])-4, substr($allele,1)."/-", 1, undef];
                }
                elsif($allele ne "*"){
                    warn("WARNING: invalid pileup indel genotype: $line\n") unless defined $config->{quiet};
                    push @return, ['non-variant'];
                }
            }
        }
        return \@return;
    }
    
    # VCF: 20      14370   rs6054257 G     A      29    0       NS=58;DP=258;AF=0.786;DB;H2          GT:GQ:DP:HQ
    elsif(
        ($config->{format} =~ /vcf/i) ||
        (
            $data[0] =~ /(chr)?\w+/ &&
            $data[1] =~ /\d+/ &&
            $data[3] =~ /^[ACGTN-]+$/ &&
            $data[4] =~ /^([\.ACGTN-]+\,?)+$/
        )
    ) {
        
        # non-variant line in VCF, return dummy line
        if($data[4] eq '.') {
            return [['non-variant']];
        }
        
        # get relevant data
        my ($chr, $start, $end, $ref, $alt) = ($data[0], $data[1], $data[1], $data[3], $data[4]);
        
        if(defined $config->{gp}) {
            $chr = undef;
            $start = undef;
            
            foreach my $pair(split /\;/, $data[7]) {
                my ($key, $value) = split /\=/, $pair;
                if($key eq 'GP') {
                    ($chr,$start) = split /\:/, $value;
                    $end = $start;
                }
            }
            
            unless(defined($chr) and defined($start)) {
                warn "No GP flag found in INFO column" unless defined $config->{quiet};
                return [['non-variant']];
            }
        }
        
        # adjust end coord
        $end += (length($ref) - 1);
        
        # find out if any of the alt alleles make this an insertion or a deletion
        my ($is_indel, $is_sub, $ins_count, $total_count);
        foreach my $alt_allele(split /\,/, $alt) {
            $is_indel = 1 if $alt_allele =~ /D|I/;
            $is_indel = 1 if length($alt_allele) != length($ref);
            $is_sub = 1 if length($alt_allele) == length($ref);
            $ins_count++ if length($alt_allele) > length($ref);
            $total_count++;
        }
        
        # multiple alt alleles?
        if($alt =~ /\,/) {
            if($is_indel) {
                
                my @alts;
                
                if($alt =~ /D|I/) {
                    foreach my $alt_allele(split /\,/, $alt) {
                        # deletion (VCF <4)
                        if($alt_allele =~ /D/) {
                            push @alts, '-';
                        }
                        
                        elsif($alt_allele =~ /I/) {
                            $alt_allele =~ s/^I//g;
                            push @alts, $alt_allele;
                        }
                    }
                }
                
                else {
                    $ref = substr($ref, 1);
                    $ref = '-' if $ref eq '';
                    $start++;
                    
                    foreach my $alt_allele(split /\,/, $alt) {
                        $alt_allele = substr($alt_allele, 1);
                        $alt_allele = '-' if $alt_allele eq '';
                        push @alts, $alt_allele;
                    }
                }
                
                $alt = join "/", @alts;
            }
            
            else {
                # for substitutions we just need to replace ',' with '/' in $alt
                $alt =~ s/\,/\//;
            }
        }
        
        else {
            if($is_indel) {
                # deletion (VCF <4)
                if($alt =~ /D/) {
                    my $num_deleted = $alt;
                    $num_deleted =~ s/\D+//g;
                    $end += $num_deleted - 1;
                    $alt = "-";
                    $ref .= ("N" x ($num_deleted - 1)) unless length($ref) > 1;
                }
                
                # insertion (VCF <4)
                elsif($alt =~ /I/) {
                    $ref = '-';
                    $alt =~ s/^I//g;
                    $start++;
                }
                
                # insertion or deletion (VCF 4+)
                else {
                    # chop off first base
                    $ref = substr($ref, 1);
                    $alt = substr($alt, 1);
                    
                    $start++;
                    
                    if($ref eq '') {
                        # make ref '-' if no ref allele left
                        $ref = '-';
                    }
                    
                    # make alt '-' if no alt allele left
                    $alt = '-' if $alt eq '';
                }
            }
        }
        
        return [[$chr, $start, $end, $ref."/".$alt, 1, ($data[2] eq '.' ? undef : $data[2])]];
        
    }
    
    # our format
    else {
        # we allow commas as delimiter so re-split
        @data = (split /\s+|\,/, $_);
        return [\@data];
    }
}

# takes a hash of VFs and fetches consequences by pre-fetching overlapping transcripts
# from database and/or cache
sub whole_genome_fetch {
    my $config = shift;
    my $vf_hash = shift;
    my $transcript_cache = shift;
    my $include_regions = shift;
    
    my $up_down_size = MAX_DISTANCE_FROM_TRANSCRIPT;
    
    my (%vf_done, @finished_vfs, %seen_trs);
    
    # convert regions to cached sizes
    my $converted_regions = &convert_regions($config, $include_regions) if defined($config->{cache});
    
    foreach my $chr(sort {$a <=> $b} keys %$vf_hash) {
        if(defined($config->{standalone}) && !-e $config->{dir}.'/'.$chr) {
            debug("No cache found for chromsome $chr") unless defined($config->{quiet});
            next;
        }
        
        my $slice_cache;
        
        debug("Analyzing chromosome $chr") unless defined($config->{quiet});
        
        my $use_regions = defined($config->{cache}) ? $converted_regions : $include_regions;
        my ($count_from_db, $count_from_cache, $count_duplicates) = (0, 0, 0);
        
        if(!defined($transcript_cache->{$chr})) {
            
            # no regions defined (this probably shouldn't happen)
            if(!defined($use_regions->{$chr})) {
                
                # spoof regions covering whole chromosome
                my $start = 1;
                my $end = $config->{cache_region_size};
                my $slice = &get_slice($config, $chr);
                
                if(defined($slice)) {
                    while($start < $slice->end) {
                        push @{$use_regions->{$chr}}, $start.'-'.$end;
                        $start += $config->{cache_region_size};
                        $end += $config->{cache_region_size};
                    }
                }
            }
            
            # check we have defined regions
            if(defined($use_regions->{$chr})) {
                my $region_count = scalar @{$use_regions->{$chr}};
                my $counter;
                
                debug("Reading transcript data from cache and/or database") unless defined($config->{quiet});
                
                foreach my $region(sort {(split /\-/, $a)[0] <=> (split /\-/, $b)[1]} @{$use_regions->{$chr}}) {
                    &progress($config, $counter++, $region_count);
                    
                    # skip regions beyond the end of the chr
                    next if defined($slice_cache->{$chr}) && (split /\-/, $region)[0] > $slice_cache->{$chr}->length;
                    
                    # force quiet so other methods don't mess up the progress bar
                    my $quiet = $config->{quiet};
                    $config->{quiet} = 1;
                    
                    # try and load cache from disk if using cache
                    my $tmp_cache;
                    if(defined($config->{cache})) {
                        $tmp_cache = &load_dumped_transcript_cache($config, $chr, $region);
                        $count_from_cache += scalar @{$tmp_cache->{$chr}} if defined($tmp_cache->{$chr});
                    }
                    
                    # no cache found on disk or not using cache
                    if(!defined($tmp_cache->{$chr})) {
                        
                        if(defined($config->{standalone})) {
                            debug("WARNING: Could not find cache for $chr\:$region") unless defined($config->{quiet});
                            next;
                        }
                        
                        # spoof temporary region hash
                        my $tmp_hash;
                        push @{$tmp_hash->{$chr}}, $region;
                        
                        $tmp_cache = &cache_transcripts($config, $tmp_hash);
                        
                        # make it an empty arrayref that gets cached
                        # so we don't get confused and reload next time round
                        $tmp_cache->{$chr} ||= [];
                        
                        $count_from_db += scalar @{$tmp_cache->{$chr}};
                        
                        # dump to disk if writing to cache
                        &dump_transcript_cache($config, $tmp_cache, $chr, $region) if defined($config->{write_cache});
                    }
                    
                    # add loaded transcripts to main cache
                    if(defined($tmp_cache->{$chr})) {
                        while(my $tr = shift @{$tmp_cache->{$chr}}) {
                            
                            # track already added transcripts by dbID
                            my $dbID = $tr->dbID;
                            if($seen_trs{$dbID}) {
                                $count_duplicates++;
                                next;
                            }
                            $seen_trs{$dbID} = 1;
                            
                            push @{$transcript_cache->{$chr}}, $tr;
                        }
                    }
                    
                    undef $tmp_cache;
                    
                    # restore quiet status
                    $config->{quiet} = $quiet;
                    
                    # build slice cache
                    $slice_cache = &build_slice_cache($config, $transcript_cache) unless defined($slice_cache->{$chr});
                }
                
                &end_progress($config);
            }
        }
        
        # skip chr if no cache
        next unless defined($transcript_cache->{$chr});
        
        # copy slice from transcript to slice cache
        $slice_cache = &build_slice_cache($config, $transcript_cache) unless defined($slice_cache->{$chr});
        
        my $tr_count = scalar @{$transcript_cache->{$chr}};
        
        debug("Retrieved $tr_count transcripts ($count_from_cache cached, $count_from_db DB, $count_duplicates duplicates)") unless defined($config->{quiet});
        debug("Analyzing variants") unless defined($config->{quiet});
        
        my $tr_counter;
        
        while($tr_counter < $tr_count) {
            
            &progress($config, $tr_counter, $tr_count);
            
            my $tr = $transcript_cache->{$chr}->[$tr_counter++];
            
            # do each overlapping VF
            my $s = $tr->start - $up_down_size;
            my $e = $tr->end + $up_down_size;
            
            # get the chunks this transcript overlaps
            my %chunks;
            $chunks{$_} = 1 for (int($s/$config->{chunk_size})..int($e/$config->{chunk_size}));
            map {delete $chunks{$_} unless defined($vf_hash->{$chr}{$_})} keys %chunks;
            
            foreach my $chunk(keys %chunks) {
                foreach my $pos(grep {$_ >= $s && $_ <= $e} keys %{$vf_hash->{$chr}{$chunk}}) {
                    foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {
                        
                        # pinch slice from slice cache if we don't already have it
                        $_->{slice} ||= $slice_cache->{$chr} for @{$vf_hash->{$chr}{$chunk}{$pos}};
                        
                        my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
                            -transcript => $tr,
                            -variation_feature => $vf,
                            -adaptor => $config->{tva},
                            -no_ref_check => 1
                        );
                        
                        # prefetching stuff here prevents doing loads at the
                        # end and makes progress reporting more useful
                        $tv->_prefetch_for_vep;
                        
                        $vf->add_TranscriptVariation($tv);
                    }
                }
            }
        }
        
        # sort results into @finished_vfs array
        foreach my $chunk(sort {$a <=> $b} keys %{$vf_hash->{$chr}}) {
            foreach my $pos(sort {$a <=> $b} keys %{$vf_hash->{$chr}{$chunk}}) {
                
                # pinch slice from slice cache if we don't already have it
                $_->{slice} ||= $slice_cache->{$chr} for @{$vf_hash->{$chr}{$chunk}{$pos}};
                
                # add to final array
                push @finished_vfs, @{$vf_hash->{$chr}{$chunk}{$pos}};
            }
        }
        
        &end_progress($config);
        
        debug("Calculating and writing output") unless defined($config->{quiet});
        &print_consequences($config, \@finished_vfs);
        
        # clean hash
        delete $vf_hash->{$chr};
        
        delete $transcript_cache->{$chr} if defined($config->{cache});
    }
}

# gets existing VFs for a vf_hash
sub check_existing_hash {
    my $config = shift;
    my $vf_hash = shift;
    my $variation_cache;
    
    debug("Checking for existing variations") unless defined($config->{quiet});
    
    my ($chunk_count, $counter);
    $chunk_count += scalar keys %{$vf_hash->{$_}} for keys %{$vf_hash};
    
    foreach my $chr(keys %{$vf_hash}) {
        
        my %loaded_regions;
        
        foreach my $chunk(keys %{$vf_hash->{$chr}}) {
            &progress($config, $counter++, $chunk_count);
            
            # get the VFs for this chunk
            my ($start, $end);
            
            # work out start and end using chunk_size
            $start = $config->{chunk_size} * $chunk;
            $end = $config->{chunk_size} * ($chunk + 1);
            
            # using cache?
            if(defined($config->{cache})) {
                my $tmp_regions;
                push @{$tmp_regions->{$chr}}, $start.'-'.$end;
                
                my $converted_regions = &convert_regions($config, $tmp_regions);
                
                foreach my $region(@{$converted_regions->{$chr}}) {
                
                    unless($loaded_regions{$region}) {
                        my $tmp_cache = &load_dumped_variation_cache($config, $chr, $region);
                        
                        # load from DB if not found in cache
                        if(!defined($tmp_cache->{$chr})) {
                            if(defined($config->{standalone})) {
                                debug("WARNING: Could not find variation cache for $chr\:$region") unless defined($config->{quiet});
                                next;
                            }
                            
                            $tmp_cache->{$chr} = &get_variations_in_region($config, $chr, $region);
                            &dump_variation_cache($config, $tmp_cache, $chr, $region) if defined($config->{write_cache});
                        }
                        
                        # merge tmp_cache with the main cache
                        foreach my $key(keys %{$tmp_cache->{$chr}}) {
                            $variation_cache->{$chr}->{$key} = $tmp_cache->{$chr}->{$key};
                            delete $tmp_cache->{$chr}->{$key};
                        }
                        
                        # clear memory
                        undef $tmp_cache;
                        
                        # record this region as fetched
                        $loaded_regions{$region} = 1;
                    }
                }
            }
            
            # no cache, get all variations in region from DB
            else {
                $variation_cache->{$chr} = &get_variations_in_region($config, $chr, $start.'-'.$end);
            }
            
            # now compare retrieved vars with vf_hash
            foreach my $pos(keys %{$vf_hash->{$chr}->{$chunk}}) {
                foreach my $var(@{$vf_hash->{$chr}->{$chunk}->{$pos}}) {
                    my @found;
                    
                    if(defined($variation_cache->{$chr})) {
                        if(my $existing_vars = $variation_cache->{$chr}->{$pos}) {
                            foreach my $existing_var(@$existing_vars) {
                                push @found, $existing_var->[0] unless &is_var_novel($config, $existing_var, $var);
                            }
                        }
                    }
                    
                    $var->{existing} = join ",", @found;
                    $var->{existing} ||= '-';
                }
            }
        }
        
        delete $variation_cache->{$chr};
    }
    
    &end_progress($config);
}

# gets a slice from the slice adaptor
sub get_slice {
    my $config = shift;
    my $chr = shift;
    
    return undef unless defined($config->{sa}) && defined($chr);
    
    my $slice;
    
    # first try to get a chromosome
    eval { $slice = $config->{sa}->fetch_by_region('chromosome', $chr); };
    
    # if failed, try to get any seq region
    if(!defined($slice)) {
        $slice = $config->{sa}->fetch_by_region(undef, $chr);
    }
    
    return $slice;
}




# METHODS THAT DEAL WITH "REGIONS"
##################################

# scans file to get all slice bits we need
sub scan_file() {
    my $config = shift;
    
    my $in_file_handle = $config->{in_file_handle};
    
    my %include_regions;
    
    debug("Scanning input file") unless defined($config->{quiet});
    
    while(<$in_file_handle>) {
        chomp;
        
        # header line?
        next if /^\#/;
        
        # some lines (pileup) may actually parse out into more than one variant)
        foreach my $sub_line(@{&parse_line($config, $_)}) {
        
            # get the sub-line into named variables
            my ($chr, $start, $end, $allele_string, $strand, $var_name) = @{$sub_line};
            $chr =~ s/chr//ig unless $chr =~ /^chromosome$/i;
            $chr = 'MT' if $chr eq 'M';
            
            next if defined($config->{chr}) && !$config->{chr}->{$chr};
            
            $include_regions{$chr} ||= [];
            
            &add_region($start, $end, $include_regions{$chr});
        }
    }
    
    # close filehandle and recycle
    close $in_file_handle;
    $config->{in_file_handle} = &get_in_file_handle($config);
    
    # merge regions
    &merge_regions(\%include_regions);
    
    return \%include_regions;
}

# gets regions from VF hash
sub regions_from_hash {
    my $config = shift;
    my $vf_hash = shift;
    
    my %include_regions;
    
    # if using cache we just want the regions of cache_region_size
    # since that's what we'll get from the cache (or DB if no cache found)
    if(defined($config->{cache})) {
        
        my $region_size = $config->{cache_region_size};
        
        foreach my $chr(keys %$vf_hash) {
            $include_regions{$chr} = [];
            my %temp_regions;
            
            foreach my $chunk(keys %{$vf_hash->{$chr}}) {
                foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                    my ($s, $e) = ($pos - MAX_DISTANCE_FROM_TRANSCRIPT, $pos + MAX_DISTANCE_FROM_TRANSCRIPT);
                    
                    my $low = int ($s / $region_size);
                    my $high = int ($e / $region_size) + 1;
                    
                    for my $i($low..($high - 1)) {
                        $temp_regions{(($i * $region_size) + 1).'-'.(($i + 1) * $region_size)} = 1;
                    }
                }
            }
            
            @{$include_regions{$chr}} = keys %temp_regions;
        }
    }
    
    # if no cache we don't want to fetch more than is necessary, so find the
    # minimum covered region of the variations in the hash
    else {
        foreach my $chr(keys %$vf_hash) {
            $include_regions{$chr} = [];
            
            foreach my $chunk(keys %{$vf_hash->{$chr}}) {
                foreach my $pos(keys %{$vf_hash->{$chr}{$chunk}}) {
                    &add_region($_->start, $_->end, $include_regions{$chr}) for @{$vf_hash->{$chr}{$chunk}{$pos}};
                }
            }
        }
        
        # merge regions
        &merge_regions(\%include_regions);
    }
    
    return \%include_regions;
}

# adds a region to region list, expanding existing one if overlaps
sub add_region {
    my $start = shift;
    my $end = shift;
    my $region_list = shift;
    
    # fix end for insertions
    $end = $start if $end < $start;
    
    my $added = 0;
    my $i = 0;
    
    while ($i < scalar @$region_list) {
        my ($region_start, $region_end) = split /\-/, $region_list->[$i];
        
        if($start <= $region_end && $end >= $region_start) {
            my $new_region_start = ($start < $end ? $start : $end) - MAX_DISTANCE_FROM_TRANSCRIPT;
            my $new_region_end = ($start > $end ? $start : $end) + MAX_DISTANCE_FROM_TRANSCRIPT;
            
            $region_start = $new_region_start if $new_region_start < $region_start;
            $region_end = $new_region_end if $new_region_end > $region_end;
            
            $region_list->[$i] = $region_start.'-'.$region_end;
            $added = 1;
        }
        
        $i++;
    }
    
    unless($added) {
        push @{$region_list}, ($start - MAX_DISTANCE_FROM_TRANSCRIPT).'-'.($end + MAX_DISTANCE_FROM_TRANSCRIPT);
    }
}

# merges overlapping regions from scans
sub merge_regions {
    my $include_regions = shift;
    
    # now merge overlapping regions
    foreach my $chr(keys %$include_regions) {
        my $max_index = $#{$include_regions->{$chr}};
        my (@new_regions, %skip);
        
        for my $i(0..$max_index) {
            next if $skip{$i};
            my ($s, $e) = split /\-/, $include_regions->{$chr}[$i];
            
            for my $j(($i+1)..$max_index) {
                next if $skip{$j};
                my ($ns, $ne) = split /\-/, $include_regions->{$chr}[$j];
                
                if($s <= $ne && $e >= $ns) {
                    $s = $ns if $ns < $s;
                    $e = $ne if $ne > $e;
                    
                    $skip{$j} = 1;
                }
            }
            
            push @new_regions, $s.'-'.$e;
        }
        
        # replace original
        $include_regions->{$chr} = \@new_regions;
        
        $config->{region_count} += scalar @new_regions;
    }
    
    return $include_regions;
}

# converts regions as determined by scan_file to regions loadable from cache
sub convert_regions {
    my $config = shift;
    my $regions = shift;
    
    return undef unless defined $regions;
    
    my $region_size = $config->{cache_region_size};
    
    my %new_regions;
    
    foreach my $chr(keys %$regions) {
        my %temp_regions;
        
        foreach my $region(@{$regions->{$chr}}) {
            my ($s, $e) = split /\-/, $region;
            
            my $low = int ($s / $region_size);
            my $high = int ($e / $region_size) + 1;
            
            for my $i($low..($high - 1)) {
                $temp_regions{(($i * $region_size) + 1).'-'.(($i + 1) * $region_size)} = 1;
            }
        }
        
        @{$new_regions{$chr}} = keys %temp_regions;
    }
    
    return \%new_regions;
}





# CACHE METHODS
###############

# get transcripts for slices
sub cache_transcripts {
    my $config = shift;
    my $include_regions = shift;
    
    my $transcript_cache;
    my $i;
    
    debug("Caching transcripts") unless defined($config->{quiet});
    
    foreach my $chr(keys %$include_regions) {
        
        my $slice = &get_slice($config, $chr);
        
        next unless defined $slice;
        
        # prefetch some things
        $slice->is_circular;
        
        # trim bumf off the slice
        delete $slice->{coord_system}->{adaptor} if defined($config->{write_cache});
        
        # no regions?
        if(!scalar @{$include_regions->{$chr}}) {
            my $start = 1;
            my $end = $config->{cache_region_size};
            
            while($start < $slice->end) {
                push @{$include_regions->{$chr}}, $start.'-'.$end;
                $start += $config->{cache_region_size};
                $end += $config->{cache_region_size};
            }
        }
        
        my $region_count;
        
        if(scalar keys %$include_regions == 1) {
            my ($chr) = keys %$include_regions;
            $region_count = scalar @{$include_regions->{$chr}};
            debug("Caching transcripts for chromosome $chr") unless defined($config->{quiet});
        }
        
        foreach my $region(@{$include_regions->{$chr}}) {
            &progress($config, $i++, $region_count || $config->{region_count});
            
            my ($s, $e) = split /\-/, $region;
            
            # sanity check start and end
            $s = 1 if $s < 1;
            $e = $slice->end if $e > $slice->end;
            
            # get sub-slice
            my $sub_slice = $slice->sub_Slice($s, $e);
            
            # add transcripts to the cache, via a transfer to the chrom's slice
            if(defined($sub_slice)) {
                foreach my $gene(@{$sub_slice->get_all_Genes(undef, undef, 1)}) {
                    my $gene_stable_id = $gene->stable_id;
                    
                    foreach my $tr(map {$_->transfer($slice)} @{$gene->get_all_Transcripts}) {
                        $tr->{_gene_stable_id} = $gene_stable_id;
                        
                        if(defined($config->{prefetch})) {
                            $tr->{_gene} = $gene;
                            &prefetch_transcript_data($config, $tr);
                            delete $tr->{_gene};
                        }
                        
                        # strip some unnecessary data from the transcript object
                        &clean_transcript($tr) if defined($config->{write_cache});
                        
                        push @{$transcript_cache->{$chr}}, $tr;
                    }
                }
            }
        }
    }
    
    &end_progress($config);
    
    return $transcript_cache;
}

# gets rid of extra bits of info attached to the transcript that we don't need
sub clean_transcript {
    my $tr = shift;
    
    foreach my $key(qw(display_xref external_db external_display_name external_name external_status created_date status description edits_enabled modified_date)) {
        delete $tr->{$key} if defined($tr->{$key});
    }
    
    # clean all attributes but miRNA
    if(defined($tr->{attributes})) {
        my @new_atts;
        foreach my $att(@{$tr->{attributes}}) {
            push @new_atts, $att if $att->{code} eq 'miRNA';
        }
        $tr->{attributes} = \@new_atts;
    }
    
    $tr->{analysis} = {};
}

# build slice cache from transcript cache
sub build_slice_cache {
    my $config = shift;
    my $transcript_cache = shift;
    
    my %slice_cache;
    
    foreach my $chr(keys %$transcript_cache) {
        $slice_cache{$chr} = $transcript_cache->{$chr}[0]->slice;
        
        # reattach adaptor to the coord system
        $slice_cache{$chr}->{coord_system}->{adaptor} ||= $config->{csa};
    }
    
    return \%slice_cache;
}

# pre-fetches per-transcript data
sub prefetch_transcript_data {
    my $config = shift;
    my $tran = shift;
    
    # introns, translateable_seq, mapper
    $tran->{_variation_effect_feature_cache}->{introns} ||= $tran->get_all_Introns;
    $tran->{_variation_effect_feature_cache}->{translateable_seq} ||= $tran->translateable_seq;
    $tran->{_variation_effect_feature_cache}->{mapper} ||= $tran->get_TranscriptMapper;
    
    # peptide
    unless ($tran->{_variation_effect_feature_cache}->{peptide}) {
        my $translation = $tran->translate;
        $tran->{_variation_effect_feature_cache}->{peptide} = $translation ? $translation->seq : undef;
    }
    
    # codon table
    unless ($tran->{_variation_effect_feature_cache}->{codon_table}) {
        # for mithocondrial dna we need to to use a different codon table
        my $attrib = $tran->slice->get_all_Attributes('codon_table')->[0];
        
        $tran->{_variation_effect_feature_cache}->{codon_table} = $attrib ? $attrib->value : 1;
    }
    
    # gene HGNC
    if(defined $config->{hgnc}) {
        # get from gene cache if found already
        if(defined($tran->{_gene}->{_hgnc})) {
            $tran->{_gene_hgnc} = $tran->{_gene}->{_hgnc};
        }
        else {
            my @entries = grep {$_->database eq 'HGNC'} @{$tran->{_gene}->get_all_DBEntries()};
            if(scalar @entries) {
                $tran->{_gene_hgnc} = $entries[0]->display_id;
            }
            
            $tran->{_gene_hgnc} ||= '-';
            
            # cache it on the gene object too
            $tran->{_gene}->{_hgnc} = $tran->{_gene_hgnc};
        }
    }
    
    return $tran;
}

# dumps out transcript cache to file
sub dump_transcript_cache {
    my $config = shift;
    my $transcript_cache = shift;
    my $chr = shift;
    my $region = shift;
    
    debug("Dumping cached transcript data") unless defined($config->{quiet});
    
    # clean the slice adaptor before storing
    &clean_slice_adaptor($config);
    
    &strip_transcript_cache($config, $transcript_cache);
    
    $config->{reg}->disconnect_all;
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.($region || "dump").'.gz';
    
    # make directory if it doesn't exist
    if(!(-e $dir)) {
        system("mkdir -p ".$dir);
    }
    
    debug("Writing to $dump_file") unless defined($config->{quiet});
    
    # storable
    open my $fh, "| gzip -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
    nstore_fd($transcript_cache, $fh);
    close $fh;
}

# loads in dumped transcript cache to memory
sub load_dumped_transcript_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.($region || "dump").'.gz';
    
    return undef unless -e $dump_file;
    
    debug("Reading cached transcript data for chromosome $chr".(defined $region ? "\:$region" : "")." from dumped file") unless defined($config->{quiet});
    
    open my $fh, $config->{compress}." ".$dump_file." |" or return undef;
    my $transcript_cache = fd_retrieve($fh);
    close $fh;
    
    return $transcript_cache;
}

# strips cache
sub strip_transcript_cache {
    my $config = shift;
    my $cache = shift;
    
    foreach my $chr(keys %$cache) {
        foreach my $tr(@{$cache->{$chr}}) {
            foreach my $exon(@{$tr->{_trans_exon_array}}) {
                delete $exon->{adaptor};
                delete $exon->{slice}->{adaptor};
            }
            
            delete $tr->{adaptor};
            delete $tr->{slice}->{adaptor};
        }
    }
}

# cleans slice adaptor before storing in cache
sub clean_slice_adaptor{
    my $config = shift;
    
    # clean some stuff off the slice adaptor
    $config->{sa}->{asm_exc_cache} = {};
    $config->{sa}->{sr_name_cache} = {};
    $config->{sa}->{sr_id_cache} = {};
    delete $config->{sa}->{db}->{seq_region_cache};
    delete $config->{sa}->{db}->{name_cache};
}


# dump adaptors to cache
sub dump_adaptor_cache {
    my $config = shift;
    
    $config->{reg}->disconnect_all;
    
    my $dir = $config->{dir};
    my $dump_file = $dir.'/adaptors.gz';
    
    # make directory if it doesn't exist
    if(!(-e $dir)) {
        system("mkdir -p ".$dir);
	}
	
    open my $fh, "| gzip -c > ".$dump_file or die "ERROR: Could not write to dump file $dump_file";
    nstore_fd($config, $fh);
    close $fh;
}

# load dumped adaptors
sub load_dumped_adaptor_cache {
    my $config = shift;
    
    my $dir = $config->{dir};
    my $dump_file = $dir.'/adaptors.gz';
    
    return undef unless -e $dump_file;
    
    debug("Reading cached adaptor data") unless defined($config->{quiet});
    
    open my $fh, $config->{compress}." ".$dump_file." |" or return undef;
    my $cached_config = fd_retrieve($fh);
    close $fh;
    
    $config->{$_} = $cached_config->{$_} for qw(sa ga ta vfa tva mca csa);
    
    return 1;
}

# dumps cached variations to disk
sub dump_variation_cache {
    my $config = shift;
    my $v_cache = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.($region || "dump").'_var.gz';
    
    # make directory if it doesn't exist
    if(!(-e $dir)) {
        system("mkdir -p ".$dir);
    }
    
    open DUMP, "| gzip -c > ".$dump_file or die "ERROR: Could not write to adaptor dump file $dump_file";
    
    foreach my $pos(keys %{$v_cache->{$chr}}) {
        foreach my $v(@{$v_cache->{$chr}->{$pos}}) {
            my ($name, $source, $start, $end, $as, $strand) = @$v;
            
            print DUMP join " ", (
                $name,
                $source == 1 ? '' : $source,
                $start,
                $end == $start ? '' : $end,
                $as,
                $strand == 1 ? '' : $strand,
            );
            print DUMP "\n";
        }
    }
    
    close DUMP;    
}

# loads dumped variation cache
sub load_dumped_variation_cache {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my $dir = $config->{dir}.'/'.$chr;
    my $dump_file = $dir.'/'.($region || "dump").'_var.gz';
    
    return undef unless -e $dump_file;
    
    open DUMP, $config->{compress}." ".$dump_file." |" or return undef;
    
    my $v_cache;
    
    while(<DUMP>) {
        chomp;
        my ($name, $source, $start, $end, $as, $strand) = split / /, $_;
        $source ||= 1;
        $end ||= $start;
        $strand ||= 1;
        
        my @v = ($name, $source, $start, $end, $as, $strand);
        push @{$v_cache->{$chr}->{$start}}, \@v;
    }
    
    close DUMP;
    
    return $v_cache;
}

# builds a full cache for this species
sub build_full_cache {
    my $config = shift;
    my $rebuild = shift;
    
    my @slices;
    
    if($config->{build} =~ /all/i) {
        @slices = @{$config->{sa}->fetch_all('toplevel')};
    }
    else {
        foreach my $val(split /\,/, $config->{build}) {
            my @nnn = split /\-/, $val;
            
            foreach my $chr($nnn[0]..$nnn[-1]) {
                my $slice = &get_slice($config, $chr);
                push @slices, $slice if defined($slice);
            }
        }
    }
    
    foreach my $slice(@slices) {
        my $chr = $slice->seq_region_name;
        
        my $regions;
        
        # for progress
        my $region_count = int($slice->end / $config->{cache_region_size}) + 1;
        my $counter = 0;
        
        # initial region
        my ($start, $end) = (1, $config->{cache_region_size});
        
        debug((defined($config->{rebuild}) ? "Rebuild" : "Creat")."ing cache for chromosome $chr") unless defined($config->{quiet});
        
        while($start < $slice->end) {
            
            &progress($config, $counter++, $region_count);
            
            # store quiet status
            my $quiet = $config->{quiet};
            $config->{quiet} = 1;
            
            # store transcripts
            $regions->{$chr} = [$start.'-'.$end];
            my $tmp_cache = ($rebuild ? &load_dumped_transcript_cache($config, $chr, $start.'-'.$end) : &cache_transcripts($config, $regions));
            $tmp_cache->{$chr} ||= [];
            
            &dump_transcript_cache($config, $tmp_cache, $chr, $start.'-'.$end);
            undef $tmp_cache;
            
            # store variations
            my $variation_cache;
            $variation_cache->{$chr} = &get_variations_in_region($config, $chr, $start.'-'.$end);
            $variation_cache->{$chr} ||= {};
            
            &dump_variation_cache($config, $variation_cache, $chr, $start.'-'.$end);
            undef $variation_cache;
            
            # restore quiet status
            $config->{quiet} = $quiet;
            
            # increment by cache_region_size to get next region
            $start += $config->{cache_region_size};
            $end += $config->{cache_region_size};
        }
        
        &end_progress($config);
        
        undef $regions;
    }
}

# format coords for printing
sub format_coords {
    my ($start, $end) = @_;
    
    if(!defined($start)) {
        return '-';
    }
    elsif(!defined($end)) {
        return $start;
    }
    elsif($start == $end) {
        return $start;
    }
    elsif($start > $end) {
        return $end.'-'.$start;
    }
    else {
        return $start.'-'.$end;
    }
}




# METHODS TO FIND CO-LOCATED / EXISTING VARIATIONS
##################################################

# finds an existing VF in the db
sub find_existing {
    my $config = shift;
    my $new_vf = shift;
    
    if(defined($new_vf->adaptor->db)) {
        
        my $sth = $new_vf->adaptor->db->dbc->prepare(qq{
            SELECT variation_name, source_id, seq_region_start, seq_region_end, allele_string, seq_region_strand
            FROM variation_feature
            WHERE seq_region_id = ?
            AND seq_region_start = ?
            AND seq_region_end = ?
            ORDER BY source_id ASC
        });
        
        $sth->execute($new_vf->slice->get_seq_region_id, $new_vf->start, $new_vf->end);
        
        my @v;
        for my $i(0..5) {
            $v[$i] = undef;
        }
        
        $sth->bind_columns(\$v[0], \$v[1], \$v[2], \$v[3], \$v[4], \$v[5]);
        
        my @found;
        
        while($sth->fetch) {
            push @found, $v[0] unless &is_var_novel($config, \@v, $new_vf);
        }
        
        $sth->finish();
        
        return (scalar @found ? join ",", @found : undef);
    }
    
    return undef;
}

# compare a new vf to one from the cache / DB
sub is_var_novel {
    my $config = shift;
    my $existing_var = shift;
    my $new_var = shift;
    
    my $is_novel = 1;
    
    $is_novel = 0 if $existing_var->[2] == $new_var->start && $existing_var->[3] == $new_var->end;
    
    if(defined($config->{check_alleles})) {
        my %existing_alleles;
        
        $existing_alleles{$_} = 1 for split /\//, $existing_var->[4];
        
        my $seen_new = 0;
        foreach my $a(split /\//, $new_var->allele_string) {
            reverse_comp(\$a) if $new_var->seq_region_strand ne $existing_var->[5];
            $seen_new = 1 unless defined $existing_alleles{$a};
        }
        
        $is_novel = 1 if $seen_new;
    }
    
    return $is_novel;
}

# gets all variations in a region
sub get_variations_in_region {
    my $config = shift;
    my $chr = shift;
    my $region = shift;
    
    my ($start, $end) = split /\-/, $region;
    
    my %variations;
    
    if(defined($config->{vfa}->db)) {
        my $sth = $config->{vfa}->db->dbc->prepare(qq{
            SELECT vf.variation_name, vf.source_id, vf.seq_region_start, vf.seq_region_end, vf.allele_string, vf.seq_region_strand
            FROM variation_feature vf, seq_region s
            WHERE s.seq_region_id = vf.seq_region_id
            AND s.name = ?
            AND vf.seq_region_start >= ?
            AND vf.seq_region_start <= ?
        });
        
        $sth->execute($chr, $start, $end);
        
        my @v;
        for my $i(0..5) {
            $v[$i] = undef;
        }
        
        $sth->bind_columns(\$v[0], \$v[1], \$v[2], \$v[3], \$v[4], \$v[5]);
        
        while($sth->fetch) {
            my @v_copy = @v;
            push @{$variations{$v[2]}}, \@v_copy;
        }
        
        $sth->finish();
    }
    
    return \%variations;
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

# update or initiate progress bar
sub progress {
    my ($config, $i, $total) = @_;
    
    return if defined($config->{quiet}) || defined($config->{no_progress});
    
    my $width = $config->{terminal_width};
    my $percent = int(($i/$total) * 100);
    my $numblobs = (($i/$total) * $width) - 2;
    
    # this ensures we're not writing to the terminal too much
    return if(defined($config->{prev_prog})) && $numblobs.'-'.$percent eq $config->{prev_prog};
    $config->{prev_prog} = $numblobs.'-'.$percent;
    
    printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

# end progress bar
sub end_progress {
    my $config = shift;
    return if defined($config->{quiet}) || defined($config->{no_progress});
    &progress($config, 1,1);
    print "\n";
    delete $config->{prev_prog};
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

-i | --input_file      Input file - if not specified, reads from STDIN. Files
                       may be gzip compressed.
--format               Alternative input file format - one of "pileup", "vcf"
-o | --output_file     Output file. Write to STDOUT by specifying -o STDOUT - this
                       will force --quiet [default: "variant_effect_output.txt"]
--force_overwrite      Force overwriting of output file [default: quit if file
                       exists]

-t | --terms           Type of consequence terms to output - one of "ensembl", "SO",
                       "NCBI" [default: ensembl]
 
--sift=[p|s|b]         Add SIFT [p]rediction, [s]core or [b]oth [default: off]
--polyphen=[p|s|b]     Add PolyPhen [p]rediction, [s]core or [b]oth [default: off]
--condel=[p|s|b]       Add Condel SIFT/PolyPhen consensus [p]rediction, [s]core or
                       [b]oth [default: off]

NB: SIFT, PolyPhen and Condel predictions are currently available for human only

--regulatory           Look for overlaps with regulatory regions. The script can
                       also call if a variant falls in a high information position
                       within a transcription factor binding site. Output lines have
                       a Feature type of RegulatoryFeature or MotifFeature. Requires
                       database connection. [default: off]
                       
NB: Regulatory consequences are currently available for human and mouse only

--hgnc                 If specified, HGNC gene identifiers are output alongside the
                       Ensembl Gene identifier [default: off]
--hgvs                 Output HGVS identifiers (coding and protein). Requires database
                       connection [default: off]
--protein              Output Ensembl protein identifer [default: off]
--gene                 Force output of Ensembl gene identifer - disabled by default
                       unless using --cache or --no_whole_genome [default: off]

--coding_only          Only return consequences that fall in the coding region of
                       transcripts [default: off]
--most_severe          Ouptut only the most severe consequence per variation.
                       Transcript-specific columns will be left blank. [default: off]
--summary              Ouptut only a comma-separated list of all consequences per
                       variation. Transcript-specific columns will be left blank.
                       [default: off]

--check_ref            If specified, checks supplied reference allele against stored
                       entry in Ensembl Core database [default: off]
--check_existing       If specified, checks for existing co-located variations in the
                       Ensembl Variation database [default: off]
--check_alleles        If specified, the alleles of existing co-located variations
                       are compared to the input; an existing variation will only
                       be reported if no novel allele is in the input (strand is
                       accounted for) [default: off]

--chr [list]           Select a subset of chromosomes to analyse from your file. Any
                       data not on this chromosome in the input will be skipped. The
                       list can be comma separated, with "-" characters representing
                       an interval [default: off]
--gp                   If specified, tries to read GRCh37 position from GP field in the
                       INFO column of a VCF file. Only applies when VCF is the input
                       format and human is the species [default: off]

--species              Species to use [default: "human"]
--host                 Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]
--genomes              Sets DB connection params for Ensembl Genomes [default: off]
--registry             Registry file to use defines DB connections [default: off]
                       Defining a registry file overrides above connection settings.
--db_version=[number]  Force script to load DBs from a specific Ensembl version. Not
                       advised due to likely incompatibilities between API and DB

--no_whole_genome      Run in old-style, non-whole genome mode [default: off]
--buffer_size          Sets the number of variants sent in each batch [default: 5000]
                       Increasing buffer size can retrieve results more quickly
                       but requires more memory. Only applies to whole genome mode.
                       
--cache                Enables read-only use of cache [default: off]
--dir [directory]      Specify the base cache directory to use [default: "\$HOME/.vep/"]
--write_cache          Enable writing to cache [default: off]
--build [all|list]     Build a complete cache for the selected species. Build for all
                       chromosomes with --build all, or a list of chromosomes (see
                       --chr). DO NOT USE WHEN CONNECTED TO PUBLIC DB SERVERS AS THIS
                       VIOLATES OUR FAIR USAGE POLICY [default: off]
                       
--compress             Specify utility to decompress cache files - may be "gzcat" or
                       "gzip -dc Only use if default does not work [default: zcat]
                       
--skip_db_check        ADVANCED! Force the script to use a cache built from a different
                       database than specified with --host. Only use this if you are
                       sure the hosts are compatible (e.g. ensembldb.ensembl.org and
                       useastdb.ensembl.org) [default: off]
--cache_region_size    ADVANCED! The size in base-pairs of the region covered by one
                       file in the cache. [default: 1MB]
END

    print $usage;
}
