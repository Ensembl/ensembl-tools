#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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


# Region Report tool
#     A script for sampling a given set of chromosomal regions, producing a simple
# summary of the features within those regions.
# Operates in batch and single region mode, and can serialize into GFF3 and textual format.
# 

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::IO::GFFSerializer;
use Bio::EnsEMBL::Utils::IO::ReportSerializer;
use Bio::EnsEMBL::Utils::IO::FASTASerializer;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Getopt::Long;
use IO::File;
use Bio::EnsEMBL::Utils::IO qw(:slurp);

sub usage {
    print <<'USAGE'; 
    Example: region_report.pl --species=Human --chromosome=x --start=133017695 --end=133161157
    
    By default, all genetic features are supplied. More control is given with:
        --exclude=
        --include=
      Choose from the following list for inclusions or exclusions
        g - genes, t - transcripts, p - translations, 
        v - variations, s - structural variants, c- constrained elements,
        r - regulatory features, q - raw sequence (off by default)

      e.g. --exclude=rcv means everything except these three types of feature.
    Also configure database options with --host --port --db_version --user and --password
    Batch processing can be performed with --file option to read a list of regions from a file
    
    For more help, run region_report.pl --help
    
USAGE
}

## Parse command line input

my %config;
GetOptions ( \%config, 
       'example',                 # trigger for test case
       'species:s', 
       'chromosome:s',
       'strand:i',
       'start:i',                 # genomic coordinate
       'end:i',                   # genomic coordinate
       'include:s',               # features to include in report
       'exclude:s',               # features to exclude from report
             
       'host=s',                  # database host
       'port=s',                  # database port
       'user=s',                  # database user name
       'password=s',              # database password
       'db_version=i',            # Ensembl database version to use e.g. 62
       
       'secondaryhost=s',         # secondary database host
       'secondaryport=s',         # secondary database port
       'secondaryuser=s',         # secondary database user name
       'secondarypassword=s',     # secondary database password

       'input=s',                 # Input file containing many regions to parse in batch.
       'output=s',                # Target file to receive report

       'gff3',                    # Format selection for output
       'report',                  # Text mode output 
             
       'verbose',
       'help',
       'h',
       'usage',
       );

if ($config{'help'} or $config{'h'}) { help();}
if ($config{'usage'}) { usage(); exit;}
# Prep for the coming storm
my $registry = 'Bio::EnsEMBL::Registry';
my @regions;
my $serializer;
my $output_fh;

if ($config{'output'}) {
    $output_fh = IO::File->new($config{'output'},'w') or die "Unable to open requested output file: ".$config{'output'}." $!\n";
}
# Select the correct serializer based on user choice


# Build up a stack of regions to process into the one output file.
if ($config{'input'}) {
    print STDERR "[Region report] Reading regions from ",$config{'input'},"\n" if ($config{'verbose'});
    @regions = @{slurp_to_array($config{'input'},1)};
}

if ($config{'example'}) {
    $config{'chromosome'} = 6;
    $config{'start'} = 133017695;
    $config{'end'} = 133161157;
    $config{'strand'} = 1;
    $config{'species'} = "Human";
}
elsif (!$config{'species'}) {
    usage;
    die "No species specified\n";
}

if ($config{'start'} and $config{'end'} and $config{'chromosome'}) {
    push @regions,$config{'chromosome'}.":".$config{'start'}."..".$config{'end'};
}

# Translate include/exclude instructions from user:
my %feature_types = ( 'g' => 1, 't' => 1, 'p' => 0, 'v' => 1, 's' => 1, 'c' => 1, 'r' => 1, 'q' => 0);
if ($config{'include'}) {
    %feature_types = ( 'g' => 0, 't' => 0, 'p' => 0, 'v' => 0, 's' => 0, 'c' => 0, 'r' => 0, 'q' => 0);
    if ($config{'verbose'} and $config{'include'} =~ /[^gtpvscrq]/) { print STDERR "Unrecognised feature requested in include= option\n";}
    my @options = split (//,$config{'include'});
    foreach (@options) { $feature_types{$_} = 1};
}
if ($config{'exclude'}) {
    if ($config{'verbose'} and $config{'exclude'} =~ /[^gtpvscrq]/) { print STDERR "Unrecognised feature requested in exclude= option\n";}
    my @options = split (//,$config{'exclude'});
    foreach (@options) { $feature_types{$_} = 0};
}
if ($config{'verbose'}) {
    print STDERR "[Region report] # Selected options: ";
    foreach (keys(%feature_types)) {print STDERR $_;}
    print STDERR "\n[Region report] #                   ";
    foreach (keys(%feature_types)) {if ($feature_types{$_} == 1) {print STDERR "x"} else {print STDERR " "} ;}
    print STDERR "\n";
}
# Supply defaults for DB options

if (not defined ($config{'host'}) ) { $config{'host'} = 'ensembldb.ensembl.org';}
if (not defined ($config{'user'}) ) { $config{'user'} = 'anonymous';}
if (not defined ($config{'db_version'})) { $config{'db_version'} = software_version();}
if (not defined ($config{'port'})) { $config{'port'} = 5306;}
if (not defined ($config{'password'})) { $config{'password'} = "";}

print STDERR "[Region report] # Connecting to ",$config{'host'},":",$config{'port'}," as ",$config{'user'},"\n" if ($config{'verbose'});
my @dbs;
push(@dbs, {
  -HOST => $config{host},
  -PORT => $config{port},
  -USER => $config{'user'},
  -VERBOSE => 0,
  -DB_VERSION => $config{db_version},
});
$dbs[0]->{-PASS} = $config{password} if $config{password};

if($config{secondaryhost}) {
  if (not defined ($config{'secondaryhost'}) ) { $config{'secondaryhost'} = 'ensembldb.ensembl.org';}
  if (not defined ($config{'secondaryuser'}) ) { $config{'secondaryuser'} = 'anonymous';}
  if (not defined ($config{'secondaryport'})) { $config{'secondaryport'} = 5306;}
  if (not defined ($config{'secondarypassword'})) { $config{'secondarypassword'} = "";}
  push(@dbs, {
    -HOST => $config{secondaryhost},
    -PORT => $config{secondaryport},
    -USER => $config{secondaryuser},
    -VERBOSE => 0,
    -DB_VERSION => $config{db_version},
  });
  $dbs[1]->{-PASS} = $config{secondarypassword} if $config{secondarypassword};
}

$registry->load_registry_from_multiple_dbs(@dbs);


# Create a serializer for output formatting.
my $ontology_adaptor; #useful for GFF serializer.

if ($config{'gff3'} or $config{'report'}) {
    if ($config{'gff3'}) {
        $serializer = "Bio::EnsEMBL::Utils::IO::GFFSerializer";
        $ontology_adaptor = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
        $serializer = new $serializer($ontology_adaptor,$output_fh);
    }
    if ($config{'report'}) {
        $serializer = "Bio::EnsEMBL::Utils::IO::ReportSerializer";
        $serializer = new $serializer($output_fh);
    }
}
else { #default to GFF for output
    print STDERR "[Region report] Defaulting to GFF3 output format\n" if ($config{'verbose'});
    $serializer = "Bio::EnsEMBL::Utils::IO::GFFSerializer";
    $ontology_adaptor = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
    $serializer = new $serializer($ontology_adaptor,$output_fh);
}

# Now iterate over supplied regions and turn them into Slice objects

my $slice_adaptor;
if ($config{'special_db'}) {
    my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                    -species => $config{'species'},
                    -group   => 'core',
                    -user    => $config{'user'},
                    -dbname  => $config{'special_db'},
                    -host    => $config{'host'},
                    -password => $config{'password'},
                    -port => $config{'port'},
                  );
     $slice_adaptor = $db->get_SliceAdaptor();
     
}
else {
    $slice_adaptor = $registry->get_adaptor( $config{'species'}, 'Core', 'Slice' );
}



my @slices;
foreach my $region (@regions) {
  $region =~ s/\s//g;
  next unless $region; #empty line
    print STDERR "[Region report] Requested region: $region\n" if ($config{'verbose'});
    
    eval {
        my $slice = $slice_adaptor->fetch_by_toplevel_location($region);
        push @slices, $slice;
    };
    if ($@) {
        print STDERR $@."\n";
        exit 2; # exit
    }
    
}
print STDERR "[Region report] Found ".scalar(@slices)." slices\n" if ($config{'verbose'});
my $dba = $slice_adaptor->db();
$serializer->print_main_header(\@slices, $dba);

# Explore all slices for features, check for validity etc.

foreach my $slice (@slices) {

    unless ($slice->seq_region_name and $slice->start and $slice->end)  {&usage; die "Please specify a species, chromosome, and start and end coordinates. The supplied slices could not be found\n";}
    my $running_total;
    my $transcript_list;

    if ($slice->is_circular) {
        $serializer->print_feature( $slice);
        # an explicit region anchor is not necessary under most circumstances, except where a circular sequence needs describing
    }
    if ($feature_types{'g'} ) {
        my $gene_list = $slice->get_all_Genes();
        $serializer->print_metadata("Genes in the region ".$slice->start."-".$slice->end." of chromosome ".$slice->seq_region_name.":");
        $serializer->print_feature_list($gene_list);
    }
    if ($feature_types{'t'} and not $config{'report'}) {
        $transcript_list = $slice->get_all_Transcripts();
        $serializer->print_metadata(" Transcripts in the region ".$slice->start."-".$slice->end." of chromosome ".$slice->seq_region_name.":");
        $serializer->print_feature_list($transcript_list);
    }
    if ($feature_types{'r'}) {
# Funcgen stuff requires a separate adaptor. Get all the regulatory features.
        my $funcgen_feature_adaptor;
        $funcgen_feature_adaptor = $registry->get_adaptor( $config{'species'}, 'funcgen', 'regulatoryfeature');
        
        $serializer->print_metadata("Regulatory features:");
# 10 kilobases equates to about 30MB of ram allocated. The code default is 1 megabase! Choose a smaller chunk size
        if ($config{'report'}) {
            my $feature_list = $funcgen_feature_adaptor->fetch_all_by_Slice($slice);
            $serializer->print_feature_list($feature_list);
        } else {
            my $iterator = $funcgen_feature_adaptor->fetch_Iterator_by_Slice($slice,undef,60000);
            $serializer->print_feature_Iterator($iterator);
        }
    }
    
    
    if ($feature_types{'s'} ) {        
        my $structural_variant_list;
        $structural_variant_list = $slice->get_all_StructuralVariationFeatures();
        $serializer->print_metadata("Structural Variations:");
        $serializer->print_feature_list($structural_variant_list)
    }
    if ($feature_types{'v'} ) {
        my $variation_feature_adaptor;
        $variation_feature_adaptor = $registry->get_adaptor( $config{'species'}, 'variation', 'variationfeature');
        $serializer->print_metadata("Variation Features:");
        if ($config{'report'}) {
            my $feature_list = $variation_feature_adaptor->fetch_all_by_Slice($slice);
            $serializer->print_feature_list($feature_list);
        } else {
            my $iterator = $variation_feature_adaptor->fetch_Iterator_by_Slice($slice,undef,60000);
            $serializer->print_feature_Iterator($iterator);
        }
    }
    if ($feature_types{'c'} ) {
        my $method_link_adaptor;
        $method_link_adaptor = $registry->get_adaptor( 'Multi', 'compara', 'MethodLinkSpeciesSet' );
        my $method_list = $method_link_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSTRAINED_ELEMENT","mammals");
        my $constrained_element_adaptor = $registry->get_adaptor('Multi','compara','ConstrainedElement');
        my $element_list = $constrained_element_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($method_list,$slice);
        $serializer->print_metadata("Constrained Elements across 34 eutherian mammals:");
        $serializer->print_feature_list($element_list);
    }
}
# FASTA dumper. Use with caution, it makes for big files.
if ($feature_types{'q'} ) {
    $serializer->print_metadata("#FASTA\n"); #Print method provides first #
    foreach my $slice (@slices) {
        my $fasta_serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($output_fh);
        $fasta_serializer->print_Seq($slice);
    }
}

unless ($serializer->printed_something) {
    exit 1; # output file does not contain any features.
}
exit 0; # normal execution

sub help {
        
    print <<'MANUAL';
    ### The Region Report Tool ###
    
    The Region Report tool produces summaries of small sections of EnsEMBL data in
    readable form. It writes a list of all features in a chosen genomic slice given
    the following parameters:
    
        Species
        Chromosome
        Start position
        End position
        Strand
        
    In addition you may choose which features to print out from this selection:
    
        Gene
        Transcript
        Translation
        Variation (EnsEMBL Variation)
        Structural Variant (EnsEMBL Variation)
        Constrained Element (EnsEMBL Compara)
        Regulatory Feature (EnsEMBL Funcgen)
        FASTA sequence dump
        
    At present the Region Report tool can produce both human-readable text and 
    GFF3 formats.
    
    ### Note on large region reports ###
    
    You may feel tempted to use this tool to dump large amounts of GFF, but we ask
    that you refrain from doing so if using EnsEMBL database servers. In particular,
    Variation features have very high feature density. As a matter of courtesy, we 
    ask that you aim to keep your selected regions under 200kb. Web users will
    experience limits depending on their choice of features.
    
    ### Dependencies ###
    
    To use the Region Report tool, the client must have the installed the following:
    
    ensembl-core
    ensembl-funcgen
    ensembl-variation
    ensembl-compara
    
    Use only stable releases and their corresponding data to avoid potential problems.
    
    Also ensure that the modules folder from each bundle is added to your Perl library 
    path ($PERL5LIB on linux).
    
    ### Output destination ###
    
    By default the Region Report tool prints to the console, and can be redirected or
    piped in the normal way. Named file output can also be declared by using the 
    --output=filename command.
    
    ### Using a different database server ###
    
    Should you wish to use a server other than EnsEMBL's, further command line options
    can be specified, e.g.
    
    region_report.pl --host=local-mirror-server --user=me --password=go
    
    Sometimes it is necessary to split data across several servers. In these
    situations the options --secondaryhost, --secondarypassword, --secondaryuser and
    --secondaryport can be used to find features from two databases at once.
    
    ### Operating in batch mode ###
    
    If you have many different regions to explore at once, they can be run in a batch.
    Create a text file containing the regions of interest, one per line like so:
    
        5:100000-110000
        7:1-150000
        x:17000000-17100000
        y:10000..20000
        
    The first number or letter is the chromosome, followed by the start and end of the
    region of interest. The region can be specified as one would when using the
    EnsEMBL genome browser.
    
    Save your file and use it for input, using the --input=batch_file command. The
    output will be combined into one large single report.
    
    
    ### The future ###
    
    This tool is under development. For feature requests and bug reports, please
    contact the mailing list http://lists.ensembl.org/mailman/listinfo/dev or http://www.ensembl.org/Help/Contact
    
MANUAL
}
