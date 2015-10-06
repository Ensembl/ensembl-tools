use strict;
use warnings;
BEGIN { $| = 1;
	use Test::More;
	#plan tests => 6;
}

use FindBin qw($RealBin);
use File::Copy;
use Data::Dumper;
use lib $RealBin.'/../';
use Bio::EnsEMBL::Variation::Utils::VEP;

#use Bio::EnsEMBL::Test::TestUtils;

# setup variables
my ($input, $output, $full_output, $expected, @lines);
my $tmpfile = "$$\_test_vep_input";

# configure script
open CONF, "$RealBin\/test.conf" or die "ERROR: Could not read from conf file $RealBin\/test.conf\n";
my $config;
while(<CONF>) {
  chomp;
  my ($k, $v) = split;
  $config->{$k} = $v;
}
close CONF;

# find where ensembl-variation is installed
my $mod_path = 'Bio/EnsEMBL/Variation/Utils/VEP.pm';
my $var_path = $INC{$mod_path};
die("ERROR: Could not find path to ensembl-variation modules\n") unless $var_path;
$var_path =~ s/$mod_path//;
my $data_path = $var_path.'/t/testdata/';
die("ERROR: Could not find test data path $data_path\n") unless -d $data_path;

my $ver    = $config->{version};
my $cver   = $config->{cache_version};
my $ass    = $config->{assembly};
my $sp     = $config->{species};
my $script = $RealBin.'/../variant_effect_predictor.pl';
my $perl   = '/usr/bin/env perl';
my $inc    = '-I ~/Variation/modules/ -I $RealBin\../';
my $bascmd = "$perl $inc $script";
my $cmd    = "$bascmd -force -offline -dir_cache $data_path/vep-cache/ -i $tmpfile -o stdout -cache_version $cver -assembly $ass -species $sp";

# unzip fasta
`gzip -dc $data_path/vep-cache/$sp/$cver\_$ass/test.fa.gz > $data_path/vep-cache/$sp/$cver\_$ass/test.fa` if(-e "$data_path/vep-cache/$sp/$cver\_$ass/test.fa.gz");


## BASIC TEST
#############

$output = `$cmd --help`;
ok($output =~ /ENSEMBL VARIANT EFFECT PREDICTOR/, "help message");


## HEADERS FOR OUTPUT FORMATS
#############################

# setup input
$input = '21 25587758 rs116645811 G A . . .';
input($input);

# vcf
$output = (grep {/^\#/} (split "\n", `$cmd -vcf`))[-1];
$expected = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
ok($output eq $expected, "output header - vcf") or diag("Expected\n$expected\n\nGot\n$output");

# default format
$full_output = `$cmd`;
$output = (grep {/^\#/} (split "\n", $full_output))[-1];
$expected =
  "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\t".
  "Feature_type\tConsequence\tcDNA_position\tCDS_position\t".
  "Protein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";

ok($output eq $expected, "output header - default") or diag("Expected\n$expected\n\nGot\n$output");


## output formats
input(qq{21 25606454 25606454 G/C + test});

# vcf
$output = `$cmd --vcf | grep -v '#'`;
$expected = '21\s25606454\stest\sG\sC\s.\s.\sCSQ=C|ENSG00000154719|ENST00000419219|Transcript|missense_variant|284|275|92|A/G|gCc/gGc|||-1,C|ENSG00000154719|ENST00000352957|Transcript|missense_variant|317|275|92|A/G|gCc/gGc|||-1,C|ENSG00000154719|ENST00000307301|Transcript|missense_variant|317|275|92|A/G|gCc/gGc|||-1';
ok($output =~ /$expected/, "VCF output") or diag("Expected\n$expected\n\nGot\n$output");

# json
eval q{ use JSON; };

if(!$@) {
  $output = `$cmd --json`;
  ok($output =~ /^\{.*"most_severe_consequence":"missense_variant".*\}\n?$/, "JSON output");
}
else {
  print STDERR "# could not find JSON perl module, skipping test\n";
}

# GVF
$output = `$cmd --gvf`;
$expected = "21\tUser\tSNV\t25606454\t25606454";
chomp $output;
ok($output =~ /$expected/, "GVF detail") or diag("Expected\n$expected\n\nGot\n$output");

my $exp_gvf_hash = {
  ID => { 1 => 1 },
  Variant_seq => { 'C' => 1 },
  Dbxref => { 'User:test' => 1},
  Reference_seq => { 'G' => 1 },
  Variant_effect => {
    'missense_variant 0 mRNA ENST00000419219' => 1,
    'missense_variant 0 mRNA ENST00000352957' => 1,
    'missense_variant 0 mRNA ENST00000307301' => 1,
  },
};

my $out_hash;
foreach my $chunk(split(";", (split("\t", $output))[-1])) {
  my ($k, $v) = split("=", $chunk);
  $out_hash->{$k}->{$v} = 1;
}

is_deeply($out_hash, $exp_gvf_hash, "GVF consequence");

# custom fields
$output = `$cmd --fields Uploaded_variation,Feature,Consequence --pick | grep -v '#'`;
$expected = 'test\sENST00000307301\smissense_variant';
ok($output =~ /$expected/, "custom fields") or diag("Expected\n$expected\n\nGot\n$output");

# convert VCF
$output = `$cmd --convert vcf -o stdout | grep -v '#'`;
$expected = '21\s25606454\stest\sG\sC';
ok($output =~ /$expected/, "convert to VCF") or diag("Expected\n$expected\n\nGot\n$output");

# convert pileup
$output = `$cmd --convert pileup -o stdout | grep -v '#'`;
$expected = '21\s25606454\sG\sC';
ok($output =~ /$expected/, "convert to pileup") or diag("Expected\n$expected\n\nGot\n$output");

# individual
$input = qq{##fileformat=VCFv4.0
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT A B
21 25587758 rs116645811 G A . . . GT 1|1 0|0};
input($input);
$output = `$cmd --individual all`;
ok($output =~ /IND=A/ && $output !~ /IND=B/ && $output =~ /ZYG=HOM/, "individual");

# html output
move('stdout.html', 'stdout.html.bak'.$$) if -e 'stdout.html';
input(qq{21 25606454 25606454 G/C + test});
$output = `$cmd --html`;

ok(-e 'stdout.html', "HTML output exists");

open HTML, "stdout.html";
@lines = <HTML>;
ok((grep {/\<title\>VEP output\<\/title\>/} @lines), "HTML output content");

unlink('stdout.html');
move('stdout.html.bak'.$$, 'stdout.html') if -e 'stdout.html.bak'.$$;



## CONSEQUENCE TYPES
####################

# all consequence types should be covered by API tests
# just test a few to be sure

# missense
$output = (grep {/missense/} (split "\n", $full_output))[0];
$expected =
  "ENSG00000154719\tENST00000307301\tTranscript\tmissense_variant\t".
  "1043\t1001\t334\tT/M\taCg/aTg";

ok($output =~ /$expected/, "consequence type - missense") or diag("Expected\n$expected\n\nGot\n$output");

# upstream
$output = (grep {/upstream/} (split "\n", $full_output))[0];
$expected =
  "ENSG00000260583\tENST00000567517\tTranscript\tupstream_gene_variant\t".
  "-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=4432";

ok($output =~ /$expected/, "consequence type - upstream") or diag("Expected\n$expected\n\nGot\n$output");

# intron
$output = (grep {/intron/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000352957\tTranscript\tintron_variant";

ok($output =~ /$expected/, "consequence type - intron_variant") or diag("Expected\n$expected\n\nGot\n$output");


## VARIANT TYPES
################

# do all in one run to save time
$input = qq{21 25587758 25587758 C/T - rev
21 25587758 25587758 G/- + del
21 25587759 25587758 -/T + ins
21 25587758 25587760 GTA/TCC + multi
21 25587758 25587760 GTA/TC + unbalanced
21 25585656 25607517 DEL + svd
21 25585656 25607517 DUP + svp
};
input($input);
$full_output = `$cmd`;

# reverse strand
$output = (grep {/rev.+missense/} (split "\n", $full_output))[0];
$expected =
  "ENSG00000154719\tENST00000307301\tTranscript\tmissense_variant\t".
  "1043\t1001\t334\tT/M\taCg/aTg";
ok($output =~ /$expected/, "variant type - reverse strand") or diag("Expected\n$expected\n\nGot\n$output");

# deletion
$output = (grep {/del.+frameshift/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tframeshift_variant\t1043\t1001\t334";
ok($output =~ /$expected/, "variant type - deletion") or diag("Expected\n$expected\n\nGot\n$output");

# insertion
$output = (grep {/ins.+frameshift/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tframeshift_variant\t1042-1043\t1000-1001\t334";
ok($output =~ /$expected/, "variant type - insertion") or diag("Expected\n$expected\n\nGot\n$output");

# multi-bp sub
$output = (grep {/multi.+missense/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tmissense_variant\t1041-1043\t999-1001\t333-334\tTT/TE\tacTACg/acGGAg";
ok($output =~ /$expected/, "variant type - multi-bp sub") or diag("Expected\n$expected\n\nGot\n$output");

# unbalanced sub
$output = (grep {/unbalanced.+frameshift/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\tframeshift_variant\t1041-1043\t999-1001\t333-334";
ok($output =~ /$expected/, "variant type - unbalanced sub") or diag("Expected\n$expected\n\nGot\n$output");

# sv - deletion
$output = (grep {/svd.+ENST00000307301/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\ttranscript_ablation";
ok($output =~ /$expected/, "variant type - structural variation - deletion") or diag("Expected\n$expected\n\nGot\n$output");

# sv - duplication
$output = (grep {/svp.+ENST00000307301/} (split "\n", $full_output))[0];
$expected = "ENSG00000154719\tENST00000307301\tTranscript\ttranscript_amplification";
ok($output =~ /$expected/, "variant type - structural variation - duplication") or diag("Expected\n$expected\n\nGot\n$output");


## OPTIONS
##########

## pathogenicity
input(qq{21 25606454 25606454 G/C +});
$full_output = `$cmd --sift b --polyphen b`;

# sift
$output = (grep {/ENST00000307301/} (split "\n", $full_output))[0];
$expected = "deleterious";
ok($full_output =~ /$expected/, "sift") or diag("Expected\n$expected\n\nGot\n$output");

# polyphen
$expected = "probably_damaging";
ok($output =~ /$expected/, "polyphen") or diag("Expected\n$expected\n\nGot\n$output");

# humdiv
$full_output = `$cmd --polyphen b --humdiv`;
$expected = 'probably_damaging\(0.998\)';
ok($full_output =~ /$expected/, "polyphen - humdiv") or diag("Expected\n$expected\n\nGot\n$full_output");


## hgvs
$full_output = `$cmd --hgvs`;
$output = (grep {/ENST00000419219/} (split "\n", $full_output))[0];
$expected = 'HGVSc=ENST00000419219.1:c.275C>G';
ok($output =~ /$expected/, "HGVSc") or diag("Expected\n$expected\n\nGot\n$output");
$expected = 'HGVSp=ENSP00000404426.1:p.Ala92Gly';
ok($output =~ /$expected/, "HGVSp") or diag("Expected\n$expected\n\nGot\n$output");


## regulation
input(qq{21 25487468 25487468 A/T +});
$full_output = `$cmd --regulatory`;

# reg feat
$output = (grep {/ENSR00000612061/} (split "\n", $full_output))[0];
$expected = 'regulatory_region_variant';
ok($output =~ /$expected/, "regfeat") or diag("Expected\n$expected\n\nGot\n$output");

# motif
$output = (grep {/MotifFeature/} (split "\n", $full_output))[0];
my %tmp_hash = split /\;|\=/, (split /\t/, $output)[-1];
$expected = {
  MOTIF_POS => 11,
  MOTIF_NAME => 'EBF1:MA0154.2',
  HIGH_INF_POS => 'N',
  MOTIF_SCORE_CHANGE => -0.022,
  STRAND => 1,
  IMPACT => 'MODIFIER',
};
is_deeply(\%tmp_hash, $expected, "motif");

# cell type
$full_output = `$cmd --cell_type HUVEC`;
$output = (grep {/ENSR00000612061/} (split "\n", $full_output))[0];
$expected = 'promoter_flanking_region';
ok($output =~ /$expected/, "cell type") or diag("Expected\n$expected\n\nGot\n$output");

# allele number
input('21 25587758 rs116645811 G A,C . . .');
@lines = split("\n", `$cmd --allele_number | grep -v '#'`);
ok(
  (grep {/T\/M/} grep {/ALLELE_NUM=1/} @lines) &&
  (grep {/T\/R/} grep {/ALLELE_NUM=2/} @lines),
  "allele number"
) or diag "Got\n".join("\n", @lines);

# total length
$output = `$cmd --total_length`;
ok($output =~ /1043\/1199\s+1001\/1062\s+334\/353/, "total length");

# no escape
input('21 25587736 esc C T . . .');
$output = `$cmd --hgvs`;
my $no_esc = `$cmd --hgvs --no_escape`;
ok($no_esc =~ /p\.=/ && $output =~ /p\.\%3D/, "no escape");


## colocated stuff
input(qq{21 25584436 25584436 C/A +});
$output = `$cmd --check_existing --gmaf --maf_1kg --pubmed`;

# colocated ID
$expected = 'rs2282471';
ok($output =~ /$expected/, "colocated ID") or diag("Expected\n$expected\n\nGot\n$output");

# gmaf
$expected = 'GMAF=T:0.1538';
ok($output =~ /$expected/, "GMAF") or diag("Expected\n$expected\n\nGot\n$output");

# population freqs
$expected = 'AFR_MAF=T:0.04';
ok($output =~ /$expected/, "AFR MAF") or diag("Expected\n$expected\n\nGot\n$output");
$expected = 'AMR_MAF=T:0.15';
ok($output =~ /$expected/, "AMR MAF") or diag("Expected\n$expected\n\nGot\n$output");
$expected = 'ASN_MAF=T:0.26';
ok($output =~ /$expected/, "ASN MAF") or diag("Expected\n$expected\n\nGot\n$output");
$expected = 'EUR_MAF=T:0.15';
ok($output =~ /$expected/, "EUR MAF") or diag("Expected\n$expected\n\nGot\n$output");

# pubmed
$expected = 'PUBMED=22272099';
ok($output =~ /$expected/, "pubmed") or diag("Expected\n$expected\n\nGot\n$output");

# check alleles
$output = `$cmd --check_alleles`;
$expected = 'rs2282471';
ok($output !~ /$expected/, "colocated ID") or diag("Expected not to find\n$expected\n\nGot\n$output");

# somatic
input(qq{21 25585742 25585742 G/A +});
$output = `$cmd --check_existing`;
ok($output =~ /COSM162567.*SOMATIC\=1/, "somatic") or diag("Expected\n$expected\n\nGot\n$output");


## frequency filtering
$input = qq{21 25000248 25000248 C/G + test1
21 25000264 25000264 A/G + test2};
input($input);

# filter common
$output = `$cmd --filter_common`;
ok($output !~ /test1/ && $output =~ /test2/, "filter common") or diag("Got\n$output");

# specific filter
$output = `$cmd --check_frequency --freq_pop 1kg_asn --freq_freq 0.04 --freq_gt_lt lt --freq_filter include`;
ok($output !~ /test1/ && $output =~ /test2/, "check frequency 1");

$output = `$cmd --check_frequency --freq_pop 1kg_asn --freq_freq 0.04 --freq_gt_lt lt --freq_filter exclude`;
ok($output =~ /test1/ && $output !~ /test2/, "check frequency 2");

$output = `$cmd --check_frequency --freq_pop 1kg_asn --freq_freq 0.04 --freq_gt_lt gt --freq_filter include`;
ok($output =~ /test1/ && $output !~ /test2/, "check frequency 3");


## external IDs etc
input('21 25587758 rs116645811 G A . . .');
$full_output = `$cmd --ccds --canonical --protein --uniprot --symbol --biotype --tsl --xref_refseq --numbers --domains --variant_class`;
$output = (grep {/ENST00000352957/} (split "\n", $full_output))[0];

# gene symbol
$expected = 'SYMBOL=MRPL39';
ok($output =~ /$expected/, "gene symbol") or diag("Expected\n$expected\n\nGot\n$output");

# hgnc id
$expected = 'HGNC_ID=HGNC:14027';
ok($output =~ /$expected/, "HGNC ID") or diag("Expected\n$expected\n\nGot\n$output");

# biotype
$expected = 'BIOTYPE=protein_coding';
ok($output =~ /$expected/, "biotype") or diag("Expected\n$expected\n\nGot\n$output");

# CCDS
$expected = 'CCDS=CCDS13573.1';
ok($output =~ /$expected/, "CCDS") or diag("Expected\n$expected\n\nGot\n$output");

# protein
$expected = 'ENSP=ENSP00000284967';
ok($output =~ /$expected/, "protein ID") or diag("Expected\n$expected\n\nGot\n$output");

# swissprot
$expected = 'SWISSPROT=Q9NYK5';
ok($output =~ /$expected/, "SWISSPROT") or diag("Expected\n$expected\n\nGot\n$output");

# uniparc
$expected = 'UNIPARC=UPI00001AEE66';
ok($output =~ /$expected/, "UniParc") or diag("Expected\n$expected\n\nGot\n$output");

# refseq xref
$expected = 'RefSeq=NM_017446.3';
ok($output =~ /$expected/, "RefSeq xref") or diag("Expected\n$expected\n\nGot\n$output");

# TSL
$expected = 'TSL=1';
ok($output =~ /$expected/, "TSL") or diag("Expected\n$expected\n\nGot\n$output");

# numbers
$output = (grep {/ENST00000307301/} (split "\n", $full_output))[0];
$expected = 'EXON=10/11';
ok($output =~ /$expected/, "exon number") or diag("Expected\n$expected\n\nGot\n$output");

# protein domains
$expected = 'DOMAINS=Low_complexity_\(Seg\):Seg';
ok($output =~ /$expected/, "protein domains") or diag("Expected\n$expected\n\nGot\n$output");


$expected = 'VARIANT_CLASS=SNV';
ok($output =~ /$expected/, "variant class") or diag("Expected\n$expected\n\nGot\n$output");

# everything
$output = `$cmd --everything | grep -v '#'`;
ok(
  $output =~ /SYMBOL/ &&
  $output =~ /BIOTYPE/ &&
  $output =~ /TSL/ &&
  $output =~ /CCDS/ &&
  $output =~ /ENSP/ &&
  $output =~ /SWISSPROT/ &&
  $output =~ /INTRON/ &&
  $output =~ /EXON/ &&
  $output =~ /HGVSc/ &&
  $output =~ /HGVSp/ &&
  $output =~ /GMAF/ &&
  $output =~ /AFR_MAF/ &&
  $output =~ /AA_MAF/ &&
  $output =~ /CANONICAL/ &&
  $output =~ /CCDS/ &&
  $output =~ /SIFT/ &&
  $output =~ /PolyPhen/ &&
  $output =~ /DOMAINS/ &&
  $output =~ /VARIANT_CLASS/,
  "everything"
) or diag("Got\n$output");

## pick-type options
$full_output = `$cmd --pick`;
@lines = grep {!/^\#/} (split "\n", $full_output);

ok(scalar @lines == 1, "pick - one line");
ok($lines[0] =~ /ENST00000307301/, "pick - correct transcript");

# per gene
$full_output = `$cmd --per_gene`;
@lines = grep {!/^\#/} (split "\n", $full_output);

ok(scalar @lines == 2, "per_gene");

# summary
$output = `$cmd --summary`;
$output =~ m/\s([\w,]*missense[\w,]*)\s/;
my %cons = map {$_ => 1} split(',', $1 || '');
ok($cons{missense_variant} && $cons{intron_variant} && $cons{upstream_gene_variant}, "summary") or diag("Got\n$output");

# most severe
$output = `$cmd --most_severe`;
$expected = '\smissense_variant\s';
ok($output =~ /$expected/, "most severe") or diag("Expected\n$expected\n\nGot\n$output");

input(qq{21 25487468 25487468 A/T +});
$output = `$cmd --regulatory --most_severe`;
$expected = '\sTF_binding_site_variant\s';
ok($output =~ /$expected/, "most severe regulatory") or diag("Expected\n$expected\n\nGot\n$output");

# pick order
input('21 25716535 25716535 A/G +');

my $o1 = `$cmd --pick --pick_order rank,length | grep -v '#'`;
my $o2 = `$cmd --pick | grep -v '#'`;

ok($o1 ne $o2 && $o1 =~ /ENST00000480456/ && $o2 =~ /ENST00000400532/, "pick order") or diag("--pick_order rank,length: $o1\ndefault order: $o2\n");

# pick allele
input('21 25716535 25716535 A/G/T +');
$full_output = `$cmd --pick_allele`;
@lines = grep {!/^\#/} (split "\n", $full_output);
ok(scalar @lines == 2, "pick allele");


## other

# show_cache_info
$output = {map {$_ => 1} split("\n", `$cmd --show_cache_info`)};
$expected = {
  'ClinVar	201410' => 1,
  'sift	sift5.2.2' => 1,
  'regbuild	13.0' => 1,
  'genebuild	2014-07' => 1,
  'assembly	GRCh38.p2' => 1,
  'COSMIC	71' => 1,
  'ESP	20140509' => 1,
  'polyphen	2.2.2' => 1,
  'dbSNP	138' => 1,
  'gencode	GENCODE' => 1,
  'HGMD-PUBLIC	20142' => 1
};
is_deeply($output, $expected, "show cache info");

# check ref
$input = qq{21 25587758 test1 G A . . .
21 25587758 test2 C A . . .};
input($input);

$output = `$cmd --check_ref`;
ok($output =~ /test1/ && $output !~ /test2/, "check ref");

# coding only
input('21 25587758 rs116645811 G A . . .');
$output = `$cmd --coding_only`;
@lines = grep {!/^\#/} split("\n", $output);
ok(scalar @lines == 1, "coding only") or diag("expected 1 line, got ".(scalar @lines));

# no intergenic
input('21 25482183 intergenic1 A G');
$output = `$cmd --no_intergenic`;
ok($output !~ /intergenic/, "no intergenic");

# chr
$input = qq{20 25587758 test1 G A . . .
21 25587758 test2 C A . . .};
input($input);
$output = `$cmd --chr 21`;
ok($output !~ /test1/ && $output =~ /test2/, "chr");

# fork
$input = qq{21      25607440        rs61735760      C       T       .       .       .
21      25606638        rs3989369       A       G       .       .       .
21      25606478        rs75377686      T       C       .       .       .
21      25603925        rs7278284       C       T       .       .       .
21      25603910        rs7278168       C       T       .       .       .
21      25603832        rs116331755     A       G       .       .       .
21      25592893        rs1057885       T       C       .       .       .
21      25592860        rs10576 T       C       .       .       .
21      25592836        rs1135638       G       A       .       .       .
21      25587758        rs116645811     G       A       .       .       .};
input($input);
$output = `$cmd --fork 2 --vcf | grep -v '#' | cut -f 2`;
$expected = "25607440\n25606638\n25606478\n25603925\n25603910\n25603832\n25592893\n25592860\n25592836\n25587758";
ok($output =~ $expected, "fork order") or diag("Expected\n$expected\n\nGot\n$output");

# merged cache
input(qq{21 25606454 25606454 G/C +});
@lines = grep {!/^\#/} split("\n", `$cmd --merged --symbol`);
my $ens = join("\n", grep {/SOURCE=Ensembl/} @lines);
my $ref = join("\n", grep {/SOURCE=RefSeq/} @lines);
ok($ens =~ /ENST00000419219/, "merged cache Ensembl") or diag("Got\n$ens");
ok($ref =~ /NM_080794/, "merged cache RefSeq") or diag("Got\n$ref");
ok($ref =~ /NM_080794.+MRPL39/, "merged cache RefSeq gene symbol") or diag("Got\n$ref");

# stats file
$output = `$bascmd -force -offline -dir_cache $data_path/vep-cache -i $tmpfile -o $tmpfile\.out -cache_version $cver -assembly $ass -species $sp`;
ok(-e $tmpfile.'.out_summary.html', "stats file");
open STATS, $tmpfile.'.out_summary.html';
@lines = <STATS>;
close STATS;
ok((grep {/\<title\>VEP summary\<\/title\>/} @lines), "stats file header");
ok((grep {/\'missense_variant\'\,3/} @lines), "stats missense");

# text stats file
$output = `$bascmd -force -offline -dir_cache $data_path/vep-cache -i $tmpfile -o $tmpfile\.out -cache_version $cver -assembly $ass -species $sp -stats_text`;
ok(-e $tmpfile.'.out_summary.txt', "stats txt file");
open STATS, $tmpfile.'.out_summary.txt';
@lines = <STATS>;
close STATS;
ok((grep {/Overlapped transcripts\s+3/} @lines), "stats txt file overlapped transcripts");

unlink($tmpfile.'.out');
unlink($tmpfile.'.out_summary.txt');
unlink($tmpfile.'.out_summary.html');


## CUSTOM FILES
if(`which tabix` =~ /tabix/) {
  $input = qq{21      25592860        rs10576 T       C       .       .       .
21      25592836        rs1135638       G       A       .       .       .
21      25587758        rs116645811     G       A       .       .       .};
  input($input);
  
  # bed
  $output = `$cmd --custom $data_path/test.bed.gz,testbed,bed,overlap,0`;
  ok($output =~ /bed1/ && $output =~ /bed2/, "custom bed");
  
  # vcf
  $output = `$cmd --custom $data_path/test.vcf.gz,testvcf,vcf,exact,0,ATTR`;
  ok($output =~ /vcf1/ && $output =~ /vcf2/ && $output =~ /vcf3/, "custom vcf");
  ok($output =~ /testvcf_ATTR=test/, "custom vcf INFO");
}
else {
  print STDERR "# tabix not found, skipping custom file tests\n";
}


## PLUGINS
input('21 25587758 rs116645811 G A . . .');
$output = `$cmd --dir_plugins $RealBin\/testdata/ --plugin TestPlugin`;
ok($output =~ /## TestPlugin\s+:\s+Test plugin/, "plugin header") or diag("Got\n$output");
ok($output =~ /TestPlugin=ENST00000567517/, "plugin data") or diag("Got\n$output");


## DATABASE
if(can_connect('ensembldb.ensembl.org')) {
  
  my $dbcmd = "$bascmd -force -database -i $tmpfile -o stdout -db $ver -assembly $ass -species $sp";
  $dbcmd =~ s/\-+cache//;
  
  # ID as input
  input('rs116645811');
  $output = `$dbcmd`;
  ok($output =~ /ENST00000567517/, "DB - ID input") or diag("Got\n$output");
  
  # use cache and database
  input('rs2282471');
  my $tmp_dbcmd = $dbcmd;
  $tmp_dbcmd =~ s/\-+database//;
  $output = `$tmp_dbcmd -force -cache -dir_cache $data_path/vep-cache/ -cache_version $cver -i $tmpfile -o stdout -db $ver -assembly $ass -species $sp -maf_1kg`;
  ok($output =~ /AFR_MAF=T:0.04/, "DB - with cache") or diag("Got\n$output");
  
  # HGVS as input
  input('21:g.25606454G>C');
  $output = `$dbcmd`;
  ok($output =~ /missense_variant/, "DB - HGVS input - genomic") or diag("Got\n$output");
  
  input('ENST00000419219:c.275C>G');
  $output = `$dbcmd`;
  ok($output =~ /missense_variant/, "DB - HGVS input - coding") or diag("Got\n$output");
  
  input('ENSP00000355627:p.Ser206Phe');
  $output = `$dbcmd`;
  ok($output =~ /missense_variant/, "DB - HGVS input - protein") or diag("Got\n$output");
  
  input('AGT:c.803T>C');
  $output = `$dbcmd`;
  ok($output =~ /missense_variant/, "DB - HGVS input - gene symbol") or diag("Got\n$output");
  
  input('LRG_101t1:c.1019T>C');
  $output = `$dbcmd`;
  ok($output =~ /missense_variant/, "DB - HGVS input - LRG") or diag("Got\n$output");
  
  input('NM_080794.3:c.275C>G');
  $output = `$dbcmd --refseq`;
  ok($output =~ /missense_variant/, "DB - HGVS input - RefSeq") or diag("Got\n$output");
  
  # check sv
  $input = '21 25491461 test C G . . .';
  input($input);
  $output = `$dbcmd --check_sv`;
  ok($output =~ /esv2669708/, "DB - check SV") or diag("Expected esv2669708\n\nGot\n$output")
}
else {
  print STDERR "# could not contact ensembldb.ensembl.org, skipping DB tests\n";
}


## BUILD



## INPUTS THAT CAUSE ERROR

# HGVS error
# 6 32191658 . TAGCAGCAGCAGC TAGCAGCAGCAGCAGC,T,TAGCAGCAGC 

#




## FINISHED
###########

unlink($tmpfile) if -e $tmpfile;
unlink("$data_path/vep-cache/$sp/$ver\_$ass/test.fa");
unlink("$data_path/vep-cache/$sp/$ver\_$ass/test.fa.index");
done_testing();


## SUBROUTINES
##############

# writes string $input to file in $tmpfile
sub input {
  my $data = shift;
  open OUT, ">$tmpfile" or die "ERROR: Could not write to temporary file $tmpfile\n";
  print OUT $data;
  close OUT;
}

sub can_connect {
  my $host = shift;
  my $port = shift || 3306;

  eval { use DBI; };
  return 0 if $@;

  my $dsn = sprintf(
    "DBI:mysql(RaiseError=>1):host=%s;port=%i",
    $host,
    $port
  );
  my $dbh = DBI->connect($dsn, 'anonymous', '');

  if($@) {
    return 0;
  }

  else {
    return 1;
  }
}
