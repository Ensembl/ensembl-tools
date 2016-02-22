use strict;
use warnings;
BEGIN { $| = 1;
	use Test::More;
	#plan tests => 6;
}

use FindBin qw($RealBin);
use Data::Dumper;
use File::Copy;
use lib $RealBin.'/../';
use Bio::EnsEMBL::Variation::Utils::VEP;

# find where ensembl-variation is installed
my $mod_path = 'Bio/EnsEMBL/Variation/Utils/VEP.pm';
my $var_path = $INC{$mod_path};
die("ERROR: Could not find path to ensembl-variation modules\n") unless $var_path;
$var_path =~ s/$mod_path//;
my $data_path = $var_path.'/t/testdata/';
die("ERROR: Could not find test data path $data_path\n") unless -d $data_path;

# configure script
open CONF, "$RealBin\/test.conf" or die "ERROR: Could not read from conf file $RealBin\/test.conf\n";
my $config;
while(<CONF>) {
  chomp;
  my ($k, $v) = split;
  $config->{$k} = $v;
}
close CONF;

my $ver = $config->{cache_version};
my $ass = $config->{assembly};
my $sp  = $config->{species};

my $script = $RealBin.'/../convert_cache.pl';
my $perl   = '/usr/bin/env perl';
my $bascmd = "$perl $script";
my $cmd    = "$bascmd -version $ver\_$ass -species $sp -dir $data_path\/vep-cache/ -quiet";

# backup info.txt
copy("$data_path\/vep-cache/$sp/$ver\_$ass/info.txt", "$data_path\/vep-cache/$sp/$ver\_$ass/info.txt.bak") or die "ERROR: Failed to backup info.txt\n";

ok(-e $script, "script exists");

# check for tabix and bgzip
unless(`which tabix` =~ /tabix/ && `which bgzip` =~ /bgzip/) {
  print STDERR "# tabix and/or bgzip not found, skipping convert_cache tests\n";
  finish_script();
}

# run script
ok(!system($cmd), "run script");

# check info.txt
open INFO, "$data_path\/vep-cache/$sp/$ver\_$ass/info.txt";
my @lines = <INFO>;
close INFO;

ok((grep {/var_type\s+tabix/} @lines), "info.txt - var_type tabix");

my ($cols) = grep {/variation_cols/} @lines;
ok($cols, "info.txt - variation_cols");
ok($cols =~ /chr/, "info.txt - first col chr");

my %col_nums;
my @split = split(',', $cols);
$col_nums{$split[$_]} = $_ for 0..$#split;

# check files exist
ok(-e "$data_path\/vep-cache/$sp/$ver\_$ass/21/all_vars.gz", "all_vars.gz");
ok(-e "$data_path\/vep-cache/$sp/$ver\_$ass/21/all_vars.gz.tbi", "all_vars.gz.tbi");

open VARS, "gzip -dc $data_path\/vep-cache/$sp/$ver\_$ass/21/all_vars.gz |";
my (%col_counts, $pos_non_ints);

while(<VARS>) {
  chomp;
  $col_counts{tr/\t/\t/}++;
  @split = split;
  
  $pos_non_ints++ unless $split[$col_nums{start}] =~ /^[0-9]+$/;
}

close VARS;

# find max col count
my $max = (sort {$a <=> $b} keys %col_counts)[-1] + 1;
ok($max == scalar keys %col_nums, "all_vars.gz - column number") or diag("Max $max\nKeys ".(scalar keys %col_nums));
ok(!$pos_non_ints, "all_vars.gz - start int");

eval q{ use Sereal::Decoder; };
if($@) {
  print STDERR "# Sereal module not installed, skipping Sereal tests\n";
  finish_script();
}

unlink("$data_path\/vep-cache/$sp/$ver\_$ass/info.txt");
copy("$data_path\/vep-cache/$sp/$ver\_$ass/info.txt.bak", "$data_path\/vep-cache/$sp/$ver\_$ass/info.txt");
unlink("$data_path\/vep-cache/$sp/$ver\_$ass/21/all_vars.gz");
unlink("$data_path\/vep-cache/$sp/$ver\_$ass/21/all_vars.gz.tbi");

ok(!system($cmd.' --sereal'), "run script with --sereal");

# check info.txt
open INFO, "$data_path\/vep-cache/$sp/$ver\_$ass/info.txt";
@lines = <INFO>;
close INFO;

ok((grep {/serialiser_type\s+sereal/} @lines), "info.txt - serialiser_type sereal");

ok(-e "$data_path\/vep-cache/$sp/$ver\_$ass/21/25000001-26000000.sereal", "transcript file exists");
ok(-e "$data_path\/vep-cache/$sp/$ver\_$ass/21/25000001-26000000_reg.sereal", "regfeat file exists");

my $decoder = Sereal::Decoder->new();

# test transcript file
open IN, "$data_path\/vep-cache/$sp/$ver\_$ass/21/25000001-26000000.sereal";
my $obj = $decoder->decode(join('', <IN>));
close IN;

ok($obj, "parsed transcript file");
ok($obj->{'21'}, "transcript hash index");
ok($obj->{'21'}->[0]->isa('Bio::EnsEMBL::Transcript'), "transcript isa");
ok($obj->{'21'}->[0]->{stable_id} eq 'ENST00000441009', "transcript stable_id");

# test regfeat file
open IN, "$data_path\/vep-cache/$sp/$ver\_$ass/21/25000001-26000000_reg.sereal";
$obj = $decoder->decode(join('', <IN>));
close IN;

ok($obj, "parsed regfeat file");
ok($obj->{'21'}, "regfeat hash index 1");
ok($obj->{'21'}->{RegulatoryFeature}, "regfeat hash index 2");
ok($obj->{'21'}->{MotifFeature}, "regfeat hash index 3");
ok($obj->{'21'}->{RegulatoryFeature}->[0]->isa('Bio::EnsEMBL::Funcgen::RegulatoryFeature'), "regfeat isa");
ok($obj->{'21'}->{RegulatoryFeature}->[0]->{stable_id} eq 'ENSR00001565774', "regfeat stable_id");

finish_script();


sub finish_script {
  # restore cache to previous state
  unlink("$data_path\/vep-cache/$sp/$ver\_$ass/info.txt");
  move("$data_path\/vep-cache/$sp/$ver\_$ass/info.txt.bak", "$data_path\/vep-cache/$sp/$ver\_$ass/info.txt");
  unlink("$data_path\/vep-cache/$sp/$ver\_$ass/21/all_vars.gz");
  unlink("$data_path\/vep-cache/$sp/$ver\_$ass/21/all_vars.gz.tbi");
  unlink("$data_path\/vep-cache/$sp/$ver\_$ass/21/25000001-26000000.sereal");
  unlink("$data_path\/vep-cache/$sp/$ver\_$ass/21/25000001-26000000_reg.sereal");

  done_testing();

  exit(0);
}
