use strict;
use warnings;
BEGIN { $| = 1;
	use Test::More;
	#plan tests => 6;
}

use FindBin qw($Bin);
use Data::Dumper;
use File::Copy;
use lib $Bin.'/../';
use Bio::EnsEMBL::Variation::Utils::VEP;

# find where ensembl-variation is installed
my $mod_path = 'Bio/EnsEMBL/Variation/Utils/VEP.pm';
my $var_path = $INC{$mod_path};
die("ERROR: Could not find path to ensembl-variation modules\n") unless $var_path;
$var_path =~ s/$mod_path//;
my $data_path = $var_path.'/t/testdata/';
die("ERROR: Could not find test data path $data_path\n") unless -d $data_path;

# configure script
open CONF, "$Bin\/test.conf" or die "ERROR: Could not read from conf file $Bin\/test.conf\n";
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

my $script = $Bin.'/../convert_cache.pl';
my $perl   = '/usr/bin/env perl';
my $bascmd = "$perl $script";
my $cmd    = "$bascmd -version $ver\_$ass -species $sp -dir $data_path\/vep-cache/ -quiet";

ok(-e $script, "script exists");

# check for tabix and bgzip
unless(`which tabix` =~ /tabix/ && `which bgzip` =~ /bgzip/) {
  print STDERR "# tabix and/or bgzip not found, skipping convert_cache tests\n";
  done_testing();
  exit(0);
}

# backup info.txt
copy("$data_path\/vep-cache/$sp/$ver\_$ass/info.txt", "$data_path\/vep-cache/$sp/$ver\_$ass/info.txt.bak");

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


# restore cache to previous state
unlink("$data_path\/vep-cache/$sp/$ver\_$ass/info.txt");
move("$data_path\/vep-cache/$sp/$ver\_$ass/info.txt.bak", "$data_path\/vep-cache/$sp/$ver\_$ass/info.txt");
unlink("$data_path\/vep-cache/$sp/$ver\_$ass/21/all_vars.gz");
unlink("$data_path\/vep-cache/$sp/$ver\_$ass/21/all_vars.gz.tbi");

done_testing();
