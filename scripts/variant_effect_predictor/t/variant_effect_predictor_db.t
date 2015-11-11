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
  ok($output =~ /AFR_MAF=T:0.0356/, "DB - with cache") or diag("Got\n$output");
  
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
