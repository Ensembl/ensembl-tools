use strict;
use warnings;
BEGIN { $| = 1;
	use Test::More;
	#plan tests => 6;
}

use FindBin qw($RealBin);
use Data::Dumper;

# configure script
my $script = $RealBin.'/../filter_vep.pl';
my $perl   = '/usr/bin/env perl';
my $cmd = "$perl $script";

my $txt = "$RealBin\/testdata/filter.txt";
my $vcf = "$RealBin\/testdata/filter.vcf";

my ($output, @lines);

my $opcmd = "$cmd -i $RealBin\/testdata/filter.txt";

# ontology
if(can_connect('ensembldb.ensembl.org')) {
  @lines = grep {!/^\#/} split("\n", `$opcmd -ontology -f "Consequence is coding_sequence_variant"`);
  ok(
    (grep {/missense_variant/} @lines) &&
    (grep {/synonymous_variant/} @lines),
    "ontology"
  );
}
else {
  print STDERR "# could not contact ensembldb.ensembl.org, skipping DB tests\n";
}

done_testing();

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
