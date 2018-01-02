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


# This is a program that takes a file with a list of stable IDs (not
# exon stable IDs), and outputs a comma-separated list of the history of
# each of these stable IDs.  The history ends with either the current
# release or at the point when the stable ID was retired.

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use IO::File;
use Getopt::Long;
use DBI qw(:sql_types);

my ( $filename, $species, $host, $user, $port, $pass );
my $help = '';

if ( !GetOptions( 'file|f=s'    => \$filename,
                  'species|s=s' => \$species,
                  'host=s'      => \$host,
                  'user=s'      => \$user,
                  'port=s'      => \$port,
                  'pass=s'      => \$pass,
                  'help|h!'     => \$help )
     || !defined($species)
     || $help )
{
  print <<END_USAGE;

Usage:
  $0 --species=species [ --file=filename ]

  $0 --help


    --species / -s  Name of species.

    --file    / -f  (Optional) Name of file containing a list of stable
                    IDs.  The default is to read the list from standard
                    input.

    --host          (Optional) Database host (defaults to public ensembldb)

    --user          (Optional) User for the database connection

    --port          (Optional) Port number for the database connection

    --pass          (Optional) Password for the database connection

    --help    / -h  To see this text.

Example usage:

  $0 -s mouse -f idlist.txt

  # Same as the above:
  $0 -s mouse < idlist.txt

  sort -u longlist.txt | $0 -s human

END_USAGE

  exit(1);
} ## end if ( !GetOptions( 'file|f=s'...))

$filename ||= '-';

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db( -host => $host || 'ensembldb.ensembl.org',
                                  -user => $user || 'anonymous',
                                  -port => $port || 3306,
                                  -pass => $pass || '' );

my $adaptor = $registry->get_DBAdaptor( $species, 'Core' );

my $in = IO::File->new($filename);

if ( !defined($in) ) {
  die( sprintf( "Could not open file '%s' for reading", $filename ) );
}

# We could do what we want to do with the API, but this is simpler and
# quicker, at the moment.  As always, when using plain SQL against our
# databases, the user should not be surprised to see the code break when
# we update the schema...
my $statement = q(
SELECT  old_version,
        old_release,
        new_stable_id, new_version,
        new_release,
        score
FROM    stable_id_event
  JOIN  mapping_session USING (mapping_session_id)
WHERE   old_stable_id = ?
ORDER BY old_version ASC, CAST(new_release AS UNSIGNED)
);

my $sth = $adaptor->dbc()->db_handle()->prepare($statement);

while ( my $line = $in->getline() ) {

  # Strip off any comment (from '#' to the end of the line).
  $line =~ s/\s*#.*$//;

  # Strip off any padding white spaces, new line chars and carriage returns
  $line =~ s/^\s*|\s*$//g;

  # Skip lines containing only whitespace.
  next unless $line;

  # if comma or space used as a delimiter, split the string into multiple ids
  foreach my $stable_id (grep $_, split /[\s\,]+/, $line) {

    print("Old stable ID, New stable ID, Release, Mapping score\n");

    $sth->bind_param( 1, $stable_id, SQL_VARCHAR );

    $sth->execute();

    my ( $version, $release, $new_stable_id, $new_version, $new_release, $score );

    $sth->bind_columns( \( $version, $release, $new_stable_id, $new_version, $new_release, $score ) );

    while ( $sth->fetch() ) {

      if ( defined($new_stable_id) ) {

        printf( "%s.%s, %s.%s, %s, %s\n", $stable_id, $version, $new_stable_id, $new_version, $new_release, $score );

      } elsif ( !defined($new_stable_id) ) {

        printf( "%s.%s, <retired>, %s, %s\n", $stable_id, $version, $new_release, $score );

      }
    }

    print("\n");
  }

} ## end while ( my $stable_id = $in...)
