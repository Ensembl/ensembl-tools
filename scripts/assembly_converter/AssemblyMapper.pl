#!/usr/bin/env perl

# A program to map slices from old assemblies to the latest assembly.

use strict;
use warnings;

use IO::File;
use Getopt::Long;

use Bio::EnsEMBL::Registry;

my ( $filename, $species );
my $help = '';

if ( !GetOptions( 'file|f=s'    => \$filename,
                  'species|s=s' => \$species,
                  'help|h!'     => \$help )
     || !( defined($filename) && defined($species) )
     || $help )
{
  print <<END_USAGE;

Usage:
  $0 --species=species --file=filename

  $0 --help


    --species / -s  Name of species.  Mappings are currently only
                    available for mouse and human.

    --file / -f     Name of file containing a list of slices to map to
                    the most recent assembly.  The format of the data
                    is the same as the output of the name() method on a
                    Slice object:

                      coord_system:version:seq_region_name:start:end:strand

                    For example:

                      chromosome:NCBI36:X:1:10000:1

                    NB:  Mappings are only available for chromosomes to
                    the latest assembly, and only from old assemblies to
                    the latest assembly, not between old assemblies or
                    *from* the latest assembly.

                    If the strand is missing, the positive ("1") strand
                    will be used.

                    If the start is missing, it is taken to be "1".  If
                    the end is missing, it is taken to be the end of the
                    seq_region.

    --help    / -h  To see this text.

Example usage:

  $0 -s mouse -f slices.txt

END_USAGE

  exit(1);
} ## end if ( !GetOptions( 'file|f=s'...))

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db( '-host' => 'ensembldb.ensembl.org',
                                  '-user' => 'anonymous' );

my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );

my $in = IO::File->new($filename);

if ( !defined($in) ) {
  die( sprintf( "Could not open file '%s' for reading", $filename ) );
}

while ( my $line = $in->getline() ) {
  chomp($line);

  # Strip off any comment (from '#' to the end of the line).
  $line =~ s/\s*#.*$//;

  # Skip lines containing only whitespace.
  if ( $line =~ /^\s*$/ ) { next }

  # We're assuming that the line will be in the same format as what's
  # outputted by the name() method for a Slice object.
  if ( $line !~ /^(\w+)?:(\w+)?:(\w+):(\d+)?:(\d+)?:(-?\d+)?$/ ) {
    printf( "Malformed line:\n%s\n", $line );
    next;
  }

  my ( $old_cs_name, $old_version, $old_sr_name,
       $old_start,   $old_end,     $old_strand
  ) = ( $1, $2, $3, $4, $5, $6 );

  # Get a slice for the old region (the region in the input file).
  my $old_slice =
    $slice_adaptor->fetch_by_region(
                                $old_cs_name, $old_sr_name, $old_start,
                                $old_end,     $old_strand,  $old_version
    );

  # Complete possibly missing info.
  $old_cs_name ||= $old_slice->coord_system_name();
  $old_sr_name ||= $old_slice->seq_region_name();
  $old_start   ||= $old_slice->start();
  $old_end     ||= $old_slice->end();
  $old_strand  ||= $old_slice->strand();
  $old_version ||= $old_slice->coord_system()->version();

  printf( "# %s\n", $old_slice->name() );

  # Project the old slice to the current assembly and display
  # information about each resulting segment.
  foreach my $segment ( @{ $old_slice->project('chromosome') } ) {
    # We display the old slice info followed by a comma and then the new
    # slice (segment) info.
    printf( "%s:%s:%s:%d:%d:%d,%s\n",
            $old_cs_name,
            $old_version,
            $old_sr_name,
            $old_start + $segment->from_start() - 1,
            $old_start + $segment->from_end() - 1,
            $old_strand,
            $segment->to_Slice()->name() );
  }
  print("\n");

} ## end while ( my $line = $in->getline...)
