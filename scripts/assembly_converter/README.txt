############################
#                          #
# Assembly   Converter     #
#                          #
############################

Copyright (c) 1999-2011 The European Bioinformatics Institute and
Genome Research Limited.  All rights reserved.

This software is distributed under a modified Apache license.
For license details, please see

http://www.ensembl.org/info/about/code_licence.html

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>
  
Summary
=======

A program to map slices from old assemblies to the latest assembly.


Documentation
=============

For a summary of command line flags, run:

  perl AssemblyMapper.pl --help

Usage:

  AssemblyMapper.pl --species=species --file=filename

    --species / -s  Name of species.

    --genomes / -g  Automatically sets DB params for e!Genomes

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

  AssemblyMapper.pl -s human -f assemblymapper.in

An example input file (assemblymapper.in) is provided with this script.
