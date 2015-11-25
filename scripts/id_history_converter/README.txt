############################
#                          #
# ID History Converter     #
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

This is a program that takes a file with a list of stable IDs (for genes, 
transcripts and/or translations), and outputs a comma-separated list of 
the history of  each of these stable IDs. The history ends with either 
the current release or at the point when the stable ID was retired.


Documentation
=============

For a summary of command line flags, run:

perl IDmapper.pl --help

Usage:
  IDmapper.pl --species=species [ --file=filename ]

  IDmapper.pl --help


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

  IDmapper.pl -s mouse -f idlist.txt

  # Same as the above:
  IDmapper.pl -s mouse < idlist.txt

  sort -u longlist.txt | IDmapper.pl -s human

An example input file (idmapper.in) is provided with this script.
