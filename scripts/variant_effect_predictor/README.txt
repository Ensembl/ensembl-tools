############################
#                          #
# Variant Effect Predictor #
#                          #
############################

Copyright (c) 1999-2011 The European Bioinformatics Institute and
Genome Research Limited.  All rights reserved.

This software is distributed under a modified Apache license.
For license details, please see

http://www.ensembl.org/info/about/code_licence.html

Please email comments or questions to the public Ensembl
developers list at <dev@ensembl.org>.

Questions may also be sent to the Ensembl help desk at
<helpdesk@ensembl.org>
  
Quickstart
==========

Install API and cache files, run in offline mode:

perl INSTALL.pl
perl variant_effect_predictor.pl --offline


Documentation
=============

For a summary of command line flags, run:

perl variant_effect_predictor.pl --help

For full documentation see

http://www.ensembl.org/info/docs/variation/vep/vep_script.html



Changelog
=========

New in version 2.5 (May 2012)
-----------------------------

- SIFT and PolyPhen predictions now available for RefSeq transcripts

- retrieve cell type-specific regulatory consequences

- consequences can be retrieved based on a single individual's genotype in
  a VCF input file

- find overlapping structural variants

- Condel support removed from main script and moved to a plugin

New in version 2.4 (February 2012)
----------------------------------
- offline mode and new installer script make it easy to use the VEP without
  the usual dependencies

- output columns configurable using the --fields flag

- VCF output support expanded, can now carry all fields

- output affected exon and intron numbers with --numbers

- output overlapping protein domains using --domains

- enhanced support for LRGs

- plugins now work on variants called as intergenic


New in version 2.3 (December 2011)
----------------------------------

- Add custom annotations from tabix-indexed files (BED, GFF, GTF, VCF, bigWig)

- Add new functionality to the VEP with user-written plugins

- Filter input on consequence type


Version 2.2 (September 2011)
----------------------------

- SIFT, PolyPhen and Condel predictions and regulatory features now accessible
  from the cache

- Support for calling consequences against RefSeq transcripts

- Variant identifiers (e.g. dbSNP rsIDs) and HGVS notations supported as input
  format

- Variants can now be filtered by frequency in HapMap and 1000 genomes
  populations

- Script can be used to convert files between formats (Ensembl/VCF/Pileup/HGVS
  to Ensembl/VCF/Pileup)

- Large amount of code moved to API modules to ensure consistency between web
  and script VEP
  
- Memory usage optimisations

- VEP script moved to ensembl-tools CVS module

- Added --canonical, --per_gene and --no_intergenic options


Version 2.1 (June 2011)
-----------------------

- ability to use local file cache in place of or alongside connecting to an
  Ensembl database

- significant improvements to speed of script

- whole-genome mode now default (no disadvantage for smaller datasets)

- improved status output with progress bars

- regulatory region consequences now reinstated and improved

- modification to output file - Transcript column is now Feature, and is
  followed by a Feature_type column
  
- full documentation now online


Version 2.0 (April 2011)
------------------------

Version 2.0 of the Variant Effect Predictor script (VEP) constitutes a complete
overhaul of both the script and the API behind it. It requires at least version
62 of the Ensembl API to function. Here follows a summary of the changes:

- support for SIFT, PolyPhen and Condel non-synonymous predictions in human

- per-allele and compound consequence types

- support for Sequence Ontology (SO) and NCBI consequence terms

- modified output format
  - support for new output fields in Extra column
  - header section containing information on database and software versions
  - codon change shown in output
  - CDS position shown in output
  - option to output Ensembl protein identifiers
  - option to output HGVS nomenclature for variants
  
- support for gzipped input files
  
- enhanced configuration options, including the ability to read configuration
  from a file

- verbose output now much more useful

- whole-genome mode now more stable

- finding existing co-located variations now ~5x faster
