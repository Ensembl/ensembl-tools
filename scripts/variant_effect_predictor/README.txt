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
developers list at <ensembl-dev@ebi.ac.uk>.

Questions may also be sent to the Ensembl help desk at
<helpdesk@ensembl.org>
  

Documentation
=============

For a summary of command line flags, run:

perl variant_effect_predictor.pl --help

For full documentation see

http://www.ensembl.org/info/docs/variation/vep/vep_script.html



Changelog
=========

Version 2.1
-----------

- ability to use local file cache in place of or alongside connecting to an
  Ensembl database

- significant improvements to speed of script

- whole-genome mode now default (no disadvantage for smaller datasets)

- improved status output with progress bars

- regulatory region consequences now reinstated and improved

- modification to output file - Transcript column is now Feature, and is
  followed by a Feature_type column
  
- full documentation now online


Version 2.0
-----------

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
