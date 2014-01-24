##############################
#                            #
#   The Region Report Tool   #
#                            #
##############################

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

A script for sampling a given set of chromosomal regions, producing a simple
summary of the features within those regions.
Operates in batch and single region mode, and can serialize into GFF3 and 
textual format.
 
Documentation
=============

The Region Report tool produces summaries of small sections of EnsEMBL data in
readable form. It writes a list of all features in a chosen genomic slice given
the following parameters:

	Species
	Chromosome
	Start position
	End position
	Strand
	
In addition you may choose which features to print out from this selection:

	Gene
	Transcript
	Translation
	Variation (EnsEMBL Variation)
	Structural Variant (EnsEMBL Variation)
	Constrained Element (EnsEMBL Compara)
	Regulatory Feature (EnsEMBL Funcgen)
	FASTA sequence dump
	
At present the Region Report tool can produce both human-readable text and 
GFF3 formats.
	
# Note on large region reports #

You may feel tempted to use this tool to dump large amounts of GFF, but we ask
that you refrain from doing so if using EnsEMBL database servers. 
In particular, Variation features have very high feature density. As a matter 
of courtesy, we ask that you aim to keep your selected regions under 200kb. 
Web users will experience limits depending on their choice of features.

# Dependencies #

The client must have the installed the following:

ensembl-core
ensembl-funcgen
ensembl-variation
ensembl-compara

Use only stable releases and their corresponding data to avoid problems.

Also ensure that the modules folder from each bundle is added to your Perl
library path ($PERL5LIB on linux).

# Output destination #

By default the Region Report tool prints to the console, and can be redirected 
or piped in the normal way. Named file output can also be declared by using the 
--output=filename command.

# Using a different database server #
	
Should you wish to use a server other than EnsEMBL's, further command line options
can be specified, e.g.

region_report.pl --host=local-mirror-server --user=me --password=go

### Operating in batch mode ###

If you have many different regions to explore at once, they can be run in a 
batch. Create a text file containing the regions of interest, one per line:

	5:100000-110000
	7:1-150000
	x:17000000-17100000
	y:10000..20000
	
The first number or letter is the chromosome, followed by the start and end of 
the region of interest. The region can be specified as one would when using the
EnsEMBL genome browser.

Save your file and use it for input, using the --input=batch_file command. The
output will be combined into one large single report.