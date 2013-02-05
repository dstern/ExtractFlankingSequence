README

ExtractFlankingSequence : Grab reference genome sequence 5' flanking or at 5' end of 
                          short reads mapped to a reference genome

Extract1stNBasesFromMappedReads.py: Grab 1st N bp from reads mapped to reference genome

CompareExtractedFiles.py: Compare results from the two scripts listed above

DEPENDENCIES

Python 2.7
pyfasta
numpy

USAGE

python ExtractFlankingSequence.py sam_file   reference_sequence.fasta   #_bp_to_extract   output_file_name
python Extract1stNBasesFromMappedReads.py <sam file> <reference sequence.fasta> <# bp to extract> <output file name> <optional sequence prefix>
python CompareExtracts.py <extract1.fa> <extract2.fa> <output file name>


e.g.

python ExtractFlankingSequence.py hits.sam dmel-4-chromosome-r5.33.fasta -20 flanks.fasta
python ExtractFlankingSequence.py hits.sam dmel.genome.fasta 20 out.fasta GATGGCAT
python ExtractFlankingSequence.py flanks.fasta out.fasta compare.fa


David L. Stern
Janelia Farm Research Campus
5 February 2013


 * Copyright 2013 Howard Hughes Medical Institute.
 * All rights reserved.
 * Use is subject to Janelia Farm Research Campus Software Copyright 1.1
 * license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html ).
