README

Extract sequences 5' flanking or at 5' end of short reads mapped to a reference genome

DEPENDENCIES

Python 2.7
pyfasta
numpy

USAGE

python ExtractFlankingSequence.py sam_file   reference_sequence.fasta   #_bp_to_extract   output_file_name

if # bp to extract < 0 extracts 5' flanking sequence

if # bp to extract > 0 extracts 5' mapped region

e.g.

python ExtractFlankingSequence.py hits.sam dmel-4-chromosome-r5.33.fasta -20 flanks.fasta


David L. Stern
Janelia Farm Research Campus
6 Dec. 2012


 * Copyright 2012 Howard Hughes Medical Institute.
 * All rights reserved.
 * Use is subject to Janelia Farm Research Campus Software Copyright 1.1
 * license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html ).
