"""
Extract sequences 5' flanking or at 5' end of short reads mapped to a reference genome

DEPENDENCIES

Python 2.7
pyfasta
numpy

USAGE

python ExtractFlankingSequence.py <sam file> <reference sequence.fasta> <# bp to extract> <output file name>
# bp to extract < 0 extracts 5' flanking sequence
# bp to extract > 0 extracts 5' mapped region

e.g.

python ExtractFlankingSequence.py hits.sam dmel-4-chromosome-r5.33.fasta 20 flanks.fasta


David L. Stern
Janelia Farm Research Campus
6 Dec. 2012

/*
 * Copyright 2012 Howard Hughes Medical Institute.
 * All rights reserved.
 * Use is subject to Janelia Farm Research Campus Software Copyright 1.1
 * license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html ).
 */

"""
from pyfasta import Fasta
import numpy as np
import sys
import getopt
import string
import re

def main():
        #parse command line options
        try:
                opts, arg = getopt.getopt(sys.argv[1:],"h", ["help"])
        except getopt.error, msg:
                print msg
                print "for help use --help"
                sys.exit(2)
        # process options
        for o, a in opts:
                if o in ("-h", "--help"):
                        print __doc__
                        sys.exit(0)
        if len(arg) < 4:
                print "Usage: python ExtractFlankingSequence.py <sam file> <reference sequence.fasta> <# bp to extract> <output file name>\n# bp to extract < 0 extracts 5' flanking sequence\n# bp to extract > 0 extracts 5' mapped region"                
                sys.exit(0)
        #process arguments

        SamFile = arg[0]
        RefGenome = arg[1]
        NumBP = int(arg[2])
        OutFile = arg[3]
        
        print "sam file = %s" %(SamFile)
        print "reference genome = %s" %(RefGenome)
        print "Num bp to extract = %s" %(NumBP)
        print "Output File = %s" %(OutFile)
        
        f = Fasta(RefGenome)
        
        sf = open(SamFile,'r')
        positions = []
        sequences = []
        i = 0
        for line in sf:
                if line[0] == '@':
                        continue
                elements = line.split()
                if elements[2] == '*':#discard reads that aren't mapped
                        continue
                if elements[4] < 20:#discard mapped reads with low quality
                        continue
                else:
                        i = i+1
                        if i % 10000 == 0:                       
                                print "%s lines done" %(i)
                        if int(elements[1]) == 0: # check read direction
                                forward = True
                        elif int(elements[1]) == 16:
                                forward = False                                                        
                        bp = int(elements[3])
                        ch = str(elements[2])
                        StartDiff,insertion = idFromCIGAR(elements[5])
                        recalbp = bp - StartDiff#recalibrate start if initial positions were masked
                        if recalbp + NumBP < 1: #discard reads that map too close to either end of the chromosome
                                continue
                        elif recalbp - NumBP > len(f[ch]) and forward == False:
                                continue                 
                        
                        if int(NumBP) < 0:#grab flanking sequence
                                if forward == True:
                                        seqn = f[ch][recalbp+NumBP - 1:recalbp - 1]#these arrays count from 0
                                elif forward == False:
                                        seqn = f[ch][recalbp + insertion + len(elements[9]) -1: recalbp + insertion + len(elements[9]) - NumBP -1]
                                        seqn = rc(seqn) #reverse complement sequence
                        elif int(NumBP) > 0:#grab internal sequence
                                if forward == True:
                                        seqn = f[ch][recalbp-1:recalbp + NumBP - 1]
                                elif forward == False:
                                        seqn = f[ch][recalbp + insertion + len(elements[9]) - NumBP - 1 : recalbp + insertion + len(elements[9]) - 1]
                                        seqn = rc(seqn) #reverse complement sequence
                        
                        # save seqn in array                                
                        positions.append("ch" + ch + "_" + str(bp))
                        sequences.append(seqn)
                                
        sf.close()
                
        of = open(OutFile,'w')
        i = 0
        for seq in sequences:
                of.write(">" + positions[i] +"\n")
                of.write(seq + "\n")
                i= i+1
        of.close()        
          
def rc(dna):

        complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
        rcseq = dna.translate(complements)[::-1]
        return rcseq          

def idFromCIGAR(CIGAR):
        position = 0
        Insertion = 0
        StartDiff = 0
        parsedCIGAR = re.findall(r"[^\W\d_]+|\d+",CIGAR)
        #if sequence masked before first mapped position, then return this # bp
        if parsedCIGAR[1] == 'S':
                StartDiff = int(parsedCIGAR[0])
        for element in CIGAR:
                if element == 'D':
                        Insertion = Insertion + int(CIGAR[position - 1])
                elif element == 'I':
                        Insertion = Insertion - int(CIGAR[position - 1])
                position = position + 1
        return StartDiff,Insertion
  
if __name__ == "__main__":
        sys.exit(main())