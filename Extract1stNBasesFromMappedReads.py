#!/usr/bin/env python
"""
Extract 1st N bases from reads mapped to a reference genome

DEPENDENCIES

Python 2.7
pyfasta
numpy

USAGE

python Extract1stNBasesFromMappedReads.py <sam file> <reference sequence.fasta> <# bp to extract> <output file name> <optional sequence prefix>

e.g.

python ExtractFlankingSequence.py hits.sam dmel.genome.fasta 20 out.fasta GATGGCAT


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
        if len(arg) < 3:
                print "\nUsage: python Extract1stNBasesFromMappedReads.py <sam file> <reference sequence.fasta> <# bp to extract> <output file name> <optional sequence prefix>\n"                
                sys.exit(0)
        #process arguments

        SamFile = arg[0]
        RefGenome = arg[1]
        NumBP = int(arg[2])
        OutFile = arg[3]
        
        if len(arg) == 5:
                prefix = arg[4]
        else:
                prefix = ''
                
        
        print "sam file = %s" %(SamFile)
        print "Num bp to extract = %s" %(NumBP)
        print "Output File = %s" %(OutFile)
        
        f = Fasta(RefGenome,key_fn=lambda key: key.split()[0])

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
                if int(elements[4]) < 20:#discard mapped reads with low quality
                        continue
                if len(elements[9]) < NumBP:
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
                        elif int(NumBP) > 0:#grab internal sequence
                                if forward == True:
                                        seqn = str(elements[9][0:NumBP])
                                elif forward == False:
                                        seqn = str(elements[9][-NumBP:])
                                        seqn = rc(seqn) #reverse complement sequence
                        
                        # save seqn in array                                
                        positions.append("ch" + ch + "_" + str(bp))
                        sequences.append(prefix + seqn)
                                
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
        for element in parsedCIGAR:
                if element == 'D':
                        Insertion = Insertion + int(parsedCIGAR[position - 1])
                elif element == 'I':
                        Insertion = Insertion - int(parsedCIGAR[position - 1])
                position = position + 1
        return StartDiff,Insertion
  
if __name__ == "__main__":
        sys.exit(main())