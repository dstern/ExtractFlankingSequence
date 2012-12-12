#!/usr/bin/env python
from __future__ import division
"""
Compare two extraction files

DEPENDENCIES

Python 2.7
numpy

USAGE

python CompareExtracts.py <extract1.fa> <extract2.fa> <output file name>

e.g.

python ExtractFlankingSequence.py a.fa b.fa out.fa


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
                print "Usage: python CompareExtracts.py <extract1.fa> <extract2.fa> <output file name>"                
                sys.exit(0)
        #process arguments

        afile= arg[0]
        bfile = arg[1]
        OutFile = arg[2]
        
        print "first file = %s" %(afile)
        print "second file = %s" %(bfile)
        print "Output File = %s" %(OutFile)
		
	f = open(afile)
	g = open(bfile)
	a = f.readlines()
	b = g.readlines()
	f.close()
	g.close()
	
	#get seqn length
	seqlen = len(a[1].strip())
	
	#make arrays for sequence motides
	c = np.zeros((len(a)/2,seqlen),dtype=np.character)
	d = np.zeros((len(b)/2,seqlen),dtype=np.character)
	
	#grab only sequences, ignore identifier lines
	j =  a[1:len(a):2]
	k =  b[1:len(b):2]
	
	#dump data in numpy arrays
	x = 0
	for i in j:
		c[x,:] = list(i[0:seqlen])
		x+=1
	
	y = 0
	for i in k:
		d[y,:] = list(i[0:seqlen])
		y+=1	
	
	#compare numpy arrays as logical
	z = (c == d)
	matches = sum(z)
	AllProp = matches / len(z)	
	
	#now determine matches for A, T, C, G

	AMatches = np.zeros(seqlen)
	TMatches = np.zeros(seqlen)
	CMatches = np.zeros(seqlen)
	GMatches = np.zeros(seqlen)
	
	#can't figure out how to calculate on all columns simultaneously, so grab one column at a time
	for i in range(seqlen):
		bpCol = c[:,i]
		boolCol = z[:,i]
		TMatches[i] = sum(boolCol[np.where(bpCol == 'T')])
		AMatches[i] = sum(boolCol[np.where(bpCol == 'A')])
		CMatches[i] = sum(boolCol[np.where(bpCol == 'C')])
		GMatches[i] = sum(boolCol[np.where(bpCol == 'G')])
		
	TProp = TMatches / len(z)
	AProp = AMatches / len(z)
	CProp = CMatches / len(z)
	GProp = GMatches / len(z)
			
	np.savetxt(OutFile,(AllProp,TProp,AProp,CProp,GProp))
          
if __name__ == "__main__":
        sys.exit(main())