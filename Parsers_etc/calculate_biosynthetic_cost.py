# -*- coding: utf-8 -*-

import string
import sys
import os
from os import listdir
from os.path import isfile, join
from collections import Counter
from Bio import Phylo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

cost_dict = {'A' : 11.7,'R' : 27.3,'N' : 14.7,'D' : 12.7,'C' : 24.7,'E' : 16.3,'Q' : 15.3,\
			'G' : 11.7,'H' : 38.3,'I' : 32.3,'L' : 27.3,'K' : 30.3,'M' : 34.3,'F' : 52, \
			'P' : 20.3, 'S' : 11.7, 'T' : 18.7, 'W' : 74.3, 'Y' : 50.0, 'V' : 23.3}


prot = SeqIO.parse(sys.argv[1], 'fasta')
for seq in prot :
	cost = 0
	length=0
	for i in seq.seq._data : 
		if i not in ['*', 'X'] :
			cost+=(cost_dict[i])		
			length+=1
	cost = float(cost/length)
	print seq.id,round(cost,3)
