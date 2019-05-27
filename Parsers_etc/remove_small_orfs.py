#!/usr/bin/python
import sys,os
from Bio import SeqIO, AlignIO, Phylo, Alphabet
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein, IUPAC
from Bio.SeqRecord import SeqRecord

if len(sys.argv) != 4 :
	exit("wrong number of arguments")

length_cutoff = int(sys.argv[2])
outfile = sys.argv[3]
infile = sys.argv[1]
outfile_handler = open(outfile, 'w')
seqs_removed=0
for seq in SeqIO.parse(open(infile),"fasta") :
	if len(seq.seq._data) > length_cutoff :
		SeqIO.write(seq, outfile_handler, 'fasta')
	else :
		seqs_removed+=1
if seqs_removed>0 :		
	print "removed "+str(seqs_removed)+" short or missing sequences"
else :
	print "no sequences removed"				
outfile_handler.close()
