import sys
import os
from Bio import SeqIO

outfile_handle = open(sys.argv[3], 'w')
infile_handle = open(sys.argv[1], 'r')
def_file_handle = open(sys.argv[2], 'r')
def_lines = def_file_handle.readlines()
for seq in SeqIO.parse(infile_handle, 'fasta') :
	for line in def_lines :
		if seq.id == line.split()[1] :
			seq.id = '\t'.join(line.split()[1:])
			seq.description = ""
			SeqIO.write(seq, outfile_handle, 'fasta')
			break

outfile_handle.close()			
		
