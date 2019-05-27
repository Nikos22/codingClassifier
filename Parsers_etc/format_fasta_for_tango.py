import string, sys
import numpy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

counter=0
glon_counter=0
outfile = open(sys.argv[2]+"_"+str(glon_counter)+"_.txt", 'w')
for seq_record in SeqIO.parse(sys.argv[1], "fasta") :
	if counter==1000 :
		glon_counter+=1
		outfile = open(sys.argv[2]+"_"+str(glon_counter)+".txt", 'w')
		counter=0
	data = seq_record.seq._data
	new = []
	for i in range(len(data)) :
		if data[i] == '*' or data[i] == 'X':
			pass
		else :
			new.append(data[i])	
	outfile.write(seq_record.id+' N N 7 298 0.1 '+''.join(new)	+'\n')	
	counter+=1
