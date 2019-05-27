import string, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

lengths=[]
for seq_record in SeqIO.parse(sys.argv[1], "fasta") :
	print seq_record.id, len(seq_record.seq._data)
	lengths.append(len(seq_record.seq._data))

print "mean", sys.argv[1], float(sum(lengths)/len(lengths))

