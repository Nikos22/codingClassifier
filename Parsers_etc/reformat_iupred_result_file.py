import string, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


for seq_record in SeqIO.parse(sys.argv[1], "fasta") :
	disordered_residues_counter = 0
	for res in seq_record.seq._data :
		if ord(res) > 96 :
			disordered_residues_counter+=1
	print seq_record.id,disordered_residues_counter

