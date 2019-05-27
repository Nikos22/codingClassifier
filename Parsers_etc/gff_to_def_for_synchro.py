import sys
import os
import roman

from Bio import SeqIO
from BCBio import GFF

def convert(*arg) :
	infile_handler = open(arg[0])
	inlines = infile_handler.readlines()
	outfile_handler_def = open(arg[1]+'.def', 'w')
	outfile_handler_ch = open(arg[1]+'.ch', 'w')
	chrom_dict = {}
	glob_counter=0
	outfile_handler_def.write("type	name		chr	start	end	strand	sens	IDg/chr	IDg/all	IDf/all\n")
	for line in inlines :	
		if line.split('\t')[2] != "" :
			glob_counter+=1
			chrom = line.split()[0].zfill(3)
			if chrom in chrom_dict.keys() :
				chrom_dict[chrom]+=1
			else :	
				chrom_dict[chrom]=1
			f_type = line.split('\t')[2]	
			f_name = line.split('\t')[8].split(';')[5][5:]
			chrom_id = str(chrom_dict[chrom]).zfill(5)
			glob_id = str(glob_counter).zfill(5)
			strand = line.split('\t')[6]
			start = line.split('\t')[3]
			end = line.split('\t')[4]
			if strand=="+" :
				strand_letter = "f"
			else :	
				strand_letter = "t"
			outfile_handler_def.write(str(f_type)+\
								      "\t"+str(f_name)+\
								      '\t'+str(chrom)+\
								      '\t'+str(start)+\
								      '\t'+str(end)+\
								      '\t'+strand+\
								      '\t'+strand_letter+\
								      '\t'+str(chrom_id)+\
								      '\t'+str(glob_id)+\
								      '\t'+str(glob_id)+'\n')			
	outfile_handler_ch.write('\t'.join([x for x in sorted(chrom_dict.keys())])+'\t\n')
	cum = []
	for x in range(len(sorted(chrom_dict))) :
		key = sorted(chrom_dict)[x]
		cum_feat = 1
		for y in range(int(key)+1) :
			if str(y).zfill(3) in chrom_dict.keys() :
				cum_feat+=int(chrom_dict[str(y).zfill(3)])
		cum.append(cum_feat)
	outfile_handler_ch.write('\t'.join([str(x) for x in cum])+'\t\n')		
	outfile_handler_ch.write('\t'.join([str(x) for x in cum])+'\t\n')		
	print chrom_dict
	#+str('\t'.join([str(chrom_dict[x]) for x in range(len(sorted(chrom_dict)]))))+'\n'+str('\t'.join([str(chrom_dict[x]) for x in sorted(chrom_dict)]))+'\n')
	outfile_handler_ch.close()
	outfile_handler_def.close()
			
	
	
if __name__ == "__main__" :
	if len(sys.argv) != 3 :
		exit("Wrong number of arguments")
	convert(sys.argv[1], sys.argv[2])	
