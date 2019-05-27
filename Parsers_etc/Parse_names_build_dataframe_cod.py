#!/usr/bin/python
#this file pastes values of properties of genes in the first argmunet that are found 
#in different files (the rest of the arguments) into the same file
#it assumes that each gene is represented as one line and its name is
#found somewhere followed by a separator that is ['"', '\t', '_',' '] 

import string, os, sys

def get_names(infile) :
	outlist = []
	with open(infile) as f :
		for line in f :
			outlist.append([line.split()[0]])
	return outlist		
				

def build(*args) :
	names_file = args[0][0]
	names_list = get_names(names_file)
	handles_list = []
	outfile = open(args[0][-1], 'w')
	for i in args[0][1:-1] :
		temp_handle = open(i)
		handles_list.append(temp_handle)
	all_files_lines = [f.readlines() for f in handles_list]
	for names_ind in range(len(names_list)) :
		print "Parsing : "+str(names_list[names_ind])
		ident = names_list[names_ind][0]
		for file_lines in all_files_lines :
			found_ind = False
			for line in file_lines :
				if line.find(ident+'\t') > -1 or\
				line.find(ident+'_') > -1 or\
				line.find(ident+' ') > -1 or\
				line.find(ident+'"') > -1:
					names_list[names_ind]+=line.split()		
					found_ind = True
					break
			if found_ind==False :
				names_list[names_ind]+=["NA" for i in line.split()]		
		outfile.write('\t'.join(names_list[names_ind])+'\n')
							
							
if __name__ == "__main__" :
	build(sys.argv[1:])
