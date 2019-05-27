import string, sys


infile = open(sys.argv[1], 'r')
lines = infile.readlines()
prot_id = lines[0][:-1]
struct = {'C':0, 'H':0, 'E':0}		
		
for line in lines :
	if len(line.split()) == 6 :
		struct[line.split()[2]]+=1
		
print prot_id,struct['C'],struct['H'],struct['E'],(struct['C']+struct['H']+struct['E'])					
