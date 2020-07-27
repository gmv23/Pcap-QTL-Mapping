#!/usr/bin/python

import sys, os

###################### INDEX BOTH GENOMES #########################

if os.path.isdir("./index2"):
	pass
else:
	os.mkdir("./index2")

def get_base(string):
	return string[string.rfind("/")+1::]

index_call1 = "makeblastdb -in " + sys.argv[1] + " -dbtype prot -out ./index2/" + get_base(sys.argv[1])
index_call2 = "makeblastdb -in " + sys.argv[2] + " -dbtype prot -out ./index2/" + get_base(sys.argv[2])

os.system(index_call1)
os.system(index_call2)

#################### RECIPROCAL QUERIES #############################

def get_base2(string):
	return string[string.rfind("/")+1:string.rfind("."):]

homs1_file = "./index2/" + get_base2(sys.argv[1]) + "_homs.txt"
homs2_file = "./index2/" + get_base2(sys.argv[2]) + "_homs.txt"

query_call1 = "blastp -query " + sys.argv[1] + " -db ./index2/" + get_base(sys.argv[2]) + " -outfmt 6 -out " + homs1_file
query_call2 = "blastp -query " + sys.argv[2] + " -db ./index2/" + get_base(sys.argv[1]) + " -outfmt 6 -out " + homs2_file 

os.system(query_call1)
os.system(query_call2)

#################### PRODUCING INTERMEDIATE HOMOLOG FILES FROM BLAST OUTPUT #################

homs1 = open(homs1_file, "r")
homs2 = open(homs2_file, "r")

def extract_homs(homs, out_name):

	out = open(out_name, "w")

	matches = {}

	for line in homs:
		line = line.strip().split("\t")
		query, hom, e = [line[0],line[1], float(line[10])]	
		if query not in matches.keys():
			matches[query] = [e, hom]
		elif e < matches[query][0]:
			matches[query] = [e, hom]
		elif e == matches[query][0]:
			matches[query].append(hom)			

	for item in sorted(matches.items()):
		out.write(item[0] + "\t" + "%.3e" %(item[1][0]) + "\t" + "\t".join(item[1][1::]) + "\n" ) 

	out.close()

extract_homs(homs1, "best_homologs1.txt")
extract_homs(homs2, "best_homologs2.txt")


########################### FINDING RECIPROCAL BEST HOMOLOGS #####################

d1 = {}
d2 = {}

homo1 = open("best_homologs1.txt", "r")
homo2 = open("best_homologs2.txt", "r")

for line in homo1:
	line = line.strip().split("\t")
	d1[line[0]] = line[2::]

for line in homo2:
	line = line.strip().split("\t")
	d2[line[0]] = line[2::]

homo1.close()
homo2.close()

recips = {}

for a in d1.keys():
	for b in d1[a]:
		if b in d2.keys() and a in d2[b]:
			if a not in recips.keys():
				recips[a] = [b]
			else:
				recips[a].append(b)

####################### PRINTING OUTPUT ###########################

out = open("reciprocal_best.txt", "w")

for item in sorted(recips.items()):
	out.write(item[0] + "\t" + "\t".join(item[1]) + "\n")

out.close()

