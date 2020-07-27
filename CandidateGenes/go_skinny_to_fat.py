#!/usr/bin/py

import sys

in_file, out_file = sys.argv[1::]

last_gene = open(in_file, "r").readline().strip().split("\t")[0]
go_list = []

fh = open(in_file, "r")
out = open(out_file, "w")

for line in fh:
	line = line.strip().split("\t")
	gene = line[0]
	term = line[3]	
	if gene == last_gene:
		go_list.append(term)
	else:
		out.write(last_gene + "\t" + ";".join(go_list) + "\n")
		last_gene = gene
		go_list = [term]
fh.close()
out.close()
