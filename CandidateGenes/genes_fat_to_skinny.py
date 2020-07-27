#!/usr/bin/python

import sys,re

in_file, out_file = sys.argv[1::]

fh = open(in_file, 'r')
out = open(out_file, 'w')

for line in fh:
	line = line.strip().split("\t")
	genes = line[3].strip().split(",")
	for gene in genes:
		out.write(line[0] + '\t' + gene + '\n')

fh.close()
out.close()
