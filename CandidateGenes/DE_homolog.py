#!/usr/bin/py

#Cross reference file with list of reciprocal best homologs between squash and melon
#With file on melon differential expression
#Print list of squash genes, whether they have RBH in melon, and whether the RBH is differentially expressed

import re, sys

in_file, out_file = sys.argv[1::]

fh = open('melon_expression.csv', 'r')

fh.readline()

degs = []
for line in fh:
	line = line.strip().split(",")
	degs.append(line[1])
fh.close()

melon_homos = {}
fh = open('reciprocal_TEST.txt')
for line in fh:
	line = line.strip().split("\t")
	squash_match = re.match('(^.*?)\.1$', line[0])
	squash_gene = squash_match.groups(0)[0]
	melon_homos[squash_gene] = []
	for melon_id in line[1::]:
		melon_match = re.match('(^.*?)\.2\.1$', melon_id)
		melon_gene = melon_match.groups(0)[0]
		melon_homos[squash_gene].append(melon_gene)
fh.close()	


sorted_keys = sorted(melon_homos.keys())

fh = open(in_file, 'r')
out = open(out_file, 'w')

for line in fh:
	line = line.strip().split()
	gene = line[0]
	if gene in melon_homos.keys():
		deg_match = 'FALSE'
		for homo in melon_homos[gene]:
			if homo in degs:
				deg_match = 'TRUE'
	else:
		deg_match = 'NA'
	out.write('\t'.join(line) + '\t' + deg_match + '\n')
fh.close()
out.close()


