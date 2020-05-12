#!/usr/bin/python

# Parse vcf file and output file with allele counts per site in each sample
# Along with genotype call on parents

import sys

in_file, out_file = sys.argv[1::]

#Get_geno function
def get_geno(x):
	if x == "0/0":
		geno = "0"
	elif x == "0/1" or x == "1/0":
		geno = "1"
	elif x == "1/1":
		geno = "2"
	else:
		geno = "NA"
	return geno

fh = open(in_file, "r")

parents = ['B02_D25', 'C02_15_6015']
samples = {}
chroms = []
poss = []

for line in fh:

	line = line.strip().split("\t")

	#Skip all headers
	if line[0][0:2] == "##":
		continue

	#Line with field names
	elif line[0] == "#CHROM":

		#Put sample names into dictionary
		sample_names = line[9::]
		column = 9
		for sample_name in sample_names:
			samples[column] = [sample_name]
			column += 1

	#Pull out information from each data field
	else:
		#Pull out chrom and pos information	
		chroms.append(line[0])
		poss.append(line[1])
		
		#Loop through samples
		for column in range(9,len(line)):
			fields = line[column].strip().split(":")
			alleles = fields[2].strip().split(",")
			A1, A2 = alleles[0:2]	
			#Determine if column corresponds to parent or not; if so get genotype and not just allele counts
			if samples[column][0] in parents:
				geno = get_geno(fields[0])
				extracted = (geno, A1, A2)
			else:
				extracted = (A1, A2)
			#Append tuple with extracted info to sample dictionary
			samples[column].append(extracted)

fh.close()

keys = samples.keys()
last_key = keys[len(keys)-1]

out = open(out_file, "w")

#First write header line
out.write("CHROM" + "\t" + "POS" + "\t")
for key in keys:
	sample_name = samples[key][0]
	if sample_name in parents:
		out.write("%s_geno\t%s_A1\t%s_A2" %(sample_name, sample_name, sample_name))
	else:
		out.write("%s_A1\t%s_A2" %(sample_name, sample_name))
	if key != last_key:
		out.write("\t")

out.write("\n")

#Now write the rest of the lines
for i in range(0, len(chroms)):
	out.write(chroms[i] + "\t" + poss[i] + "\t")
	for key in keys:
		out.write("\t".join(samples[key][i + 1])) #The first element of the list is the sample name
		if key != last_key:
			out.write("\t")
	out.write("\n")
out.close()
