#!/usr/bin/python2.7
import os
import argparse


from .utils import CompletePath


# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="Convert a vcf-file generated with Varscan into fasta format",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--vcf',
		required=True,
		action=CompletePath,
		default=None,
		help='The vcf file that you want to transform into fasta format'
	)
	parser.add_argument(
		'--fasta_out',
		required=True,
		action=CompletePath,
		default=None,
		help='The name of the output fasta file'
	)
	return parser.parse_args()

# Preparation for calling input variables and files
args = get_args()


vcf = args.vcf
out_fasta = args.fasta_out

with open(vcf) as f:
    content = [x.strip('\n') for x in f.readlines()]

loci = []
seq_dict = {}
for line in content:
    if '#' not in line:
        element = line.split('\t')
            # Create a list of all loci names
        seq_name = element[0]
        if seq_name not in loci:
            loci.append(seq_name)
        basecall = element[3]
        if element[4] != '.':
            basecall = element[4]
        seq_dict.setdefault(element[0],[])
        seq_dict[element[0]].append(basecall)

# Join all basecalls for each key (=locus) into one sequence and deposit in new dictionary
concat_basecalls = {}
for key, value in seq_dict.items():
	concat_basecalls[key] = "".join(value)
#print concat_basecalls

with open(out_fasta, "wb") as f:
	for k, v in concat_basecalls.items():
		f.write(">" + k+ "\n")
		f.write(v+ "\n")
