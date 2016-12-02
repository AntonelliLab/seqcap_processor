#!/usr/local/opt/python/bin/python
#author: Tobias Hofmann, tobiashofmann@gmx.net

import os
import sys
import re
import glob
import shutil
import argparse
import csv
import random


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Input %%%

# Complete path function
class CompletePath(argparse.Action):
    """give the full path of an input file/folder"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="This script will create consensus sequences from pairs of allele sequences, thereby turning allele alignments into consensus alignments.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--input',
		required=True,
		action=CompletePath,
		default=None,
		help='The directory containing fasta alignments'
	)
	parser.add_argument(
		'--config',
		required=True,
		help='A configuration file containing the full paths to the following programs: samtools, bcftools, vcfutils, emboss, picard. Also the paths to either clc-assembly-cell or bwa, depending on which of these two mapping softwares is chosen (see --mapper)'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be safed'
	)

	return parser.parse_args()

	
# Get arguments
args = get_args()
# Set working directory
work_dir = args.input
out_dir = args.output	
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
config = args.config

	
	
	
	
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Functions %%%
	
	
def read_fasta(fasta):
	name, seq = None, []
	for line in fasta:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

		
		
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Workflow %%%	

# Read in config file
with open(config, 'r') as c:
	conf_dict = {}
	reader = csv.reader(c, delimiter='\t')
	reader = list(reader)
	for row in reader:
		conf_dict.setdefault(row[1],[])
		conf_dict[row[1]].append(row[0])


# Create a list of all fasta files
fasta_files = []	
for fasta in os.listdir(work_dir):
	if fasta.endswith(".fasta") or fasta.endswith(".fa"):
		fasta_files.append(fasta)


for fasta in fasta_files:

	# Create an output consensus fasta file for each allele alignment
	fasta_cons_name = re.sub("allele","consensus",fasta)
	out_fasta = open(os.path.join(out_dir, fasta_cons_name), 'w')
	
	# Create a dictionary for each fasta file, where both allele sequences are assigned to the same key for each sample
	seq_dict = {}
	for seq_pair in conf_dict.values():
		key = conf_dict.keys()[conf_dict.values().index(seq_pair)]
		#print seq_pair[0], "=", key
		#print seq_pair[1], "=", key
		with open("%s/%s" %(work_dir,fasta)) as f:
			for name, seq in read_fasta(f):
				name = re.sub('>', '', name)
				if name in seq_pair:
					name = key
					seq_dict.setdefault(name,[])
					seq_dict[name].append(seq)
	# Create a consensus dict for each fasta file with the correct new header name as key and the consensus sequence of the two alleles as value
	consensus_dict = {}
	for header in seq_dict:
		consensus_dict.setdefault(header,[])
		sequence = seq_dict[header]
		allele0 = sequence[0]
		allele1 = sequence[1]
		# Find those positions where the two alleles differ from each other and make a random pick of one of the versions, simulationg a consensus sequence
		for id, base in enumerate(allele0):
			if base != allele1[id]:
				variation = [base,allele1[id]]
				base = random.choice(variation)
			consensus_dict[header].append(base)
	# Write the consensus dictionary into a fasta output file
	for cons_header in consensus_dict:
		cons_sequence =  "".join(consensus_dict[cons_header])
		cons_header = ">%s" %cons_header
		out_fasta.write(cons_header+"\n")
		out_fasta.write(cons_sequence+"\n")
	
	out_fasta.close()
		
		

	
