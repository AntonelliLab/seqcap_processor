# encoding: utf-8
'''
Concatenate mutliple alignments (MSAs) into one supermatrix
'''

import os
import sys
import glob
import shutil
import argparse
import configparser
import subprocess
import subprocess
import re
#from cogent import LoadSeqs, DNA
from .utils import CompletePath


# Get arguments
def add_arguments(parser):
	parser.add_argument(
		'--input',
		required=True,
		action=CompletePath,
		default=None,
		help='The directory containing the fasta-alignment-files'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be safed'
	)


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


def main(args):
	# Set working directory
	work_dir = args.input
	out_dir = args.output
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	# Create a dictionary with the name-pattern as key and all file-names sharing that name-pattern
	fasta_dict = {}
	for fasta in os.listdir(work_dir):
		if fasta.endswith(".fasta") or fasta.endswith(".fa"):
			name_pattern = 'all_fastas'
			fasta_dict.setdefault(name_pattern,[]).append(fasta)
		#else:
		#	print "didn't work for", fasta
	print(('Found %i alignments. Concatenating all...'%len(fasta_dict['all_fastas'])))
	# Get the list of taxa names (headers) for each locus, key is out-file, values are in-files
	for key, value in fasta_dict.items():
		list_headers=[]
		# Each k is a separate fasta input file belonging to the same locus ()to be joined)
		for k in sorted(value):
			with open("%s/%s" %(work_dir,k)) as f:
				for name, seq in read_fasta(f):
					if not name in list_headers:
						list_headers.append(name)


		# Find the missing taxa in each fasta input file and simulate a sequence of correct length (only "n")
		in_fasta = os.path.join(work_dir, fasta)
		# Each k is a separate fasta input file belonging to the same locus ()to be joined)
		all_seq_dict = {}
		for k in sorted(value):
			taxa_names_single = []
			present_seq = []
			length_alignment = ""
			with open("%s/%s" %(work_dir,k)) as f:
				for name, seq in read_fasta(f):
					taxa_names_single.append(name)
					present_seq.append((name,seq))
					length_alignment = len(seq)
			# Make a list of all missing taxa in each fasta input file
			missing_taxa = []
			for header in list_headers:
				if header not in taxa_names_single:
					missing_taxa.append(header)
			simulated_seq = []
			for mistax in missing_taxa:
				fake_string = "n" * length_alignment
				simulated_seq.append((mistax,fake_string))
			all_seq = sorted(simulated_seq+present_seq)

			for seq_header, sequence in all_seq:
				all_seq_dict.setdefault(seq_header,[]).append(sequence)

		out_fasta = open(os.path.join(out_dir, "%s.fasta" %key), 'w')
		for seqname, sequences in all_seq_dict.items():
			final_sequence = "".join(sequences)
			final_sequence = final_sequence.replace('\n','')
			out_fasta.write(seqname+"\n")
			out_fasta.write(final_sequence+"\n")

		out_fasta.close()
		print(('Concatenation finished. Supermatrix printed to %s'%os.path.join(out_dir, "%s.fasta" %key)))
