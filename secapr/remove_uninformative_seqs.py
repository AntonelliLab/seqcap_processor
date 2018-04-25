#author: Tobias Andermann, tobias.andermann@bioenv.gu.se

import os
import sys
import re
import glob
import shutil
import argparse
from Bio import SeqIO

from .utils import CompletePath


# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="Set the maximum fraction of missing data that you want to allow in an alignment and drop all sequences above this threshold.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--alignment',
		required=True,
		action=CompletePath,
		default=None,
		help='The alignment in fasta format.'
	)
	parser.add_argument(
		'--maximum_missing',
		type=float,
		default=0.8,
		help='Define the maximal fraction of missing data that you want to allow. All sequences below this threshold will be exported into a new alignment.'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be safed.'
	)
	return parser.parse_args()
args = get_args()

# Set working directory
out_dir = args.output
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Get other input variables
alignment = args.alignment
max_mis = args.maximum_missing


def manage_homzygous_samples(fasta,threshold,output):
	fasta_alignment = SeqIO.parse(open(fasta),'fasta')
	with open('%s/cleaned_alignment_all_sequences_less_than_%f_missing_data.fasta' %(output,threshold), 'w') as outfile:
		final_seqs = {}
		for sample in fasta_alignment:
			header = sample.description
			sequence = sample.seq
			chars = list(sequence)
			bad_chars = []
			for char in chars:
				if char not in ['A','C','T','G','a','c','t','g']:
					bad_chars.append(char)
			sequence_length = float(len(chars))
			count_bad_chars = float(len(bad_chars))
			fraction = float(count_bad_chars/sequence_length)
			if fraction <= threshold:
				final_seqs.setdefault(header,[]).append(sequence)
			else:
				print "Dropped sequence for", header
		for seqname, seq in final_seqs.iteritems():
			sequence = str(seq[0])
			outfile.write(">"+seqname+"\n")
			outfile.write(sequence+"\n")
	outfile.close

manage_homzygous_samples(alignment,max_mis,out_dir)
