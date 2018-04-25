# encoding: utf-8
#author: Tobias Andermann, tobias.andermann@bioenv.gu.se
"""
This script automates the complete secapr pipeline,
"""

import os
import sys
import glob
import fnmatch
import shutil
import ConfigParser
from .utils import CompletePath
import subprocess
import pandas as pd
import numpy as np

def add_arguments(parser):
	parser.add_argument(
		'--input',
		required=True,
		action=CompletePath,
		default=None,
		help='The directory containing cleaned fastq files'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where all intermediate and final data files will be stored'
	)
	parser.add_argument(
		'--reference',
		required=True,
		action=CompletePath,
		default=None,
		help='Provide a reference library (FASTA) containing sequences for the genes of interest (required to find contigs matching targeted regions).'
	)
	parser.add_argument(
		'--setting',
		choices=["relaxed", "medium", "conservative"],
		default="medium",
		help='The setting you want to run SECAPR on. "relaxed" uses very non-restrictive default values (use when samples are expected to differ considerably from provided reference or are covering wide evolutionary range, e.g. different families or orders). "conservative" is very restrictive and can be used when samples are closely related and match provided reference very well.'
	)
	parser.add_argument(
		'--cores',
		type=int,
		default=1,
		help='Number of computational cores for parallelization of computation.'
	)


def main(args):
	# Get arguments
	input_folder = args.input
	output_folder = args.output
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)
	reference = args.reference
	setting = args.setting
	cores = args.cores
	central_log_file = os.path.join(output_folder, "secapr_log.txt")
	print('Starting SECAPR in automated mode with setting "%s"' %setting)
	with open(central_log_file, 'w') as log_err_file:
		initiation_cmd = ['echo "Running SECAPR automated with setting %s"' %setting]
		a = subprocess.Popen(initiation_cmd, stdout=log_err_file, stderr=log_err_file, shell=True)
		a.communicate()

	if setting == 'medium':
		# assemble_reads
		print('Running assembly ...')
		assembly_out = os.path.join(output_folder,'assembly/')
		command = ["secapr assemble_reads --input %s --output %s --assembler abyss --cores %i" %(input_folder,assembly_out,cores)]
		with open(central_log_file, 'a') as log_err_file:
			b = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			b.communicate()

		# find_target_contigs
		print('Extracting target contigs ...')
		find_target_out = os.path.join(output_folder,'target_contigs/')
		command = ["secapr find_target_contigs --contigs %s --reference %s --output %s --keep-duplicates" %(assembly_out,reference,find_target_out)]
		with open(central_log_file, 'a') as log_err_file:
			c = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			c.communicate()

		# align_sequences
		print('Generating contig MSAs')
		target_contig_fasta = os.path.join(find_target_out,'extracted_target_contigs_all_samples.fasta')
		contig_align_out = os.path.join(output_folder,'alignments/contig_MSAs')
		command = ["secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i" %(target_contig_fasta,contig_align_out,cores)]
		with open(central_log_file, 'a') as log_err_file:
			c = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			c.communicate()


	elif setting == 'relaxed':
		pass


	elif setting == 'conservative':
		pass

# assembly
# find target contigs (needs fasta reference)
# make contig alignments
# add missing sequences
# remapping
# phasing
# make phased alignments
# add missing sequences
