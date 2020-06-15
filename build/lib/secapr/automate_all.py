# encoding: utf-8
#author: Tobias Andermann, tobias.andermann@bioenv.gu.se
"""
This script automates the complete secapr pipeline, producing MSAs (allele, contig and BAM-consensus) from FASTQ files
"""

import os
import sys
import glob
import fnmatch
import shutil
import configparser
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
		'--assembler',
		choices=["abyss", "trinity"],
		default="abyss",
		help='The assembler to use for de-novo assembly (default=abyss).'
	)
	parser.add_argument(
		'--cores',
		type=int,
		default=1,
		help='Number of computational cores for parallelization of computation.'
	)


def log_command(command,logfile):
	with open(logfile, 'a') as log_file:
		initiation_cmd = ['echo "%s"' %command]
		a = subprocess.Popen(initiation_cmd, stdout=log_file, stderr=log_file, shell=True)
		a.communicate()

def main(args):
	# Get arguments
	input_folder = args.input
	output_folder = args.output
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)
	reference = args.reference
	setting = args.setting
	cores = args.cores
	assembler = args.assembler
	central_log_file = os.path.join(output_folder, "secapr_log.txt")
	command_log_file = os.path.join(output_folder, "secapr_command_list.txt")
	print(('Starting SECAPR in automated mode with setting "%s"' %setting))
	with open(central_log_file, 'w') as log_err_file:
		initiation_cmd = ['echo "Running SECAPR automated with setting %s"' %setting]
		a = subprocess.Popen(initiation_cmd, stdout=log_err_file, stderr=log_err_file, shell=True)
		a.communicate()
	with open(command_log_file, 'w') as command_file:
		initiation_cmd = ['echo "Running SECAPR automated with setting %s"' %setting]
		a = subprocess.Popen(initiation_cmd, stdout=command_file, stderr=command_file, shell=True)
		a.communicate()

	if setting == 'medium':
		# assemble_reads
		print('Running assembly ...')
		assembly_out = os.path.join(output_folder,'assembly/')
		command = ["secapr assemble_reads --input %s --output %s --assembler %s --cores %i" %(input_folder,assembly_out,assembler,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			b = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			b.communicate()
		
		# find_target_contigs
		print('Extracting target contigs ...')
		find_target_out = os.path.join(output_folder,'target_contigs/')
		command = ["secapr find_target_contigs --contigs %s --reference %s --output %s --keep-duplicates" %(assembly_out,reference,find_target_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			d = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			d.communicate()

		# align_sequences
		print('Generating contig MSAs...')
		target_contig_fasta = os.path.join(find_target_out,'extracted_target_contigs_all_samples.fasta')
		contig_align_out = os.path.join(output_folder,'alignments/contig_MSAs')
		command = ["secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i" %(target_contig_fasta,contig_align_out,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			e = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			e.communicate()	

		# add_missing_sequences
		print('Filling alignments with ambiguities for missing sequences...')
		contig_align_complete_out = os.path.join(output_folder,'alignments/contig_MSAs_filled')
		command = ["secapr add_missing_sequences --input %s --output %s" %(contig_align_out,contig_align_complete_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			f = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			f.communicate()

		# reference_assembly
		print('Running reference-assembly, mapping reads to contig alignment consensus...')
		reference_assembly_out = os.path.join(output_folder,'reference_assembly')
		command = ["secapr reference_assembly --reads %s --reference_type alignment-consensus --reference %s --output %s --min_coverage 4" %(input_folder,contig_align_out,reference_assembly_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			g = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			g.communicate()

		# align_sequences
		print('Generating MSAs from BAM-consensus sequences...')
		target_bam_cons_fasta = os.path.join(reference_assembly_out,'joined_unphased_fastas.fasta')
		bam_consensus_align_out = os.path.join(output_folder,'alignments/bam_consensus_MSAs')
		command = ["secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i" %(target_bam_cons_fasta,bam_consensus_align_out,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			h = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			h.communicate()

		# add_missing_sequences
		print('Filling alignments with ambiguities for missing sequences...')
		bam_consensus_align_complete_out = os.path.join(output_folder,'alignments/bam_consensus_MSAs_filled')
		command = ["secapr add_missing_sequences --input %s --output %s" %(bam_consensus_align_out,bam_consensus_align_complete_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			i = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			i.communicate()

		# phase_reads
		print ('Phasing reads to produce allele sequences...')
		phasing_out = os.path.join(output_folder,'phasing_results')
		command = ['secapr phase_alleles --input %s --output %s --min_coverage 3' %(reference_assembly_out,phasing_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			j = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			j.communicate()

		# align_sequences
		print('Generating MSAs from phased allele data...')
		target_allele_fasta = os.path.join(phasing_out,'joined_allele_fastas.fasta')
		allele_align_out = os.path.join(output_folder,'alignments/phased_allele_MSAs')
		command = ['secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i' %(target_allele_fasta,allele_align_out,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			k = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			k.communicate()

		# add_missing_sequences
		print('Filling alignments with ambiguities for missing sequences...')
		allele_align_complete_out = os.path.join(output_folder,'alignments/phased_allele_MSAs_filled')
		command = ["secapr add_missing_sequences --input %s --output %s" %(allele_align_out,allele_align_complete_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			l = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			l.communicate()
		




	elif setting == 'relaxed':
		# assemble_reads
		print('Running assembly ...')
		assembly_out = os.path.join(output_folder,'assembly/')
		command = ["secapr assemble_reads --input %s --output %s --assembler %s --cores %i" %(input_folder,assembly_out,assembler,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			b = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			b.communicate()
		
		# find_target_contigs
		print('Extracting target contigs ...')
		find_target_out = os.path.join(output_folder,'target_contigs/')
		command = ["secapr find_target_contigs --contigs %s --reference %s --output %s --min-coverage 70 --min-identity 70 --keep-duplicates" %(assembly_out,reference,find_target_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			d = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			d.communicate()

		# align_sequences
		print('Generating contig MSAs...')
		target_contig_fasta = os.path.join(find_target_out,'extracted_target_contigs_all_samples.fasta')
		contig_align_out = os.path.join(output_folder,'alignments/contig_MSAs')
		command = ["secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i" %(target_contig_fasta,contig_align_out,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			e = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			e.communicate()	

		# add_missing_sequences
		print('Filling alignments with ambiguities for missing sequences...')
		contig_align_complete_out = os.path.join(output_folder,'alignments/contig_MSAs_filled')
		command = ["secapr add_missing_sequences --input %s --output %s" %(contig_align_out,contig_align_complete_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			f = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			f.communicate()

		# reference_assembly
		print('Running reference-assembly, mapping reads to contig alignment consensus...')
		reference_assembly_out = os.path.join(output_folder,'reference_assembly')
		command = ["secapr reference_assembly --reads %s --reference_type alignment-consensus --reference %s --output %s --min_coverage 4 --k 30 --b 1" %(input_folder,contig_align_out,reference_assembly_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			g = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			g.communicate()

		# align_sequences
		print('Generating MSAs from BAM-consensus sequences...')
		target_bam_cons_fasta = os.path.join(reference_assembly_out,'joined_unphased_fastas.fasta')
		bam_consensus_align_out = os.path.join(output_folder,'alignments/bam_consensus_MSAs')
		command = ["secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i" %(target_bam_cons_fasta,bam_consensus_align_out,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			h = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			h.communicate()

		# add_missing_sequences
		print('Filling alignments with ambiguities for missing sequences...')
		bam_consensus_align_complete_out = os.path.join(output_folder,'alignments/bam_consensus_MSAs_filled')
		command = ["secapr add_missing_sequences --input %s --output %s" %(bam_consensus_align_out,bam_consensus_align_complete_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			i = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			i.communicate()

		# phase_reads
		print ('Phasing reads to produce allele sequences...')
		phasing_out = os.path.join(output_folder,'phasing_results')
		command = ['secapr phase_alleles --input %s --output %s --min_coverage 3' %(reference_assembly_out,phasing_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			j = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			j.communicate()

		# align_sequences
		print('Generating MSAs from phased allele data...')
		target_allele_fasta = os.path.join(phasing_out,'joined_allele_fastas.fasta')
		allele_align_out = os.path.join(output_folder,'alignments/phased_allele_MSAs')
		command = ['secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i' %(target_allele_fasta,allele_align_out,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			k = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			k.communicate()

		# add_missing_sequences
		print('Filling alignments with ambiguities for missing sequences...')
		allele_align_complete_out = os.path.join(output_folder,'alignments/phased_allele_MSAs_filled')
		command = ["secapr add_missing_sequences --input %s --output %s" %(allele_align_out,allele_align_complete_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			l = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			l.communicate()






	elif setting == 'conservative':
		# assemble_reads
		print('Running assembly ...')
		assembly_out = os.path.join(output_folder,'assembly/')
		command = ["secapr assemble_reads --input %s --output %s --assembler %s --cores %i" %(input_folder,assembly_out,assembler,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			b = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			b.communicate()
		
		# find_target_contigs
		print('Extracting target contigs ...')
		find_target_out = os.path.join(output_folder,'target_contigs/')
		command = ["secapr find_target_contigs --contigs %s --reference %s --output %s --min-coverage 90 --min-identity 90 --keep-duplicates" %(assembly_out,reference,find_target_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			d = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			d.communicate()

		# align_sequences
		print('Generating contig MSAs...')
		target_contig_fasta = os.path.join(find_target_out,'extracted_target_contigs_all_samples.fasta')
		contig_align_out = os.path.join(output_folder,'alignments/contig_MSAs')
		command = ["secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i" %(target_contig_fasta,contig_align_out,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			e = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			e.communicate()	

		# add_missing_sequences
		print('Filling alignments with ambiguities for missing sequences...')
		contig_align_complete_out = os.path.join(output_folder,'alignments/contig_MSAs_filled')
		command = ["secapr add_missing_sequences --input %s --output %s" %(contig_align_out,contig_align_complete_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			f = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			f.communicate()

		# reference_assembly
		print('Running reference-assembly, mapping reads to contig alignment consensus...')
		reference_assembly_out = os.path.join(output_folder,'reference_assembly')
		command = ["secapr reference_assembly --reads %s --reference_type alignment-consensus --reference %s --output %s --min_coverage 4 --k 60" %(input_folder,contig_align_out,reference_assembly_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			g = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			g.communicate()

		# align_sequences
		print('Generating MSAs from BAM-consensus sequences...')
		target_bam_cons_fasta = os.path.join(reference_assembly_out,'joined_unphased_fastas.fasta')
		bam_consensus_align_out = os.path.join(output_folder,'alignments/bam_consensus_MSAs')
		command = ["secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i" %(target_bam_cons_fasta,bam_consensus_align_out,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			h = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			h.communicate()

		# add_missing_sequences
		print('Filling alignments with ambiguities for missing sequences...')
		bam_consensus_align_complete_out = os.path.join(output_folder,'alignments/bam_consensus_MSAs_filled')
		command = ["secapr add_missing_sequences --input %s --output %s" %(bam_consensus_align_out,bam_consensus_align_complete_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			i = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			i.communicate()

		# phase_reads
		print ('Phasing reads to produce allele sequences...')
		phasing_out = os.path.join(output_folder,'phasing_results')
		command = ['secapr phase_alleles --input %s --output %s --min_coverage 3' %(reference_assembly_out,phasing_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			j = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			j.communicate()

		# align_sequences
		print('Generating MSAs from phased allele data...')
		target_allele_fasta = os.path.join(phasing_out,'joined_allele_fastas.fasta')
		allele_align_out = os.path.join(output_folder,'alignments/phased_allele_MSAs')
		command = ['secapr align_sequences --sequences %s --output %s --aligner mafft --output-format fasta --no-trim --ambiguous --cores %i' %(target_allele_fasta,allele_align_out,cores)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			k = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			k.communicate()

		# add_missing_sequences
		print('Filling alignments with ambiguities for missing sequences...')
		allele_align_complete_out = os.path.join(output_folder,'alignments/phased_allele_MSAs_filled')
		command = ["secapr add_missing_sequences --input %s --output %s" %(allele_align_out,allele_align_complete_out)]
		log_command(command[0],command_log_file)
		with open(central_log_file, 'a') as log_err_file:
			l = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file, shell=True)
			l.communicate()

# assembly
# find target contigs (needs fasta reference)
# make contig alignments
# add missing sequences
# remapping
# phasing
# make phased alignments
# add missing sequences
