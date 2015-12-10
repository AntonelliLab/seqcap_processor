#!/usr/local/opt/python/bin/python
#author: Tobias Hofmann, tobiashofmann@gmx.net

import os
import sys
import re
import glob
import shutil
import argparse
import ConfigParser
import commands
import subprocess
from Bio import SeqIO

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
		description="In case of tetraploid samples, run this script on the already phased data. It will run the phasing command across all samples again, creating 4 haplotype sequences for each sample.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--input',
		required=True,
		action=CompletePath,
		default=None,
		help='The folder containing the output from the phase_alleles.py script.'
	)
	parser.add_argument(
		'--config',
		required=True,
		help='The same configuration file that was used for the previous script (phase_alleles.py) containing the full paths to the following programs: samtools, clc-assembly-cell, bcftools, vcfutils, emboss, picard'
	)
	parser.add_argument(
		'--reference',
		required=True,
		action=CompletePath,
		default=None,
		help='Path to the reference fasta file.'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be safed.'
	)
	parser.add_argument(
		'--no_duplicates',
		action='store_true',
		default=False,
		help='Use this flag if you want to clean the mapped reads from all duplicates with Picard.'
	)
	parser.add_argument(
		'--conservative',
		action='store_true',
		default=False,
		help='Use this flag if you want to discard all base calls with limited certainty (covered by <3 reads). This will produce the ambiguity character "N" instead of that potential base call in the final sequence.'
	)	
	parser.add_argument(
		'--cores',
		type=int,
		default=1,
		help='For parallel processing you can choose the number of cores you want CLC to run on.'
	)
	return parser.parse_args()

# Preparation for calling input variables and files
args = get_args()
conf = ConfigParser.ConfigParser()
conf.optionxform = str
conf.read(args.config)

# Collect all program paths from control file
paths = conf.items('paths')

# Find each program name in list and define variable
samtools = ""	
for i in paths:
		if "samtools" in i:
			samtools = i[1]	
bcftools = ""	
for i in paths:
		if "bcftools" in i:
			bcftools = i[1]	
vcfutils = ""	
for i in paths:
		if "vcfutils" in i:
			vcfutils = i[1]
picard = ""	
for i in paths:
		if "picard" in i:
			picard = i[1]
clc = ""	
for i in paths:
		if "clc" in i:
			clc = i[1]
emboss = ""
for i in paths:
		if "emboss" in i:
			emboss = i[1]


reference = args.reference
phased_folder = args.input
out_dir = args.output



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Functions %%%

def phase_bam(sorted_bam_file,sample_output_folder):
	# First we need to get the name of the phasing output straight
	split_bam_path = re.split("/",sorted_bam_file)
	bam_file_name = split_bam_path[-1]
	phasing_basename_pre = re.sub('.bam$', '', bam_file_name)
	phasing_basename = re.sub('_sorted', '', phasing_basename_pre)
	phasing_basename_path = "%s/%s" %(sample_output_folder,phasing_basename)
	phasing_cmd = [
		samtools,
		"phase", 
		"-A",
		"-F",
		"-Q",
		"20",
		"-b",
		phasing_basename_path,
		sorted_bam_file
	]
	try:
		print "Phasing bam file.........."
		with open(os.path.join(sample_output_folder, "%s_screen_out.txt" %phasing_basename), 'w') as phasing_screen:
			ph = subprocess.Popen(phasing_cmd, stdout=phasing_screen)
			ph.communicate()
		print "Phasing completed."
	except:
		print "Phasing unsuccessful. Script terminated."
		sys.exit()
	
	allele_0_file = "%s.0.bam" %phasing_basename_path
	allele_1_file = "%s.1.bam" %phasing_basename_path
	allele_0_sorted_base = "%s/%s_0_sorted" %(sample_output_folder,phasing_basename)
	allele_1_sorted_base = "%s/%s_1_sorted" %(sample_output_folder,phasing_basename)


	# Sorting phased bam files:
	sort_phased_0 = "%s sort %s %s" %(samtools,allele_0_file,allele_0_sorted_base)
	sort_phased_1 = "%s sort %s %s" %(samtools,allele_1_file,allele_1_sorted_base)
	os.system(sort_phased_0)
	os.system(sort_phased_1)
	
	allele_0_sorted_file = "%s.bam" %allele_0_sorted_base
	allele_1_sorted_file = "%s.bam" %allele_1_sorted_base
	
	# Index the sorted bam files:
	index_allele_0 = "%s index %s" %(samtools,allele_0_sorted_file)
	index_allele_1 = "%s index %s" %(samtools,allele_1_sorted_file)
	os.system(index_allele_0)
	os.system(index_allele_1)	
	
	
	
	
#	
#	# Creating consensus sequences from bam-files
#	print "Creating consensus sequences from bam-files.........."
#	make_cons_0 = "%s mpileup -u -f %s %s | %s view -cg - | %s vcf2fq > %s_0.fq" %(samtools,reference,allele_0_sorted_file,bcftools,vcfutils,phasing_basename)
#	make_cons_1 = "%s mpileup -u -f %s %s | %s view -cg - | %s vcf2fq > %s_1.fq" %(samtools,reference,allele_1_sorted_file,bcftools,vcfutils,phasing_basename)
#	make_cons_unphased = "%s mpileup -u -f %s %s | %s view -cg - | %s vcf2fq > %s.fq" %(samtools,reference,sorted_bam_file,bcftools,vcfutils,bam_basename)
#	os.system(make_cons_0)
#	os.system(make_cons_1)
#	os.system(make_cons_unphased)
#	
#	# Converting fq into fasta files
#	make_fasta_cmd_0 = "seqtk seq -a %s_0.fq > %s_0.fasta" %(phasing_basename,phasing_basename)
#	make_fasta_cmd_1 = "seqtk seq -a %s_1.fq > %s_1.fasta" %(phasing_basename,phasing_basename)
#	make_fasta_cmd_unphased = "seqtk seq -a %s.fq > %s.fasta" %(bam_basename,bam_basename)
#	os.system(make_fasta_cmd_0)
#	os.system(make_fasta_cmd_1)
#	os.system(make_fasta_cmd_unphased)
#	
#	# Cleaning up output directory 
#	output_files = [val for sublist in [[os.path.join(i[0], j) for j in i[2]] for i in os.walk(sample_output_folder)] for val in sublist]
#	fasta_dir = "%s/final_fasta_files" %sample_output_folder
#	if not os.path.exists(fasta_dir):
#		os.makedirs(fasta_dir)
#	# check the names to make sure we're not deleting something improperly
#	#try:
#	#	assert "%s.fasta" %bam_basename in output_files
#	#except:
#	#	raise IOError("Output-files were not created properly.")
#	for file in output_files:
#		if file in ("%s.fasta" %bam_basename, "%s_0.fasta" %phasing_basename, "%s_1.fasta" %phasing_basename):
#			shutil.move(file,fasta_dir)
#
#	if args.conservative:
#		# Clean up the final fasta alignments and replace all uncertain base-calls (non-capitalized letters) with "N"
#		replace_uncertain_base_calls = "for fasta in $(ls %s/*.fasta); do sed -i -e '/>/! s=[actgn]=N=g' $fasta; done" %fasta_dir		
#		os.system(replace_uncertain_base_calls)
#	
#	# Standardize all the different ambiguity code-bases with N
#	standardize_all_ambiguities = "for fasta in $(ls %s/*.fasta); do sed -i -e '/>/! s=[ywrksmYWRKSM]=N=g' $fasta; done" %fasta_dir
#	os.system(standardize_all_ambiguities)
#	# Remove the unnecessary .fq files
#	remove_fq_file = "rm %s/*.fq" %sample_output_folder
#	remove_fq_file_2 = "rm %s/*/*.fq" %sample_output_folder
#	if args.no_duplicates:
#		os.system(remove_fq_file_2)
#	else: 
#		os.system(remove_fq_file)
#	return fasta_dir
#
#
#
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Workflow %%%

# Loop over every sample folder
for sample_folder in os.listdir(phased_folder):
	if "remapped" in sample_folder:
		phased_sample_folder = "%s/phased" %sample_folder
		path_phased_sample_folder = os.path.join(phased_folder,phased_sample_folder)
		sample_id = re.sub("_remapped","",sample_folder)
		sample_output_folder = "%s/%s_double_phased" %(out_dir,sample_id)
		if not os.path.exists(sample_output_folder):
			os.makedirs(sample_output_folder)
		#print sample_output_folder
		print "\n", "#" * 50
		for files in os.listdir(path_phased_sample_folder):
			if "allele_0.bam" in files:
				print "\n", "Rephasing allele_0 for sample", sample_id, "\n"
				path_allele_0 = os.path.join(path_phased_sample_folder,files)
				#print path_allele_0
				allele_0_fastas = phase_bam(path_allele_0,sample_output_folder)
				
				
			if "allele_1.bam" in files:
				print "\n", "Rephasing allele_1 for sample", sample_id, "\n"					
				path_allele_1 = os.path.join(path_phased_sample_folder,files)
				#print path_allele_1
				allele_1_fastas = phase_bam(path_allele_1,sample_output_folder)
				
				
				
				
#					path_bams = os.path.join(phased_sample_folder,bams)
#					allele_0_fastas = phase_bam(path_bams,sample_output_folder)
#					
#					# The following is for the case that no phased bam files were created, i.e. the individual is homozygous for all loci (happens when only looking at one locus or a very few)
#					allele0 = ""
#					allele1 = ""
#					# testing if phasing files were created
#					for file in os.listdir(allele_0_fastas):
#						if file.endswith(".fasta"):
#							if "allele_0" in file:
#								allele0 = file
#							if "allele_1" in file:
#								allele1 = file		
#					if allele0 == 0:
#						manage_homzygous_samples(allele_0_fastas,sample_id)
#						os.remove(os.path.join(allele_0_fastas,allele0))
#						os.remove(os.path.join(allele_0_fastas,allele1))
#					# Give fasta headers the correct format for phasing script
#					edit_fasta_headers(allele_0_fastas,sample_id)
#			
#				if "allele_1.bam" in bams:
#					print "\n", "Processing allele_1 for sample", sample_id, "\n"
#					path_bams = os.path.join(phased_sample_folder,bams)
#					allele_1_fastas = phase_bam(path_bams,sample_output_folder)
#					
#					# The following is for the case that no phased bam files were created, i.e. the individual is homozygous for all loci (happens when only looking at one locus or a very few)
#					allele0 = ""
#					allele1 = ""
#					# testing if phasing files were created
#					for file in os.listdir(allele_1_fastas):
#						if file.endswith(".fasta"):
#							if "allele_0" in file:
#								allele0 = file
#							if "allele_1" in file:
#								allele1 = file		
#					if allele0 == 0:
#						manage_homzygous_samples(allele_1_fastas,sample_id)
#						os.remove(os.path.join(allele_1_fastas,allele0))
#						os.remove(os.path.join(allele_1_fastas,allele1))
#					# Give fasta headers the correct format for phasing script
#					edit_fasta_headers(allele_1_fastas,sample_id)		
#			
#		print "\n", "#" * 50
#join_allele_fastas()