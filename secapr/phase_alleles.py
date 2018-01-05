
"""
Phase remapped reads form reference-based assembly into two separate alleles. Then produce consensus sequence for each allele.
"""

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
import pickle
from Bio import SeqIO
from .utils import CompletePath
from reference_assembly import bam_consensus, join_fastas

# Get arguments
def add_arguments(parser):
	parser.add_argument(
		'--input',
		required=True,
		action=CompletePath,
		default=None,
		help='Call the folder that contains the results of the reference based assembly (output of reference_assembly function, containing the bam-files).'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be safed.'
	)
	parser.add_argument(
		'--min_coverage',
		type=int,
		default=4,
		help='Set the minimum read coverage. Only positions that are covered by this number of reads will be called in the consensus sequence, otherwise the program will add an ambiguity at this position.'
	)


def phase_bam(sorted_bam_file,sample_output_folder,min_cov,reference):
	# Phasing:
	bam_basename = re.sub('.bam$', '', sorted_bam_file)
	split_sample_path = re.split("/",sorted_bam_file)
	split_file_name = split_sample_path[-1]
	phasing_file_base_pre = re.sub('.sorted.bam$', '', split_file_name)
	if 'unphased' in phasing_file_base_pre:
		phasing_file_base_pre = re.sub('_unphased','',phasing_file_base_pre)
	phasing_out_dir = "%s/phased_bam_files" %(sample_output_folder)
	if not os.path.exists(phasing_out_dir):
		os.makedirs(phasing_out_dir)
	phasing_basename = "%s/%s_allele" %(phasing_out_dir,phasing_file_base_pre)
	phasing_cmd = [
		"samtools",
		"phase",
		"-A",
		"-F",
		"-Q",
		"20",
		"-b",
		phasing_basename,
		sorted_bam_file
	]
	try:
		print ("Phasing bam file..........")
		with open(os.path.join(phasing_out_dir, "phasing_screen_out.txt"), 'w') as phasing_screen:
			ph = subprocess.Popen(phasing_cmd, stdout=phasing_screen)
			ph.communicate()
		print ("Phasing completed.")
	except:
		print ("Phasing unsuccessful. Script terminated.")
		sys.exit()

	allele_0_file = "%s.0.bam" %phasing_basename
	allele_1_file = "%s.1.bam" %phasing_basename
	allele_0_sorted_base = "%s/%s_sorted_allele_0" %(phasing_out_dir,phasing_file_base_pre)
	allele_1_sorted_base = "%s/%s_sorted_allele_1" %(phasing_out_dir,phasing_file_base_pre)
	allele_0_sorted_file = "%s.bam" %allele_0_sorted_base
	allele_1_sorted_file = "%s.bam" %allele_1_sorted_base

	# Sorting phased bam files:
	#sort_phased_0 = "samtools sort %s %s" %(allele_0_file,allele_0_sorted_base)
	#sort_phased_1 = "samtools sort %s %s" %(allele_1_file,allele_1_sorted_base)
	sort_phased_0 = "samtools sort -o %s %s" %(allele_0_sorted_file, allele_0_file)
	sort_phased_1 = "samtools sort -o %s %s" %(allele_1_sorted_file,allele_1_file)
	os.system(sort_phased_0)
	os.system(sort_phased_1)

	# Creating index file for phased bam-files:
	index_allele0 = "samtools index %s" %(allele_0_sorted_file)
	index_allele1 = "samtools index %s" %(allele_1_sorted_file)
	os.system(index_allele0)
	os.system(index_allele1)
	
	print ("Creating consensus sequences from bam-files..........")

	allele0_stem = re.split("/", allele_0_sorted_base)[-1]
	allele0_stem = re.sub('_sorted', '', allele0_stem)

	allele1_stem = re.split("/", allele_1_sorted_base)[-1]
	allele1_stem = re.sub('_sorted', '', allele1_stem)
	fasta_allele0 = bam_consensus(reference,allele_0_sorted_file,allele0_stem,sample_output_folder,min_cov)
	fasta_allele1 = bam_consensus(reference,allele_1_sorted_file,allele1_stem,sample_output_folder,min_cov)
	
	# Cleaning up output directory	
	output_files = [val for sublist in [[os.path.join(i[0], j) for j in i[2]] for i in os.walk(sample_output_folder)] for val in sublist]
	intermediate_files = "%s/intermediate_files" %sample_output_folder
	if not os.path.exists(intermediate_files):
		os.makedirs(intermediate_files)
	# check the names to make sure we're not deleting something improperly
	try:
		assert fasta_allele0 in output_files
	except:
		raise IOError("Output-files were not created properly.")
	for file in output_files:
		if file.endswith('.mpileup') or file.endswith('.vcf'):
			shutil.move(file,intermediate_files)
	return fasta_allele0,fasta_allele1


def manage_homzygous_samples(fasta_dir, sample_id):
	fasta_sequences = SeqIO.parse(open("%s/%s.sorted.fasta" %(fasta_dir,sample_id)),'fasta')
	with open('%s/%s_joined_homozygous_alleles.fasta'%(fasta_dir,sample_id), 'w') as outfile:
		for fasta in fasta_sequences:
			name = re.split(" ", fasta.description)
			name[0] += "_0"
			fasta.description = " ".join(name)
			fasta.id += "_0"
			SeqIO.write(fasta, outfile, "fasta")
			name = re.split(" ", fasta.description)
			allele_1_name = re.sub("_0$", "_1", name[0])
			name[0] = allele_1_name
			fasta.description = " ".join(name)
			allele_1_id = re.sub("_0$", "_1", str(fasta.id))
			fasta.id = allele_1_id
			SeqIO.write(fasta, outfile, "fasta")
	outfile.close


def main(args):
	min_cov = args.min_coverage
	# Set working directory
	out_dir = args.output
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	else:
		raise IOError("The directory {} already exists.  Please check and remove by hand.".format(out_dir))
	input_folder = args.input
	sample_out_list = []
	# Iterate through all sample specific subfolders
	for subfolder in os.listdir(input_folder):
		path = os.path.join(input_folder,subfolder)
		if os.path.isdir(path):
			subfolder_path = os.path.join(input_folder,subfolder)
			if subfolder_path.endswith('_remapped') or subfolder_path.endswith('_locus_selection'):
				sample = subfolder.split('_')[0]
				sample_output_folder = os.path.join(out_dir,'%s_phased' %sample)
				if not os.path.exists(sample_output_folder):
					os.makedirs(sample_output_folder)
				sample_out_list.append(sample_output_folder)
				tmp_folder = os.path.join(subfolder_path,'tmp')
				reference_pickle = os.path.join(tmp_folder,'%s_reference.pickle' %sample)
				with open(reference_pickle, 'rb') as handle:
					reference = pickle.load(handle)
				for file in os.listdir(subfolder_path):
					if file.endswith("sorted.bam"):
						sorted_bam = file
						sorted_bam_path = os.path.join(subfolder_path,sorted_bam)
						print ("#" * 50)
						print ('Processing sample %s' %sample)
						allele_fastas = phase_bam(sorted_bam_path,sample_output_folder,min_cov,reference)

						# The following is for the case that no phased bam files were created, i.e. the individual is homozygous for all loci (happens when only looking at one locus or a very few)
						allele0 = ""
						allele1 = ""
						# testing if phasing files were created
						for file in allele_fastas:
							if file.endswith(".fasta"):
								if "allele_0" in file:
									allele0 = file
								if "allele_1" in file:
									allele1 = file
						if allele0 == 0:
							manage_homzygous_samples(allele_fastas,sample_id)
							os.remove(os.path.join(allele_fastas,allele0))
							os.remove(os.path.join(allele_fastas,allele1))
						
	join_fastas(out_dir,sample_out_list)