
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
from Bio import SeqIO
from .utils import CompletePath


# Get arguments
def add_arguments(parser):
	parser.add_argument(
		'--reads',
		required=True,
		action=CompletePath,
		default=None,
		help='Call the folder that contains the trimmed reads, organized in a separate subfolder for each sample. The name of the subfolder has to start with the sample name, delimited with an underscore [_].'
	)
	parser.add_argument(
		'--reference_type',
		choices=["alignment-consensus", "sample-specific", "user-ref-lib"],
		default="user-ref-lib",
		help='Please choose which type of reference you want to map the samples to. "alignment-consensus" will create a consensus sequence for each alignment file which will be used as a reference for all samples. This is recommendable when all samples are rather closely related to each other. "sample-specific" will extract the sample specific sequences from an alignment and use these as a separate reference for each individual sample. "user-ref-lib" enables to input one single fasta file created by the user which will be used as a reference library for all samples.'
	)
	parser.add_argument(
		'--reference',
		required=True,
		action=CompletePath,
		default=None,
		help='When choosing "alignment-consensus" or "sample-specific" as reference_type, this flag calls the folder containing the alignment files for your target loci (fasta-format). In case of "user-ref-lib" as reference_type, this flag calls one single fasta file that contains a user-prepared reference library which will be applied to all samples.'
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

def bam_consensus(reference,bam_file,name_base,out_dir,min_cov,args):
	# Creating consensus sequences from bam-files
	mpileup = [
		"samtools",
		"mpileup",
		"-u",
		"-f",
		reference,
		bam_file
	]
	mpileup_file = os.path.join(out_dir, "%s.mpileup" %name_base)
	with open(mpileup_file, 'w') as mpileupfile:
		mp = subprocess.Popen(mpileup, stdout=mpileupfile, stderr=subprocess.PIPE)
		mp.communicate()
		mp.wait()
	
	vcf_file = os.path.join(out_dir, "%s.vcf" %name_base)
	#bcf_cmd = [
	#	"bcftools",
	#	"view",
	#	"-c",
	#	"-g",
	#	mpileup_file
	#]
	bcf_cmd = [
		"bcftools",
		"call",
		"-c",
		mpileup_file
	]
	with open(vcf_file, 'w') as vcffile:
		vcf = subprocess.Popen(bcf_cmd, stdout=vcffile)
		vcf.communicate()
		vcf.wait()

	fq_file = os.path.join(out_dir, "%s.fq" %name_base)
	vcfutils_cmd = [
		"vcfutils.pl",
		"vcf2fq",
		"-d",
		str(min_cov),
		vcf_file
	]
	with open(fq_file, 'w') as fqfile:
		fq = subprocess.Popen(vcfutils_cmd, stdout=fqfile)
		fq.communicate()
		fq.wait()

	# Converting fq into fasta files
	fasta_file = os.path.join(out_dir,"%s_temp.fasta" %name_base)
	make_fasta = [
		"seqtk",
		"seq",
		"-a",
		fq_file,
	]
	with open(fasta_file, 'w') as fastafile:
		fasta = subprocess.Popen(make_fasta, stdout=fastafile)
		fasta.communicate()
		fasta.wait()
	
	# Create a new fasta file for final fasta printing
	final_fasta_file = os.path.join(out_dir,"%s.fasta" %name_base)

	if "allele_0" in name_base:
		fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
		sample_id = name_base.split("_")[0]
		with open(final_fasta_file, "wb") as out_file:
			for fasta in fasta_sequences:
				name, sequence = fasta.id, fasta.seq.tostring()
				name = re.sub('_consensus_sequence','',name)
				name = re.sub('_\(modified\)','',name)
				name = re.sub(r'(.*)',r'\1_%s_0 |\1' %sample_id ,name)
				sequence = re.sub('[a,c,t,g,y,w,r,k,s,m,n,Y,W,R,K,S,M]','N',sequence)
				out_file.write(">%s\n%s\n" %(name,sequence))
	
	elif "allele_1" in name_base:
		fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
		sample_id = name_base.split("_")[0]
		with open(final_fasta_file, "wb") as out_file:
			for fasta in fasta_sequences:
				name, sequence = fasta.id, fasta.seq.tostring()
				name = re.sub('_consensus_sequence','',name)
				name = re.sub('_\(modified\)','',name)
				name = re.sub(r'(.*)',r'\1_%s_1 |\1' %sample_id ,name)
				sequence = re.sub('[a,c,t,g,y,w,r,k,s,m,n,Y,W,R,K,S,M]','N',sequence)
				out_file.write(">%s\n%s\n" %(name,sequence))

	else:
		fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
		sample_id = name_base.split("_")[0]
		with open(final_fasta_file, "wb") as out_file:
			for fasta in fasta_sequences:
				name, sequence = fasta.id, fasta.seq.tostring()
				name = re.sub('_consensus_sequence','',name)
				name = re.sub('_\(modified\)','',name)
				name = re.sub(r'(.*)',r'\1_%s |\1' %sample_id ,name)
				sequence = re.sub('[a,c,t,g,y,w,r,k,s,m,n,Y,W,R,K,S,M]','N',sequence)
				out_file.write(">%s\n%s\n" %(name,sequence))
		if not args.keep_duplicates:
			dupl_folder = "%s/with_duplicates" %out_dir
			if not os.path.exists(dupl_folder):
				os.makedirs(dupl_folder)
			standard_bam = glob.glob('%s/*.bam' %out_dir)[0]
			standard_bai = glob.glob('%s/*.bam.bai' %out_dir)[0]
			shutil.move(standard_bam,dupl_folder)
			shutil.move(standard_bai,dupl_folder)
			picard_bam = glob.glob('%s/picard/*.bam' %out_dir)[0]
			picard_bai = glob.glob('%s/picard/*.bam.bai' %out_dir)[0]
			shutil.move(picard_bam,out_dir)
			shutil.move(picard_bai,out_dir)

	os.remove(fasta_file)
	os.remove(fq_file)

	return final_fasta_file


def phase_bam(sorted_bam_file,sample_output_folder,reference,min_cov,args):
	# Phasing:
	bam_basename = re.sub('.bam$', '', sorted_bam_file)
	split_sample_path = re.split("/",sorted_bam_file)
	split_file_name = split_sample_path[-1]
	phasing_file_base_pre = re.sub('.sorted.bam$', '', split_file_name)
	phasing_out_dir = "%s/phased" %(sample_output_folder)
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
	sort_phased_0 = "samtools sort -o %s %s" %(allele_0_sorted_file, allele_0_file)
	#sort_phased_1 = "samtools sort %s %s" %(allele_1_file,allele_1_sorted_base)
	sort_phased_1 = "samtools sort -o %s %s" %(allele_1_sorted_file,allele_1_file)
	os.system(sort_phased_0)
	os.system(sort_phased_1)

	# Creating index file for phased bam-files:
	index_allele0 = "samtools index %s" %(allele_0_sorted_file)
	index_allele1 = "samtools index %s" %(allele_1_sorted_file)
	os.system(index_allele0)
	os.system(index_allele1)
	
	print ("Creating consensus sequences from bam-files..........")
	name_stem = "%s_unphased" %phasing_file_base_pre

	allele0_stem = re.split("/", allele_0_sorted_base)[-1]
	allele0_stem = re.sub('_sorted', '', allele0_stem)

	allele1_stem = re.split("/", allele_1_sorted_base)[-1]
	allele1_stem = re.sub('_sorted', '', allele1_stem)
	fasta_unphased = bam_consensus(reference,sorted_bam_file,name_stem,sample_output_folder,min_cov,args)
	fasta_allele0 = bam_consensus(reference,allele_0_sorted_file,allele0_stem,sample_output_folder,min_cov,args)
	fasta_allele1 = bam_consensus(reference,allele_1_sorted_file,allele1_stem,sample_output_folder,min_cov,args)

	# Cleaning up output directory	
	output_files = [val for sublist in [[os.path.join(i[0], j) for j in i[2]] for i in os.walk(sample_output_folder)] for val in sublist]
	fasta_dir = "%s/final_fasta_files" %sample_output_folder
	if not os.path.exists(fasta_dir):
		os.makedirs(fasta_dir)
	# check the names to make sure we're not deleting something improperly
	try:
		assert fasta_unphased in output_files
	except:
		raise IOError("Output-files were not created properly.")
	for file in output_files:
		if file in (fasta_unphased,fasta_allele0,fasta_allele1):
			shutil.move(file,fasta_dir)
	return fasta_dir


def join_fastas(out_dir):
	fasta_files = [os.path.join(dirpath, f) for dirpath, dirnames, files in os.walk(out_dir) for f in files if f.endswith('.fasta')]
	allele_fastas = []
	unphased_fastas = []
	for fasta in fasta_files:
		if "allele_0" in fasta:
			allele_fastas.append(fasta)
		elif "allele_1" in fasta:
			allele_fastas.append(fasta)
		elif "_unphased" in fasta:
			unphased_fastas.append(fasta)

	joined_allele_fastas = "%s/joined_allele_fastas.fasta" %out_dir
	with open(joined_allele_fastas, 'w') as outfile:
		for fname in allele_fastas:
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)

	joined_unphased_fastas = "%s/joined_unphased_fastas.fasta" %out_dir
	with open(joined_unphased_fastas, 'w') as outfile:
		for fname in unphased_fastas:
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)	
		
	#final_merging = "for folder in $(find %s -type d -name '*_remapped'); do cat $folder/final_fasta_files/*allele*; done > %s/joined_allele_sequences_all_samples.fasta" %(out_dir,out_dir)
	#os.system(final_merging)


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
#			allele_fastas = phase_bam(sorted_bam,sample_output_folder,reference,min_cov,args)
#		## THIS iS STILL NOT WORKING PROPERLY WHEN ONLY SINGLE FILE PRESENT:
#		# The following is for the case that no phased bam files were created, i.e. the individual is homozygous for all loci (happens when only looking at one locus or a very few)
#			allele0 = ""
#			allele1 = ""
#			# testing if phasing files were created
#			for file in os.listdir(allele_fastas):
#				if file.endswith(".fasta"):
#					if "allele_0" in file:
#						allele0 = file
#					if "allele_1" in file:
#						allele1 = file
#			if allele0 == 0:
#				manage_homzygous_samples(allele_fastas,sample_id)
#				os.remove(os.path.join(allele_fastas,allele0))
#				os.remove(os.path.join(allele_fastas,allele1))
#			print ("#" * 50)
#			join_fastas(out_dir)