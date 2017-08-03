'''
Create new reference library and map raw reads against the library (reference-based assembly)
'''


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
		'--keep_duplicates',
		action='store_true',
		default=False,
		help='Use this flag if you do not want to discard all duplicate reads with Picard.'
	)
	parser.add_argument(
		'--k',
		type=float,
		default=80,
		help='Reads with a length shorter than this threshold will not be used for mapping.'
	)
	parser.add_argument(
		'--min_coverage',
		type=int,
		default=4,
		help='Set the minimum read coverage. Only positions that are covered by this number of reads will be called in the consensus sequence, otherwise the program will add an ambiguity at this position.'
	)
#	parser.add_argument(
#		'--mapper',
#		choices=["bwa"],
#		default="bwa",
#		help='Choose your desired mapping software.' 
#	)
#	parser.add_argument(
#		'--l',
#		type=float,
#		default=0.7,
#		help='(only for CLC mapper): Define the fraction of the read that has to fulfil the similarity-threshold.'
#	)
#	parser.add_argument(
#		'--s',
#		type=float,
#		default=0.9,
#		help='(only for CLC mapper): Set a similarity threshold, defining how similar the read has to be to the reference in order to be a match.'
#	)
#	parser.add_argument(
#		'--cores',
#		type=int,
#		default=1,
#		help='For parallel processing you can choose the number of cores you want CLC to run on.'
#	)


def create_reference_fasta(reference_folder,alignments):
	# Create a list of fasta files from the input directory
	file_list = [fn for fn in os.listdir(alignments) if fn.endswith(".fasta")]
	reference_list = []
	for fasta_alignment in file_list:
		sequence_name = re.sub(".fasta","",fasta_alignment)
		orig_aln = os.path.join(alignments,fasta_alignment)
		sep_reference = "%s/%s" %(reference_folder,fasta_alignment)
		reference_list.append(sep_reference)
		cons_cmd = "cons -sequence %s -outseq %s -name %s -plurality 0.1 -setcase 0.1" %(orig_aln,sep_reference,sequence_name)
		os.system(cons_cmd)
	reference = os.path.join(reference_folder,"joined_fasta_library.fasta")
	join_fastas = "cat %s/*.fasta > %s" %(reference_folder,reference)
	os.system(join_fastas)
	return reference


def create_sample_reference_fasta(reference_folder,sample_id,alignments):
	print "Creating reference library for %s .........." %sample_id
#	get the sequence header with the correct fasta id and extract sequence
#	store these sequences in separate fasta file for each locus at out_dir/reference_seqs/sample_id
#	header of sequence remains the locus name
#	remove all "-" and "?" in sequence (not in header!!)
	sample_reference_folder = os.path.join(reference_folder,sample_id)
	if not os.path.exists(sample_reference_folder):
		os.makedirs(sample_reference_folder)
	file_list = [fn for fn in os.listdir(alignments) if fn.endswith(".fasta")]
	for fasta_alignment in file_list:
		locus_id = fasta_alignment.replace(".fasta", "")
		sample_reference_fasta = os.path.join(sample_reference_folder,fasta_alignment)
		fasta_sequences = SeqIO.parse(open("%s/%s" %(alignments,fasta_alignment)),'fasta')
		outfile = open(sample_reference_fasta, 'w')
		for fasta in fasta_sequences:
			if fasta.id == sample_id:
				sequence = re.sub('[-,?]','',str(fasta.seq))
				outfile.write(">%s\n%s\n" %(locus_id,sequence))
		outfile.close()
	reference = os.path.join(sample_reference_folder,"joined_fasta_library.fasta")
	join_fastas = "cat %s/*.fasta > %s" %(sample_reference_folder,reference)
	os.system(join_fastas)
	return reference


def mapping_bwa(forward,backward,reference,sample_id,sample_output_folder, min_length, log):
	#Indexing
	command1 = ["bwa","index",reference]
	bwa_out = os.path.join(log, "bwa_screen_out.txt")
	try:
		with open(bwa_out, 'w') as logfile:
			sp1 = subprocess.Popen(command1, shell=False, stderr = subprocess.STDOUT, stdout=logfile)
			sp1.wait()
	except:
		print ("Running bwa (%s) caused an error." %bwa)
		sys.exit()

	#Mapping
	command2 = ["bwa","mem","-k",str(min_length),reference,forward,backward]
	sam_name = "%s/%s.sam" %(sample_output_folder,sample_id)
	print ("Mapping..........")
	with open(sam_name, 'w') as out, open(bwa_out, 'a') as err:
		sp2 = subprocess.Popen(command2, stderr = err, stdout=out)
		sp2.wait()

	#Converting to bam-format with samtools
	print ("Converting to bam..........")
	raw_bam = os.path.join(sample_output_folder,"%s_raw.bam" %sample_id)
	command3 = ["samtools","view","-b","-o",raw_bam,"-S",sam_name]
	sp3 = subprocess.Popen(command3,stderr=subprocess.PIPE)
	sp3.wait()
	sorted_bam = "%s/%s.sorted.bam" %(sample_output_folder,sample_id)
	command4 = ["samtools","sort", "-o", sorted_bam, raw_bam]
	sp4 = subprocess.Popen(command4)
	sp4.wait()

	#Indexing bam files
	print ("Indexing bam..........")
	command5 = ["samtools","index",sorted_bam]
	sp5 = subprocess.Popen(command5)
	sp5.wait()
	
	#Remove some big and unnecessary intermediate files
	os.remove(sam_name)
	os.remove(raw_bam)

	return sorted_bam


def mapping_clc(forward,backward,reference,sample_id,sample_output_folder):
	print ("Mapping..........")
	cas = "%s/%s.cas" %(sample_output_folder,sample_id)
	command1 = "%s -o %s -d %s -q -p fb ss 100 1000 -i %s %s -l %d -s %d --cpus %d" %(clc_mapper,cas,reference,forward,backward,length,similarity,args.cores)
	os.system(command1)

	print ("Converting to bam..........")
	bam = "%s/%s.bam" %(sample_output_folder,sample_id)
	command2 = "%s -a %s -o %s -f 33 -u" %(clc_cas_to_sam,cas,bam)
	os.system(command2)

	print ("Sorting bam..........")
	sorted_bam = "%s/%s.sorted.bam" %(sample_output_folder,sample_id)
	#sorted = "%s/%s.sorted" %(sample_output_folder,sample_id)
	command3 = "samtools sort -o %s %s" %(sorted.bam,bam)
	os.system(command3)

	print ("Indexing bam..........")
	command4 = "samtools index %s" %(sorted_bam)
	os.system(command4)

	print ("Removing obsolete files..........")
	command5 = "rm %s %s" %(cas,bam)
	os.system(command5)

	return sorted_bam


def clean_with_picard(sample_output_folder,sample_id,sorted_bam,log):

	picard_out = "%s/%s_no_dupls_sorted.bam" %(sample_output_folder,sample_id)
	dupl_log = "%s/%s_dupls.log" %(log,sample_id)
	run_picard = [
		"picard",
		"MarkDuplicates",
		"I=%s" %sorted_bam,
		"O=%s" %picard_out,
		"M=%s" %dupl_log,
		"REMOVE_DUPLICATES=true",
		"VALIDATION_STRINGENCY=LENIENT"
	]
	print ("Removing duplicate reads with Picard..........")
	with open(os.path.join(log, "picard_screen_out.txt"), 'w') as log_err_file:
		pi = subprocess.Popen(run_picard, stderr=log_err_file)
		pi.communicate()
	print ("Duplicates successfully removed.")
	# Cleaning up a bit
	has_duplicates = "%s/including_duplicate_reads" %sample_output_folder
	if not os.path.exists(has_duplicates):
		os.makedirs(has_duplicates)
	mv_duplicates_1 = 'mv %s/*.bam %s' %(sample_output_folder,has_duplicates)
	mv_duplicates_2 = 'mv %s/*.bam.bai %s' %(sample_output_folder,has_duplicates)
	mv_final = 'mv %s/*_no_dupls_sorted.bam %s' %(has_duplicates,sample_output_folder)
	os.system(mv_duplicates_1)
	os.system(mv_duplicates_2)
	os.system(mv_final)

	print ("Indexing Picard-cleaned bam..........")
	index_picard_bam = "samtools index %s" %(picard_out)
	os.system(index_picard_bam)
	dupl_bam_name = sorted_bam.split('/')[-1]
	dupl_out = os.path.join(has_duplicates,dupl_bam_name)
	return picard_out, dupl_out


def bam_consensus(reference,bam_file,name_base,out_dir,min_cov):
	print ("Generating a consensus sequence from bam-file..........")
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
				name, sequence = fasta.id, str(fasta.seq)
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
				name, sequence = fasta.id, str(fasta.seq)
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
				name, sequence = fasta.id, str(fasta.seq)
				name = re.sub('_consensus_sequence','',name)
				name = re.sub('_\(modified\)','',name)
				name = re.sub(r'(.*)',r'\1_%s |\1' %sample_id ,name)
				sequence = re.sub('[a,c,t,g,y,w,r,k,s,m,n,Y,W,R,K,S,M]','N',sequence)
				out_file.write(">%s\n%s\n" %(name,sequence))
		tmp_folder = "%s/tmp" %out_dir
		if not os.path.exists(tmp_folder):
			os.makedirs(tmp_folder)
		mpileup = glob.glob('%s/*.mpileup' %out_dir)[0]
		vcf = glob.glob('%s/*.vcf' %out_dir)[0]
		shutil.move(vcf,tmp_folder)
		shutil.move(mpileup,tmp_folder)

	os.remove(fasta_file)
	os.remove(fq_file)

	return final_fasta_file


def main(args):
	mapper = 'bwa'
	# Set working directory
	out_dir = args.output
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	else:
		raise IOError("The directory {} already exists.  Please check and remove by hand.".format(out_dir))
	
	# Get other input variables
	alignments = args.reference
	reads = args.reads
	min_length = args.k
	min_cov = args.min_coverage
	reference = ''
	reference_folder = "%s/reference_seqs" %out_dir
	if not os.path.exists(reference_folder):
		os.makedirs(reference_folder)
	if args.reference_type == "user-ref-lib":
		reference = args.reference
		manage_reference = "cp %s %s" %(reference,reference_folder)
		os.system(manage_reference)
	elif args.reference_type == "alignment-consensus":
		reference = create_reference_fasta(reference_folder,alignments)
	for subfolder in os.listdir(reads):
		path = os.path.join(reads,subfolder)
		if os.path.isdir(path):
			subfolder_path = os.path.join(reads,subfolder)
			sample_folder = subfolder
			sample_id = re.sub("_clean","",sample_folder)
			if args.reference_type == "sample-specific":
				reference = create_sample_reference_fasta(reference_folder,sample_id,alignments)
			# Loop through each sample-folder and find read-files
			sample_output_folder = "%s/%s_remapped" %(out_dir,sample_id)
			#if os.path.exists(sample_output_folder):
			#	print "\nOutput folder %s already exists. This sample (%s) will be skipped. Please delete or specify different output directory to rerun this sample\n" %(sample_output_folder,sample_id)
			#	continue
			if not os.path.exists(sample_output_folder):
				os.makedirs(sample_output_folder)
			forward = ""
			backward = ""
			for fastq in os.listdir(subfolder_path):
				if fastq.endswith('.fastq') or fastq.endswith('.fq'):
					if sample_id in fastq and "READ1.fastq" in fastq:
						forward = os.path.join(subfolder_path,fastq)
					elif sample_id in fastq and "READ2.fastq" in fastq:
						backward = os.path.join(subfolder_path,fastq)
			if forward != "" and backward != "":
				print "\n", "#" * 50
				print "Processing sample", sample_id, "\n"
				sorted_bam = ""
				log = os.path.join(sample_output_folder,'log')
				if not os.path.exists(log):
					os.makedirs(log)
				if mapper == "bwa":
					sorted_bam = mapping_bwa(forward,backward,reference,sample_id,sample_output_folder,min_length,log)
				if not args.keep_duplicates:
					sorted_bam, dupl_bam = clean_with_picard(sample_output_folder,sample_id,sorted_bam,log)
				name_stem = '%s_bam_consensus' %sample_id
				bam_consensus_file = bam_consensus(reference,sorted_bam,name_stem,sample_output_folder,min_cov)
				if not args.keep_duplicates:
					dupl_output_folder = ('/').join(dupl_bam.split('/')[:-1])
					dupl_name_stem = '%s_with_duplicates_bam_consensus' %sample_id
					bam_consensus_with_duplicates = bam_consensus(reference,dupl_bam,dupl_name_stem,dupl_output_folder,min_cov)
