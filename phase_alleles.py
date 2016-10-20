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
		description="Map raw reads against reference sequence and phase the reads into two separate alleles. Then produce consensus sequence for each allele.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--reads',
		required=True,
		action=CompletePath,
		default=None,
		help='Call the folder that contains the trimmed reads, organized in a separate subfolder for each sample. The name of the subfolder has to start with the sample name, delimited with an underscore [_].'
	)
	parser.add_argument(
		'--mapper',
		choices=["clc", "bwa"],
		default="bwa",
		help='Choose your desired mapping software. Both mappers perform well in our experience, but you will need to purchase a licence for clc-assembly-cell in case you want to download it. We therefore recommend the cost-free alternative bwa.'
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
		'--config',
		required=True,
		help='A configuration file containing the full paths to the following programs: samtools, bcftools, vcfutils, emboss, picard. Also the paths to either clc-assembly-cell or bwa, depending on which of these two mapping softwares is chosen (see --mapper)'
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
		'--min_coverage',
		type=int,
		default=4,
		help='Set the minimum read coverage. Only positions that are covered by this number of reads will be called in the consensus sequence, otherwise the program will add an ambiguity at this position.'
	)
	parser.add_argument(
		'--k',
		type=float,
		default=80,
		help='(only for bwa mapper): Reads with a length shorter than this threshold will not be used for mapping.'
	)
	parser.add_argument(
		'--l',
		type=float,
		default=0.7,
		help='(only for CLC mapper): Define the fraction of the read that has to fulfil the similarity-threshold.'
	)
	parser.add_argument(
		'--s',
		type=float,
		default=0.9,
		help='(only for CLC mapper): Set a similarity threshold, defining how similar the read has to be to the reference in order to be a match.'
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
emboss = ""
for i in paths:
	if "emboss" in i:
		emboss = i[1]
seqtk = ""
for i in paths:
	if "seqtk" in i:
		seqtk = i[1]
clc = ""
for i in paths:
	if "clc" in i:
		clc = i[1]
bwa = ""
for i in paths:
	if "bwa" in i:
		bwa = i[1]


# Specify the different CLC commands
mark_duplicates = os.path.join(picard,"MarkDuplicates.jar")
cons = os.path.join(emboss,"bin/cons")
clc_mapper = os.path.join(clc,"clc_mapper")
clc_cas_to_sam = os.path.join(clc,"clc_cas_to_sam")
#clc_assembler = "%s/clc_assembler" %clc

# Set working directory
out_dir = args.output
if not os.path.exists(out_dir):
	os.makedirs(out_dir)

# Get other input variables
alignments = args.reference
reads = args.reads
length = args.l
similarity = args.s
min_length = args.k
min_cov = args.min_coverage

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Functions %%%


def create_reference_fasta(reference_folder):
	# Create a list of fasta files from the input directory
	file_list = [fn for fn in os.listdir(alignments) if fn.endswith(".fasta")]
	reference_list = []
	for fasta_alignment in file_list:
		sequence_name = re.sub(".fasta","",fasta_alignment)
		orig_aln = os.path.join(alignments,fasta_alignment)
		sep_reference = "%s/%s" %(reference_folder,fasta_alignment)
		reference_list.append(sep_reference)
		cons_cmd = "%s -sequence %s -outseq %s -name %s -plurality 0.1 -setcase 0.1" %(cons,orig_aln,sep_reference,sequence_name)
		os.system(cons_cmd)
	reference = os.path.join(reference_folder,"joined_fasta_library.fasta")
	join_fastas = "cat %s/*.fasta > %s" %(reference_folder,reference)
	os.system(join_fastas)
	return reference


def create_sample_reference_fasta(reference_folder,sample_id):
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
				sequence = re.sub('[-,?]','',fasta.seq)
				outfile.write(">%s\n%s\n" %(locus_id,sequence))
		outfile.close()
	reference = os.path.join(sample_reference_folder,"joined_fasta_library.fasta")
	join_fastas = "cat %s/*.fasta > %s" %(sample_reference_folder,reference)
	os.system(join_fastas)
	return reference


def mapping_bwa(forward,backward,reference,sample_id,sample_output_folder):
	#Indexing
	command1 = [bwa,"index",reference]
	bwa_out = os.path.join(sample_output_folder, "bwa_screen_out.txt")
	try:
		with open(bwa_out, 'w') as logfile:
			sp1 = subprocess.Popen(command1, shell=False, stderr = subprocess.STDOUT, stdout=logfile)
			sp1.wait()
	except:
		print "Running bwa (%s) caused an error. Please check your bwa path specification and version in the control file." %bwa
		sys.exit()

	#Mapping
	command2 = [bwa,"mem","-k",str(min_length),reference,forward,backward]
	sam_name = "%s/%s.sam" %(sample_output_folder,sample_id)
	print "Mapping.........."
	with open(sam_name, 'w') as out, open(bwa_out, 'a') as err:
		sp2 = subprocess.Popen(command2, stderr = err, stdout=out)
		sp2.wait()

	#Converting to bam-format with samtools
	print "Converting to bam.........."
	raw_bam = os.path.join(sample_output_folder,"%s_raw.bam" %sample_id)
	command3 = [samtools,"view","-b","-o",raw_bam,"-S",sam_name]
	sp3 = subprocess.Popen(command3,stderr=subprocess.PIPE)
	sp3.wait()

	bam_core = "%s/%s.sorted" %(sample_output_folder,sample_id)
	command4 = [samtools,"sort",raw_bam,bam_core]
	sp4 = subprocess.Popen(command4)
	sp4.wait()

	#Indexing bam files
	print "Indexing bam.........."
	sorted_bam = "%s/%s.sorted.bam" %(sample_output_folder,sample_id)
	command5 = [samtools,"index",sorted_bam]
	sp5 = subprocess.Popen(command5)
	sp5.wait()
	
	#Remove some big and unnecessary intermediate files
	os.remove(sam_name)
	os.remove(raw_bam)

	return sorted_bam


def mapping_clc(forward,backward,reference,sample_id,sample_output_folder):
	print "Mapping.........."
	cas = "%s/%s.cas" %(sample_output_folder,sample_id)
	command1 = "%s -o %s -d %s -q -p fb ss 100 1000 -i %s %s -l %d -s %d --cpus %d" %(clc_mapper,cas,reference,forward,backward,length,similarity,args.cores)
	os.system(command1)

	print "Converting to bam.........."
	bam = "%s/%s.bam" %(sample_output_folder,sample_id)
	command2 = "%s -a %s -o %s -f 33 -u" %(clc_cas_to_sam,cas,bam)
	os.system(command2)

	print "Sorting bam.........."
	sorted = "%s/%s.sorted" %(sample_output_folder,sample_id)
	command3 = "%s sort %s %s" %(samtools,bam,sorted)
	os.system(command3)

	print "Indexing bam.........."
	sorted_bam = "%s/%s.sorted.bam" %(sample_output_folder,sample_id)
	command4 = "%s index %s" %(samtools,sorted_bam)
	os.system(command4)

	print "Removing obsolete files.........."
	command5 = "rm %s %s" %(cas,bam)
	os.system(command5)

	return sorted_bam


def clean_with_picard(sample_output_folder,sample_id,sorted_bam):
	picard_folder = "%s/picard" %sample_output_folder
	if not os.path.exists(picard_folder):
		os.makedirs(picard_folder)
	picard_log_folder = "%s/log" %picard_folder
	if not os.path.exists(picard_log_folder):
		os.makedirs(picard_log_folder)
	picard_out = "%s/%s_no_dupls_sorted.bam" %(picard_folder,sample_id)
	dupl_log = "%s/%s_dupls.log" %(picard_log_folder,sample_id)
	run_picard = [
		"java",
		"-jar",
		mark_duplicates,
		"I=%s" %sorted_bam,
		"O=%s" %picard_out,
		"M=%s" %dupl_log,
		"REMOVE_DUPLICATES=true",
		"VALIDATION_STRINGENCY=LENIENT"
	]
	try:
		print "Removing duplicate reads with Picard.........."
		with open(os.path.join(picard_log_folder, "picard_screen_out.txt"), 'w') as log_err_file:
			pi = subprocess.Popen(run_picard, stderr=log_err_file)
			pi.communicate()
		print "Duplicates successfully removed."
	except:
		print "Running Picard caused an error. Please check your picard path specification in the control file."
		sys.exit()

	print "Indexing Picard-cleaned bam.........."
	index_picard_bam = "%s index %s" %(samtools,picard_out)
	os.system(index_picard_bam)
	return picard_out


def bam_consensus(reference,bam_file,name_base,out_dir,min_cov):
	# Creating consensus sequences from bam-files
	mpileup = [
		samtools,
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
	bcf_cmd = [
		bcftools,
		"view",
		"-c",
		"-g",
		mpileup_file
	]
	with open(vcf_file, 'w') as vcffile:
		vcf = subprocess.Popen(bcf_cmd, stdout=vcffile)
		vcf.communicate()
		vcf.wait()

	fq_file = os.path.join(out_dir, "%s.fq" %name_base)
	vcfutils_cmd = [
		vcfutils,
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
		seqtk,
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


def phase_bam(sorted_bam_file,sample_output_folder):
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
		samtools,
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
		print "Phasing bam file.........."
		with open(os.path.join(phasing_out_dir, "phasing_screen_out.txt"), 'w') as phasing_screen:
			ph = subprocess.Popen(phasing_cmd, stdout=phasing_screen)
			ph.communicate()
		print "Phasing completed."
	except:
		print "Phasing unsuccessful. Script terminated."
		sys.exit()

	allele_0_file = "%s.0.bam" %phasing_basename
	allele_1_file = "%s.1.bam" %phasing_basename
	allele_0_sorted_base = "%s/%s_sorted_allele_0" %(phasing_out_dir,phasing_file_base_pre)
	allele_1_sorted_base = "%s/%s_sorted_allele_1" %(phasing_out_dir,phasing_file_base_pre)
	allele_0_sorted_file = "%s.bam" %allele_0_sorted_base
	allele_1_sorted_file = "%s.bam" %allele_1_sorted_base

	# Sorting phased bam files:
	sort_phased_0 = "%s sort %s %s" %(samtools,allele_0_file,allele_0_sorted_base)
	sort_phased_1 = "%s sort %s %s" %(samtools,allele_1_file,allele_1_sorted_base)
	os.system(sort_phased_0)
	os.system(sort_phased_1)

	# Creating index file for phased bam-files:
	index_allele0 = "%s index %s" %(samtools,allele_0_sorted_file)
	index_allele1 = "%s index %s" %(samtools,allele_1_sorted_file)
	os.system(index_allele0)
	os.system(index_allele1)
	
	print "Creating consensus sequences from bam-files.........."
	name_stem = "%s_unphased" %phasing_file_base_pre

	allele0_stem = re.split("/", allele_0_sorted_base)[-1]
	allele0_stem = re.sub('_sorted', '', allele0_stem)

	allele1_stem = re.split("/", allele_1_sorted_base)[-1]
	allele1_stem = re.sub('_sorted', '', allele1_stem)
	fasta_unphased = bam_consensus(reference,sorted_bam_file,name_stem,sample_output_folder,min_cov)
	fasta_allele0 = bam_consensus(reference,allele_0_sorted_file,allele0_stem,sample_output_folder,min_cov)
	fasta_allele1 = bam_consensus(reference,allele_1_sorted_file,allele1_stem,sample_output_folder,min_cov)

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


def assembly_clc(forward,backward):
	print "De-novo assembly with CLC.........."
	output_fasta = "%s/%s_contigs.fasta" %(sample_output_folder,sample_id)
	distance_file = "%s/%s-estimated-distances.txt" %(sample_output_folder,sample_id)
	command = "%s -o %s -q -i %s %s -e %s -p fb ss 150 800 --cpus %d" %(clc_assembler,output_fasta,forward,backward,distance_file,args.cores)
	os.system(command)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Workflow %%%

reference = ''
reference_folder = "%s/reference_seqs" %out_dir
if not os.path.exists(reference_folder):
	os.makedirs(reference_folder)
if args.reference_type == "user-ref-lib":
	reference = args.reference
	manage_reference = "cp %s %s" %(reference,reference_folder)
	os.system(manage_reference)
elif args.reference_type == "alignment-consensus":
	reference = create_reference_fasta(reference_folder)
for subfolder in os.listdir(reads):
#	if not "A10" in subfolder:
#		exit()
	subfolder_path = os.path.join(reads,subfolder)
	sample_folder = subfolder
	sample_id = re.sub("_clean","",sample_folder)
	if args.reference_type == "sample-specific":
		reference = create_sample_reference_fasta(reference_folder,sample_id)
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
		if args.mapper == "bwa":
			sorted_bam = mapping_bwa(forward,backward,reference,sample_id,sample_output_folder)
		elif args.mapper =="clc":
			sorted_bam = mapping_clc(forward,backward,reference,sample_id,sample_output_folder)
		if not args.keep_duplicates:
			sorted_bam = clean_with_picard(sample_output_folder,sample_id,sorted_bam)
		allele_fastas = phase_bam(sorted_bam,sample_output_folder)

## THIS iS STILL NOT WORKING PROPERLY WHEN ONLY SINGLE FILE PRESENT:

				# The following is for the case that no phased bam files were created, i.e. the individual is homozygous for all loci (happens when only looking at one locus or a very few)
		allele0 = ""
		allele1 = ""
		# testing if phasing files were created
		for file in os.listdir(allele_fastas):
			if file.endswith(".fasta"):
				if "allele_0" in file:
					allele0 = file
				if "allele_1" in file:
					allele1 = file
		if allele0 == 0:
			manage_homzygous_samples(allele_fastas,sample_id)
			os.remove(os.path.join(allele_fastas,allele0))
			os.remove(os.path.join(allele_fastas,allele1))
		print "#" * 50
join_fastas(out_dir)


#			else:
#				print "\nError: Read-files for sample %s could not be found. Please check if subfolders/sample-folders are named in this pattern: 'sampleID_clean' and if the cleaned fastq files in the sample-folder end with 'READ1.fastq' and 'READ2.fastq' respectively." %sample_id
#				raise SystemExit
#	else:
#		print "\nError: Check your folder structure. The folder given at the --reads flag must contain a separate subfolder for each sample. Please check if these sample specific subfolders are named in this pattern: 'sampleID_clean' and if the cleaned fastq files in the sample-folder end with 'READ1.fastq' and 'READ2.fastq' respectively. "
#		raise SystemExit
