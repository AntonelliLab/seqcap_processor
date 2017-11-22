# encoding: utf-8
'''
Extract the contigs that match the reference database
'''

from __future__ import print_function
import re
import os
import sys
import glob
import logging
import argparse
import subprocess
import numpy as np
import pandas as pd
#import pylab as plt
from phyluce.helpers import is_dir, is_file, FullPaths
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#import pdb
log = logging.getLogger(__name__)

def add_arguments(parser):
	parser.add_argument(
		'--contigs',
		required=True,
		type=is_dir,
		action=FullPaths,
		help="The directory containing the assembled contigs in fasta format."
	)
	parser.add_argument(
		'--reference',
		required=True,
		type=is_file,
		action=FullPaths,
		help="The fasta-file containign the reference sequences (probe-order-file)."
	)
	parser.add_argument(
		'--output',
		required=True,
		action=FullPaths,
		help="The directory in which to store the extracted target contigs and lastz results."
	)
	parser.add_argument(
		'--min-coverage',
		default=80,
		type=int,
		help="The minimum percent coverage required for a match [default=80]."
	)
	parser.add_argument(
		'--min-identity',
		default=80,
		type=int,
		help="The minimum percent identity required for a match [default=80]."
	)
	parser.add_argument(
		"--regex",
		type=str,
		default=".*",
		help="A regular expression to apply to the reference sequence names as tags in the output table.",
	)
	parser.add_argument(
		"--keep-duplicates",
		action='store_true',
		default=False,
		help="Use this flag in case you want to keep those contigs that span across multiple exons.",
	)


def contig_count(contig):
	"""Return a count of contigs from a fasta file"""
	return sum([1 for line in open(contig, 'rU').readlines() if line.startswith('>')])


def new_get_probe_name(header, regex):
	match = re.search(regex, header)
	#print match
	return match.groups()[0]


def contigs_matching_exons(lastz_df):
	# make a dictionary with all contig names that match a exon locus
	exon_contig_dict = {}
	contig_exon_dict = {}
	contig_orientation_dict = {}
	contig_multi_exon_dict = {}
	for row in lastz_df.iterrows():
		locus = row[1].name2
		locus_name = re.sub('\_p[0-9]* \|.*', '', locus)
		locus_name = re.sub('^>', '', locus_name)
		contig_header = row[1].name1
		#print(contig_header)
		#contig_name = re.sub('^\>([0-9]*) .*', '\\1', contig_header)
		contig_name = re.sub('^\>([^\s]*) .*', '\\1', contig_header)
		#print(contig_name)
		exon_contig_dict.setdefault(locus_name,[])
		exon_contig_dict[locus_name].append(contig_name)
		contig_exon_dict.setdefault(contig_name,[])
		contig_exon_dict[contig_name].append(locus_name)
		orientation = row[1].strand2
		contig_orientation_dict.setdefault(contig_name,orientation)
	for contig in contig_exon_dict.keys():
		if len(contig_exon_dict[contig]) > 1:
			contig_multi_exon_dict.setdefault(contig,contig_exon_dict[contig])
	return exon_contig_dict, contig_exon_dict, contig_orientation_dict, contig_multi_exon_dict


def find_duplicates(exon_contig_dict,contig_exon_dict):
	# get exons that have multiple contigs matching them
	invalid_exon_loci = []
	exons_with_multiple_hits = []
	for exon in exon_contig_dict.keys():
		if len(exon_contig_dict[exon]) > 1:
			exons_with_multiple_hits.append(exon)
			invalid_exon_loci.append(exon)
	# get exons that match on multiple contigs
	contigs_matching_multiple_exons = []
	for contig in contig_exon_dict.keys():
		if len(contig_exon_dict[contig]) > 1:
			contigs_matching_multiple_exons.append(contig)
			for exon in contig_exon_dict[contig]:
				invalid_exon_loci.append(exon)
	return invalid_exon_loci, exons_with_multiple_hits, contigs_matching_multiple_exons


def get_list_of_valid_exons_and_contigs(exon_contig_dict,duplicate_loci,exons_with_multiple_hits,contigs_matching_multiple_exon_dict,keep_duplicates_boolean,outdir):
	# summarize all exons that should be excluded form further processing (duplicates)
	if keep_duplicates_boolean:
		# then only mark the list exons_with_multiple_hits as bad exons
		invalid_exons_unique = list(set(exons_with_multiple_hits))
		valid_contigs_matching_multiple_exon_dict = {}
		for exon in contigs_matching_multiple_exon_dict.keys():
			if not exon in invalid_exons_unique:
				valid_contigs_matching_multiple_exon_dict.setdefault(exon,contigs_matching_multiple_exon_dict[exon])
		dupl_info = pd.DataFrame.from_dict(valid_contigs_matching_multiple_exon_dict, orient='index')
		dupl_info.to_csv(os.path.join(outdir,'info_contigs_matching_multiple_exons.txt'),header=False,sep="\t")
	else:
		# remove all duplicates
		invalid_exons_unique = list(set(duplicate_loci))
	print(len(invalid_exons_unique), 'possibly paralogous exons detected - excluded from processing')
	# get list of valid contig names
	valid_contig_names = []
	for exon in exon_contig_dict:
		if exon not in invalid_exons_unique:
			contig_name = exon_contig_dict[exon]
			valid_contig_names.append(contig_name[0])
	return valid_contig_names


def extract_target_contigs(sample_id,contig_sequences,valid_contig_names,contig_exon_dict,contig_orientation_dict,subfolder):
	printed_contigs_counter = 0
	# define the output file where extracted contigs will be stored
	global_match_output_name = 'extracted_target_contigs_all_samples.fasta'
	global_match_output_file = os.path.join('/'.join(subfolder.split('/')[:-1]),global_match_output_name)
	sample_match_output_name = 'extracted_target_contigs%s.fasta'%sample_id
	sample_match_output_file = os.path.join(subfolder,sample_match_output_name)
	# extract valid contigs form contig file and print to fasta file with exon-names+ sample_id as headers
	with open(global_match_output_file, "a") as out_file:
		with open(sample_match_output_file, "w") as sample_file:
			for fasta in contig_sequences:
				if fasta.id in valid_contig_names:
					orientation = contig_orientation_dict[fasta.id]
					if orientation == '-':
						seq = fasta.seq.reverse_complement()
					else:
						seq = fasta.seq
					# get the corresponding exon locus name from the dictionary
					if len(contig_exon_dict[fasta.id])>1:
						for matching_locus in contig_exon_dict[fasta.id]:
							header = '%s_%s |%s' %(matching_locus,sample_id,matching_locus)
							new_fasta = SeqRecord(seq, id=header, name='', description='')
							out_file.write(new_fasta.format('fasta'))
							sample_file.write(new_fasta.format('fasta'))
							printed_contigs_counter += 1 						
					else:		
						header = '%s_%s |%s' %(contig_exon_dict[fasta.id][0],sample_id,contig_exon_dict[fasta.id][0])
						new_fasta = SeqRecord(seq, id=header, name='', description='')
						out_file.write(new_fasta.format('fasta'))
						sample_file.write(new_fasta.format('fasta'))
						printed_contigs_counter += 1 
	return printed_contigs_counter


def main(args):
	if not os.path.isdir(args.output):
		os.makedirs(args.output)
	else:
		raise IOError("The directory {} already exists.  Please check and remove by hand.".format(args.output))
	# Get the list of exons from reference file
	pre_regex = args.regex
	regex = re.compile("^(%s)(?:.*)" %pre_regex)
	exons = set(new_get_probe_name(seq.id, regex) for seq in SeqIO.parse(open(args.reference, 'rU'), 'fasta'))
	sorted_exon_list = sorted(list(exons))
	# Get the paths to the contig fasta files for all samples
	fasta_files = glob.glob(os.path.join(args.contigs, '*.fa*'))
	sample_ids = [os.path.basename(fasta).split('.')[0] for fasta in fasta_files]
	# Create a dataframe filled with 0's
	contig_match_df = pd.DataFrame(index=sorted_exon_list,columns=sample_ids)
	for locus in sorted_exon_list:
		contig_match_df.loc[locus] = [0]*len(sample_ids)
	# Print some log screen output
	log.info("Processing contig data")
	log.info("{}".format("-" * 65))
	# Start processing by iterating through contig files (=samples)
	for contig_file in sorted(fasta_files):
		# Get the name of the sample
		critter = os.path.basename(contig_file).split('.')[0]#.replace('-', "_")
		# Make subfolder for each sample
		subfolder = os.path.join(args.output,critter)
		if not os.path.isdir(subfolder):
			os.makedirs(subfolder)		
		# Define sample-specific lastz output file
		lastz_output = os.path.join(subfolder,'%s.lastz'%critter)
		# Print some stats to screen
		total_count_of_contig = contig_count(contig_file)
		print('%s:\n'%critter,'Total contigs: %i\nSearching for contigs with matches in reference database.'%total_count_of_contig)
		# Blast the the contigs against the reference file
		with open(lastz_output, 'w') as lastz_out_file:
			lastz_command = [
				'lastz',
				'%s[multiple,nameparse=full]'%contig_file,
				'%s[nameparse=full]'%args.reference,
                '--strand=both',
                '--seed=12of19',
                '--transition',
                '--nogfextend',
                '--nochain',
                '--gap=400,30',
                '--xdrop=910',
                '--ydrop=8370',
                '--hspthresh=3000',
                '--gappedthresh=3000',
                '--noentropy',
				'--coverage=%i'%args.min_coverage,
				'--identity=%i'%args.min_identity,
				'--ambiguous=iupac',
				'--format=general:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity'
			]
			run_lastz = subprocess.Popen(lastz_command, stdout=lastz_out_file, stderr=None)
			run_lastz.communicate()
		# load the lastz matches from the previous command
		lastz_df = pd.read_csv(lastz_output,sep='\t')
		# store the data in dictionaries for convenience
		exon_contig_dict, contig_exon_dict, contig_orientation_dict, contig_multi_exon_dict = contigs_matching_exons(lastz_df)
		# mark duplicate loci
		duplicate_loci, possible_paralogous, contigs_covering_several_loci = find_duplicates(exon_contig_dict,contig_exon_dict)
		# remove duplicate loci from the list of targeted loci and contigs
		target_contigs = get_list_of_valid_exons_and_contigs(exon_contig_dict,duplicate_loci,possible_paralogous,contig_multi_exon_dict,args.keep_duplicates,subfolder)
		# load the actual contig sequences
		contig_sequences = SeqIO.parse(open(contig_file),'fasta')
		# write those contigs that match the reference library to the file
		extracted_contig_counter = extract_target_contigs(critter,contig_sequences,target_contigs,contig_exon_dict,contig_orientation_dict,subfolder)
		# Fill the extracted target contig into the dataframe
		for contig in target_contigs:
			for exon in contig_exon_dict[contig]:
				contig_match_df.loc[exon,critter] = 1		
		print('Extracted %i contigs matching reference exons\n' %extracted_contig_counter)
		log.info("{}".format("-" * 65))
	
	contig_match_df.to_csv(os.path.join(args.output,'match_table.txt'),sep='\t',index=True,encoding='utf-8')
	# Print summary stats
	table = pd.read_csv(os.path.join(args.output,'match_table.txt'), delimiter = '\t',index_col=0)
	with open(os.path.join(args.output,'summary_stats.txt'), "w") as out_file:
		out_file.write('Total number of samples: %i\nTotal number of targeted exons: %i\n\n'%(len(table.columns),len(table)))
		complete_loci_counter = 0
		for locus in table.iterrows():
			if sum(locus[1]) == len(locus[1]):
				complete_loci_counter += 1
		out_file.write('%i exons are shared by all samples.\n\n'%complete_loci_counter)
		count_list = []
		for column in table.columns:
			count_list.append(sum(table[column]))
			out_file.write('%s: %i extracted contigs\n'%(column,sum(table[column])))
		out_file.write('mean: %f stdev: %f'%(np.mean(count_list),np.std(count_list)))

#	# Make plot of match table
#	# Get the data from the df
#	sample_labels = contig_match_df.columns
#	locus_labels = np.array(contig_match_df.index)
#	data = np.matrix(contig_match_df).T
#	# Define the figure and plot to png file
#	fig, ax = plt.subplots()
#	mat = ax.imshow(data, cmap='GnBu', interpolation='nearest')
#	plt.xticks(range(data.shape[1])[::20], locus_labels[::20],fontsize=3)
#	plt.yticks([0],fontsize=0)
#	#plt.yticks(range(data.shape[0])[::3], sample_labels[::3],fontsize=3)
#	plt.xticks(rotation=90)
#	plt.xlabel('exon',fontsize=3)
#	plt.ylabel('sample',fontsize=3)
#	fig.savefig(os.path.join(args.output,'contig_exon_matrix.png'), dpi = 500)
