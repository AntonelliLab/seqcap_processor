#!/usr/local/opt/python/bin/python
#author: Tobias Hofmann, tobiashofmann@gmx.net

'''
Clean and trim raw Illumina read files
'''



import os
import sys
import glob
import shutil
import argparse
import ConfigParser
import commands
import subprocess

from .utils import CompletePath

        
def add_arguments(parser):
    parser.add_argument(
        '--input',
        required=True,
        action=CompletePath,
        default=None,
        help='The directory containing the unzipped .fastq or .fq files (raw read files)'
    )
    parser.add_argument(
        '--config',
        required=True,
        help='A configuration file containing the adapter information and the sample names'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CompletePath,
        default=None,
        help='The output directory where results will be safed'
    )
    parser.add_argument(
        '--read_min',
        type=int,
        default=200000,
        help='Set the minimum read count threshold. Any read file containing fewer reads than this minimum threshold will not be processed further'
    )
    parser.add_argument(
        '--index',
        type=str,
        choices=("single", "double"),
        default="single",
        help="Specify if single- or double-indexed adapters were used for the library preparation (essential information in order to interpret the control-file correctly).",
    )
    '''
    parser.add_argument(
        '--trimmomatic',
        default="/usr/local/packages/anaconda2/jar/trimmomatic.jar",
        action=CompletePath,
        help='The path to the trimmomatic-0.XX.jar file.'
    )
    '''

def main(args):
    # Set working directory
    work_dir = args.input
    out_dir = args.output
    # Return the user-set or default read-threshold
    read_threshold = args.read_min
    print "\n\n[Info:] Files with a read-count of less than %d are not being processed. If required you can set a different threshold, using the --read_min flag.\n" %read_threshold
    
    adapt_index = args.index
    
    # Set conf as variable
    conf = ConfigParser.ConfigParser()
    # Read the config argument and define input as string
    conf.optionxform = str
    conf.read(args.config)
    # Call a config element
    #import ipdb; ipdb.set_trace()
    adapters = conf.items('adapters')
    barcodes = conf.items('barcodes')
    names = conf.items('names')
    
    # Read the sample name information from the config file
    names_id = []
    for element in names:
	names_id.append(element[0])
    delimiter = []
    for element in names:
        delimiter.append(element[1])
    # Add delimiter after the sample-name
    name_pattern = []
    for i in range(len(names_id)):
        name_pattern.append("%s%s" %(names_id[i],delimiter[i]))
                

    # Create the output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


    # Find samples for which both reads exist
    read_pairs = find_fastq_pairs(name_pattern, work_dir)


    # For each pair execute the quality_trim command (trimmomatic)
    for key, values in read_pairs.items():
	if len(values) > 1:
            list_values = list(values)
            clean_list = []
            for i in range(len(list_values)):
                fq_path = "/".join((work_dir, list_values[i]))
                if read_count(fq_path) >= read_threshold:
                    clean_list.append(list_values[i])
                else:
                    print "\n***The file", list_values[i], "does not contain enough reads.***\n"
                    pass
            if len(clean_list) > 1:
                r1 = ""
                r2 = ""
                for fq in clean_list:
                    pattern_r1 = ["R1","READ1","Read1","read1"]
                    pattern_r2 = ["R2","READ2","Read2","read2"]
                    if any(pat in fq for pat in pattern_r1):
                        r1 = fq
                    elif any(pat in fq for pat in pattern_r2):
                        r2 = fq
                    else:
                        print "#" * 50, "\n"
                        print "No matching read designation (R1 or R2) found for %s" %fq
                        print "#" * 50, "\n"
                # Remove the delimiter after the sample name in case it is part of the key
                if key.endswith(delimiter[0]):
                    clean_key = rchop(key,delimiter[0])
                    quality_trim(r1,r2,clean_key,work_dir,out_dir,barcodes,conf,adapt_index)
                else:
                    quality_trim(r1,r2,key,work_dir,out_dir,barcodes,conf,adapt_index)
 








                
                
def find_barcode(direction,sample_id,barcodes):
	for element in barcodes:
		tag1, tag2 = element[0].split("-")
		if direction == tag1 and sample_id == tag2:
			return element
		else:
			pass


def make_adapter_fasta(sample,sampledir,barcodes,conf,adapt_index):
	adapters = os.path.join(sampledir,"%s_adapters.fasta" %sample)
	try:
		i7_barcode = find_barcode("i7",sample,barcodes)[1]
		i7 = conf.get('adapters', 'i7')
		i7 = i7.replace("*", i7_barcode)
		i5 = conf.get('adapters', 'i5')
		if adapt_index == "single":
			try:
				i5_barcode = find_barcode("i5",sample,barcodes)[1]
			except:
				i5_barcode = None
				pass
			if not i5_barcode is None:
				print "Reads are not single-indexed. Use '--index double' in your command."
				sys.exit()
		if adapt_index == "double":
			i5_barcode = find_barcode("i5",sample,barcodes)[1]
			i5 = i5.replace("*", i5_barcode)
		with open(adapters, 'w') as outf:
			outf.write(">i5\n%s\n>i7\n%s\n" %(i5,i7))
		return adapters
	except TypeError:
		return None



def longest_common_substring(s1, s2):
	m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
	longest, x_longest = 0, 0
	for x in xrange(1, 1 + len(s1)):
		for y in xrange(1, 1 + len(s2)):
			if s1[x - 1] == s2[y - 1]:
				m[x][y] = m[x - 1][y - 1] + 1
				if m[x][y] > longest:
					longest = m[x][y]
					x_longest = x
			else:
				m[x][y] = 0
	return s1[x_longest - longest: x_longest]


def rchop(thestring, ending):
  if thestring.endswith(ending):
    return thestring[:-len(ending)]
  return thestring


def read_count(input):
	# Use bash command to count the read number
	number_reads = commands.getstatusoutput("grep -c '^+$' %s" %input)[1]
	num = int(number_reads)
	return num


def find_fastq_pairs(name_pattern,work_dir):
	# Create a sorted (by name) list of only fastq/fq files in the input directory
	included_extenstions = ['fastq','fq']
	file_list = [fn for fn in sorted(os.listdir(work_dir)) if any([fn.endswith(ext) for ext in included_extenstions])]
	# Recover the longest substring of the filename that matches an element of the sample ID list from the control file
	id_of_file = []
	for fastq in file_list:
		matches = []
		for name in name_pattern:
			matches.append(longest_common_substring(fastq,name))
		# Give the longest match
		sample_id = max(matches, key=len)
		id_of_file.append(sample_id)
	# Create a dictionary with the file names as keys and the corresponding sample IDs as values
	file_info = dict(zip(file_list, id_of_file))
	# Reverse the dictionary
	rev_file_info = {}
	for key, value in file_info.items():
		rev_file_info.setdefault(value, set()).add(key)
	# Check if the pattern defined as key represents a full element from the name_pattern list
	for key, value in rev_file_info.items():
		if key not in name_pattern:
			print "Sample", key, "not found in control file. Skipped."
			rev_file_info.pop(key, None)
		else:
			pass

	return rev_file_info


def quality_trim(r1,r2,sample_id,work_dir,out_dir,barcodes,conf,adapt_index):
	print "\n", "#" * 50
	print "Processing %s...\n" %sample_id
        # Forward and backward read file paths
	R1 = "/".join((work_dir, r1))
	R2 = "/".join((work_dir, r2))
	# Names of output files
	output = []
	output_sample_dir = "%s/%s_clean" %(out_dir,sample_id)
	if not os.path.exists(output_sample_dir):
		os.makedirs(output_sample_dir)
	for read in ["READ1", "READ1-single", "READ2", "READ2-single"]:
		output.append(os.path.join(output_sample_dir, "%s_clean-%s.fastq" %(sample_id,read)))
	# Adapters to trim
	adapter_fasta = make_adapter_fasta(sample_id,output_sample_dir,barcodes,conf,adapt_index)
	# Command for trimmomatic
	if not adapter_fasta == None:
		try:
			with open(os.path.join(output_sample_dir, "%s_stats.txt" %sample_id), 'w') as log_err_file:
				command1 = [
                                    "trimmomatic",
					"PE",
					"-phred33",
					R1,
					R2,
					output[0],
					output[1],
					output[2],
					output[3],
					"ILLUMINACLIP:%s:2:30:10" %adapter_fasta,
					"SLIDINGWINDOW:4:15",
					"LEADING:20",
					"TRAILING:20",
					"MINLEN:40"
				]
				p1 = subprocess.Popen(command1, stderr=log_err_file)
				p1.communicate()
				print "%s successfully cleaned and trimmed. Stats are printed into %s" %(sample_id, os.path.join(output_sample_dir, "%s_stats.txt" %sample_id))
				print "#" * 50, "\n"
		except:
			print "Trimmomatic was interrupted or did not start properly. Use --trimmomatic flag in command to specify the correct path to trimmomatic."
			sys.exit()
	else:
		print "***********No barcodes for %s stored in config-file. Only quality trimming (no adapter trimming) will be performed***********" %sample_id
		with open(os.path.join(output_sample_dir, "%s_stats.txt" %sample_id), 'w') as log_err_file:
			command2 = [
				"java",
				"-jar",
				trimmomatic,
				"PE",
 				"-phred33",
				R1,
				R2,
				output[0],
				output[1],
				output[2],
				output[3],
				"SLIDINGWINDOW:4:15",
				"LEADING:20",
				"TRAILING:20",
				"MINLEN:40"
			]
			p2 = subprocess.Popen(command2, stderr=log_err_file)
			p2.communicate()
			print "%s successfully cleaned. Stats are printed into %s" %(sample_id, os.path.join(output_sample_dir, "%s_stats.txt" %sample_id))
			print "#" * 50, "\n"


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
