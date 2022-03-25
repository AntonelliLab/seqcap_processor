#author: Tobias Andermann, tobias.andermann@bioenv.gu.se

'''
Clean and trim raw Illumina read files
'''

import os
import sys
import glob
import csv
import shutil
import configparser
import subprocess
import pandas as pd
from secapr.helpers import FullPaths


def add_arguments(parser):
    parser.add_argument(
        '--input',
        required=True,
        action=FullPaths,
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
        action=FullPaths,
        default=None,
        help='The output directory where results will be saved'
    )
    parser.add_argument(
        '--read_min',
        type=int,
        default=200000,
        help='Set the minimum read count threshold. Any read file containing fewer reads than this minimum threshold will not be processed further. Default: %(default)s'
    )
    parser.add_argument(
        '--index',
        type=str,
        choices=("single", "double"),
        default="single",
        help="Specify if single- or double-indexed adapters were used for the library preparation (essential information in order to interpret the control-file correctly).",
    )
    parser.add_argument(
        '--seedMismatches',
        type=int,
        default=2,
        help='Specifies the maximum mismatch count which will still allow a full match to be performed. For more information see trimmoatic tutorial. Default: %(default)s'
    )
    parser.add_argument(
        '--palindromeClipThreshold',
        type=int,
        default=30,
        help='Specifies how accurate the match between the two "adapter ligated" reads must be for PE palindrome read alignment. Default: %(default)s'
    )
    parser.add_argument(
        '--simpleClipThreshold',
        type=int,
        default=10,
        help='Specifies how accurate the match between any adapter etc. sequence must be against a read. For more information see trimmoatic tutorial. Default: %(default)s'
    )
    parser.add_argument(
        '--windowSize',
        type=int,
        default=4,
        help='Specifies the number of bases to average across. For more information see trimmoatic tutorial. Default: %(default)s'
    )
    parser.add_argument(
        '--requiredQuality',
        type=int,
        default=15,
        help='Specifies the average quality required. For more information see trimmoatic tutorial. Default: %(default)s'
    )
    parser.add_argument(
        '--leadingQuality',
        type=int,
        default=20,
        help='Specifies the minimum quality required to keep a base at the beginning of the read. For more information see trimmoatic tutorial. Default: %(default)s'
    )
    parser.add_argument(
        '--trailingQuality',
        type=int,
        default=20,
        help='Specifies the minimum quality required to keep a base at the end of a read. For more information see trimmoatic tutorial. Default: %(default)s'
    )
    parser.add_argument(
        '--cropToLength',
        type=int,
        default=250,
        help='The number of bases to keep, from the start of the read. Everything exceeding this length will be removed from the end of the read. For more information see trimmoatic tutorial. Default: %(default)s'
    )
    parser.add_argument(
        '--headCrop',
        type=int,
        default=0,
        help='The number of bases to remove from the start of the read. For more information see trimmoatic tutorial. Default: %(default)s'
    )
    parser.add_argument(
        '--minLength',
        type=int,
        default=40,
        help='Specifies the minimum length of reads to be kept. For more information see trimmoatic tutorial. Default: %(default)s'
    )

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
                print ("Reads are not single-indexed. Use '--index double' in your command.")
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
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
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
    number_reads = subprocess.getstatusoutput("gzip -cd %s | wc -l | awk '{ print $1/4 }' " %input)[1]
    num = int(number_reads)
    return num

def get_read_count_from_stats_file(stats_file):
    F = open(stats_file,'r') 
    for line in F:
        if line.startswith('Input'):
            reads_before = line.split(' ')[3]
            reads_after = line.split(' ')[6]
    return(reads_before,reads_after)

def find_fastq_pairs(name_pattern,work_dir):
    # Create a sorted (by name) list of only fastq/fq files in the input directory
    included_extenstions = ['fastq','fq','fastq.gz','fq.gz']
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
    file_info = dict(list(zip(file_list, id_of_file)))
    # Reverse the dictionary
    rev_file_info = {}
    for key, value in list(file_info.items()):
        rev_file_info.setdefault(value, set()).add(key)
    # Check if the pattern defined as key represents a full element from the name_pattern list
    for key, value in list(rev_file_info.items()):
        if key not in name_pattern:
            print(("Sample", key, "not found in control file. Skipped."))
            rev_file_info.pop(key, None)
        else:
            pass

    return rev_file_info

def quality_trim(r1,r2,sample_id,work_dir,out_dir,barcodes,conf,adapt_index,seed_mismatches,palindrome_clip_threshold,simple_clip_threshold,window_size,required_quality,leading,trailing,tail_crop,head_crop,min_length,stats_dict):
    print(('#' * 50))
    print(("Processing %s...\n" %sample_id))
    # Forward and backward read file paths
    R1 = "/".join((work_dir, r1))
    R2 = "/".join((work_dir, r2))
    # Names of output files
    output = []
    output_sample_dir = "%s/%s_clean" %(out_dir,sample_id)
    if not os.path.exists(output_sample_dir):
        os.makedirs(output_sample_dir)
    for read in ["READ1", "READ1-single", "READ2", "READ2-single"]:
        output.append(os.path.join(output_sample_dir, "%s_clean-%s.fastq.gz" %(sample_id,read)))
    # Adapters to trim
    adapter_fasta = make_adapter_fasta(sample_id,output_sample_dir,barcodes,conf,adapt_index)
    # Command for trimmomatic
    if not adapter_fasta == None:
        stats_file = os.path.join(output_sample_dir, "%s_stats.txt" %sample_id)
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
            "ILLUMINACLIP:%s:%d:%d:%d" %(adapter_fasta,seed_mismatches,palindrome_clip_threshold,simple_clip_threshold),
            "SLIDINGWINDOW:%d:%d" %(window_size,required_quality),
            "LEADING:%d" %leading,
            "TRAILING:%d" %trailing,
            "CROP:%d" %tail_crop,
            "HEADCROP:%d" %head_crop,
            "MINLEN:%d" %min_length
        ]
        with open(stats_file, 'w') as log_err_file:
            try:
                p1 = subprocess.Popen(command1, stderr=log_err_file)
                p1.communicate()
                before_reads, after_reads = get_read_count_from_stats_file(stats_file)
                stats_dict.setdefault(sample_id,[before_reads,after_reads])
                print(("%s successfully cleaned and trimmed. Stats are printed into %s" %(sample_id, os.path.join(output_sample_dir, "%s_stats.txt" %sample_id)) ))
                print(("#" * 50))
            except:
                print ("Trimmomatic was interrupted or did not start properly. You may have entered impossible values in the trimmomatic settings or trimmomatic cannot be found in the environment. Rerun again with different values for the trimmomatic flags. If that doesn't solve the problem, reinstall the secapr environment, to ensure trimmomatic being installed in the correct path.")
                sys.exit()
    else:
        print(("***********No barcodes for %s stored in config-file. Only quality trimming (no adapter trimming) will be performed***********" %sample_id))
        with open(stats_file, 'w') as log_err_file:
            command2 = [
                "trimmomatic",
                "PE",
                "-phred33",
                R1,
                R2,
                output[0],
                output[1],
                output[2],
                output[3],
                "SLIDINGWINDOW:%d:%d" %(window_size,required_quality),
                "LEADING:%d" %leading,
                "TRAILING:%d" %trailing,
                "CROP:%d" %tail_crop,
                "HEADCROP:%d" %head_crop,
                "MINLEN:%d" %min_length
            ]
            p2 = subprocess.Popen(command2, stderr=log_err_file)
            p2.communicate()
            before_reads, after_reads = get_read_count_from_stats_file(stats_file)
            stats_dict.setdefault(sample_id,[before_reads,after_reads])
            print(("%s successfully cleaned. Stats are printed into %s" %(sample_id, os.path.join(output_sample_dir, "%s_stats.txt" %sample_id)) ))
            print(("#" * 50))
    stats_df=pd.DataFrame.from_dict(stats_dict, orient='index').reset_index()
    stats_df.columns = ['sample', 'fastq_read_pairs_raw','fastq_read_pairs_cleaned']
    print(stats_df)
    return(stats_df)


def main(args):
    # Set working directory
    work_dir = args.input
    out_dir = args.output
    # Get all trimmomatic settings
    seed_mismatches = args.seedMismatches
    palindrome_clip_threshold = args.palindromeClipThreshold
    simple_clip_threshold = args.simpleClipThreshold
    window_size = args.windowSize
    required_quality = args.requiredQuality
    leading = args.leadingQuality
    trailing = args.trailingQuality
    tail_crop = args.cropToLength
    head_crop = args.headCrop
    min_length = args.minLength
    # Return the user-set or default read-threshold
    read_threshold = args.read_min
    print(("\n\n[Info:] Files with a read-count of less than %d are not being processed. If required you can set a different threshold, using the --read_min flag.\n" %read_threshold))
    adapt_index = args.index
    # Set conf as variable
    conf = configparser.ConfigParser()
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
    if len(set(delimiter)) > 1:
        quit('Multiple delimiters defined in [names] section of config file. Please choose consistent delimiter in filenames and config file!')
    # Add delimiter after the sample-name
    name_pattern = []
    for i in range(len(names_id)):
        name_pattern.append("%s%s" %(names_id[i],delimiter[i]))
    # Create the output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Find samples for which both reads exist
    read_pairs = find_fastq_pairs(name_pattern, work_dir)
    if len(read_pairs) ==0:
        sys.exit('***SECAPR-ERROR: No FASTQ files were found. Check if correct path is provided for --input flag')
    # For each pair execute the quality_trim command (trimmomatic)
    #read_count_file = open("%s/read_count_overview.txt" %out_dir, "w")
    #countlog=csv.writer(read_count_file, delimiter='\t')
    #countlog.writerow(["file","readcount"])
    stats_dict = {}
    for key, values in list(read_pairs.items()):
        if len(values) > 1:
            list_values = list(values)
            clean_list = []
            for i in range(len(list_values)):
                fq_path = "/".join((work_dir, list_values[i]))
                readcount = read_count(fq_path)
                #countlog.writerow([list_values[i],readcount])
                if readcount >= read_threshold:
                    clean_list.append(list_values[i])
                else:
                    print(('***The file %s does not contain enough reads.***' %str(list_values[i])))
                    pass
            if len(clean_list) > 1:
                r1 = ""
                r2 = ""
                for fq in clean_list:
                    pattern_r1 = ["%sR1"%str(delimiter[0]),"%sREAD1"%str(delimiter[0]),"%sRead1"%str(delimiter[0]),"%sread1"%str(delimiter[0])]
                    pattern_r2 = ["%sR2"%str(delimiter[0]),"%sREAD2"%str(delimiter[0]),"%sRead2"%str(delimiter[0]),"%sread2"%str(delimiter[0])]
                    if any(pat in fq for pat in pattern_r1):
                        r1 = fq
                    elif any(pat in fq for pat in pattern_r2):
                        r2 = fq
                    else:
                        print(('#' * 50))
                        print(("No matching read designation (R1 or R2) found for %s" %fq))
                        print(('#' * 50))
                # Remove the delimiter after the sample name in case it is part of the key
                if key.endswith(delimiter[0]):
                    clean_key = rchop(key,delimiter[0])
                    stats_df = quality_trim(r1,r2,clean_key,work_dir,out_dir,barcodes,conf,adapt_index,seed_mismatches,palindrome_clip_threshold,simple_clip_threshold,window_size,required_quality,leading,trailing,tail_crop,head_crop,min_length,stats_dict)
                else:
                    stats_df = quality_trim(r1,r2,key,work_dir,out_dir,barcodes,conf,adapt_index,seed_mismatches,palindrome_clip_threshold,simple_clip_threshold,window_size,required_quality,leading,trailing,tail_crop,head_crop,min_length,stats_dict)
    stats_df.to_csv(os.path.join(out_dir,'sample_stats.txt'),sep = '\t',index=False)

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
