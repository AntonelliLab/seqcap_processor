#author: Tobias Andermann, tobias.andermann@bioenv.gu.se

'''
Clean and trim raw Illumina read files
'''

import os
import sys
import glob
import csv
import shutil
import subprocess
import numpy as np
import pandas as pd
from secapr.helpers import FullPaths


def add_arguments(parser):
    parser.add_argument(
        '--input',
        required=True,
        action=FullPaths,
        default=None,
        help='The directory containing the .fastq or .fq files (raw read files). Files can be zipped or unzipped.'
    )
    parser.add_argument(
        '--sample_annotation_file',
        required=True,
        action=FullPaths,
        default=None,
        help='A simple comma-delimited text file containing the sample names in the first column and the name-stem of the corresponding fastq files in the second column (file should not have any column headers).'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=FullPaths,
        default=None,
        help='The output directory where the cleaned reads will be stored.'
    )
    parser.add_argument(
        '--read_min',
        type=int,
        default=400000,
        help='Set the minimum read count threshold. Any sample with fewer reads than this minimum threshold will not be processed further. Default: %(default)s'
    )
    parser.add_argument(
        '--qualified_quality_phred',
        type=int,
        default=20,
        help='Specifies how accurate the match between any adapter etc. sequence must be against a read. For more information see trimmoatic tutorial. Default: %(default)s'
    )
    parser.add_argument(
        '--unqualified_percent_limit',
        type=int,
        default=40,
        help='Set the maximum percent of low-quality nucleotides allowed. Any read with a higher percentage of unqualified (low quality) nucleotides will be discarded. Default: %(default)s'
    )
    parser.add_argument(
        '--cut_window_size',
        type=int,
        default=5,
        help='Set the size of the moving window (in nucleotides) for quality trimming. The window will start moving from the end toward the beginning of the read, applying the quality threshold set with the --cut_mean_quality flag. Default: %(default)s'
    )
    parser.add_argument(
        '--cut_mean_quality',
        type=int,
        default=20,
        help='Set quality threshold for moving window. If the mean quality across the window drops below this threshold, the nucleotides within the window are removed (cut), as well as all trailing nucleotides. Default: %(default)s'
    )
    parser.add_argument(
        '--trim_front',
        type=int,
        default=0,
        help='Remove this number of nucleotides from the beginning of each read. Default: %(default)s'
    )
    parser.add_argument(
        '--trim_tail',
        type=int,
        default=0,
        help='Remove this number of nucleotides from the end of each read. Default: %(default)s'
    )
    parser.add_argument(
        '--required_read_length',
        type=int,
        default=0,
        help='Set this value to only allow reads to pass which are equal to or longer than this threshold. Default: %(default)s'
    )
    parser.add_argument(
        '--disable_complexity_filter',
        action='store_true',
        default=False,
        help='Use this flag if you want to disable the removal of low-complexity reads, e.g. AAAAAAACCCCCCCCAAAAAAAAAAAAAAAAAAGGGGG (activated by default).'
    )
    parser.add_argument(
        '--complexity_threshold',
        type=int,
        default=10,
        help='Reads below this complexity threshold will be removed (consult fastp manual for more explanation). Default: %(default)s'
    )
    parser.add_argument(
        '--disable_poly_g_trimming',
        action='store_true',
        default=False,
        help='Use this flag if you want to disable trimming of G repeats at the end of the read. Poly-G read ends are common when working with very short fragments. (activated by default).'
    )
    parser.add_argument(
        '--poly_g_min_len',
        type=int,
        default=7,
        help='Specifies the length of the nucleotide repeat region at end of read to be trimmed. Default: %(default)s'
    )
    parser.add_argument(
        '--disable_poly_x_trimming',
        action='store_true',
        default=False,
        help='Use this flag if you want to disable trimming of nucleotide repeats of any kind at the end of the read, e.g. poly-A tails (activated by default).'
    )
    parser.add_argument(
        '--poly_x_min_len',
        type=int,
        default=7,
        help='Specifies the length of the nucleotide repeat region at end of read to be trimmed. Default: %(default)s'
    )


def longest_common_substring(s1, s2):
    m = []
    for i, letter in enumerate(list(s1)):
        if letter == list(s2)[i]:
            m.append(letter)
        else:
            break
    return ''.join(m)

# def longest_common_substring(s1, s2):
#     m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
#     longest, x_longest = 0, 0
#     for x in range(1, 1 + len(s1)):
#         for y in range(1, 1 + len(s2)):
#             if s1[x - 1] == s2[y - 1]:
#                 m[x][y] = m[x - 1][y - 1] + 1
#                 if m[x][y] > longest:
#                     longest = m[x][y]
#                     x_longest = x
#             else:
#                 m[x][y] = 0
#     return s1[x_longest - longest: x_longest]

def read_count(input):
    # Use bash command to count the read number
    number_reads = subprocess.getstatusoutput("gzip -cd %s | wc -l | awk '{ print $1/4 }' " %input)[1]
    num = int(number_reads)
    return num

def quality_trim(read1_path,read2_path,sample_id,out_dir,args,pair_index=0):
    # Get all user settings
    qualified_quality_phred = args.qualified_quality_phred
    unqualified_percent_limit = args.unqualified_percent_limit
    cut_window_size = args.cut_window_size
    cut_mean_quality = args.cut_mean_quality
    trim_front = args.trim_front
    trim_tail = args.trim_tail
    required_read_length = args.required_read_length
    disable_complexity_filter = args.disable_complexity_filter
    complexity_threshold = args.complexity_threshold
    disable_poly_g_trimming = args.disable_poly_g_trimming
    poly_g_min_len = args.poly_g_min_len
    disable_poly_x_trimming = args.disable_poly_x_trimming
    poly_x_min_len = args.poly_x_min_len
    #merge_overlapping_reads = args.merge_overlapping_reads
    #! allow for manual addition of adapters
    # Start processing
    print(('#' * 50))
    print(("Processing %s...\n" %sample_id))
    # Forward and backward read file paths
    # Names of output files
    output_sample_dir = os.path.join(out_dir,sample_id)
    if not os.path.exists(output_sample_dir):
        os.makedirs(output_sample_dir)
    outpath_r1 = os.path.join(output_sample_dir, '%s_%i_clean-READ1.fastq.gz'%(sample_id,pair_index))
    outpath_r2 = os.path.join(output_sample_dir, '%s_%i_clean-READ2.fastq.gz'%(sample_id,pair_index))
    outpath_unpaired = os.path.join(output_sample_dir, '%s_%i_clean-unpaired.fastq.gz' % (sample_id, pair_index))
    # Command for fastp trimming and cleaning
    command1 = [
        "fastp",
        "--in1",
        read1_path,
        "--in2",
        read2_path,
        "--out1",
        outpath_r1,
        "--out2",
        outpath_r2,
        "--unpaired1",
        outpath_unpaired,
        "--unpaired2",
        outpath_unpaired,
        "-h",
        os.path.join(output_sample_dir,'%s_%i_fastp.html'%(sample_id,pair_index)),
        "-j",
        os.path.join(output_sample_dir,'%s_%i_fastp.json'%(sample_id,pair_index)),
        "--qualified_quality_phred",
        qualified_quality_phred,
        "--unqualified_percent_limit",
        unqualified_percent_limit,
        "--cut_tail",
        "--cut_window_size",
        cut_window_size,
        "--cut_mean_quality",
        cut_mean_quality,
        "--trim_front1",
        trim_front,
        "--trim_tail1",
        trim_tail
    ]
    if required_read_length == 0:
        command1.append("--disable_length_filtering")
    else:
        command1 += ["--length_required", required_read_length]
    if not disable_complexity_filter:
        command1.append("--low_complexity_filter")
        command1 += ["--complexity_threshold", complexity_threshold]
    if not disable_poly_g_trimming:
        command1.append("--trim_poly_g")
        command1 += ["--poly_g_min_len", poly_g_min_len]
    if not disable_poly_x_trimming:
        command1.append("--trim_poly_x")
        command1 += ["--poly_x_min_len",poly_x_min_len]
    command1 = list(np.array(command1).astype(str))
    stats_file = os.path.join(output_sample_dir,'fastp_out_%i.txt'%pair_index)
    with open(stats_file, 'w') as log_err_file:
        p1 = subprocess.Popen(command1, stderr=log_err_file)
        p1.communicate()
        print("%s successfully cleaned and trimmed. Stats are printed into %s" %(sample_id, stats_file))
        print("#" * 50)
    return(output_sample_dir)


def main(args):
    # Set working directory
    work_dir = args.input
    out_dir = args.output
    # Print the user-set read-threshold
    read_threshold = args.read_min
    print(("\n\n[Info:] Samples with a total read-count of less than %d (forward + reverse reads) are not being processed. If required you can set a different threshold, using the --read_min flag.\n" %read_threshold))

    # Read txt-file with new sample names
    name_info_file = args.sample_annotation_file
    name_info = pd.read_csv(name_info_file,header=None)
    name_info_dict = dict(name_info.values)
    # flip keys and values in dict for easier lookup of file names
    name_info_dict = {name_info_dict[k]:k for k in name_info_dict}

    # Get all fastq files belonging to each name (can be multiple runs)
    included_extenstions = ['fastq','fq','fastq.gz','fq.gz']
    file_list = [fn for fn in sorted(os.listdir(work_dir)) if any([fn.endswith(ext) for ext in included_extenstions])]
    if len(file_list) == 0:
        sys.exit('SECAPR-ERROR: No FASTQ files were found. Check if correct path is provided for --input flag')
    sample_list = [[name_info_dict[name_stem],i] for name_stem in name_info_dict.keys() for i in file_list if i.startswith(name_stem)]
    sample_filename_dict = {}
    for sample_id,filename in sample_list:
        sample_filename_dict.setdefault(sample_id,[])
        sample_filename_dict[sample_id].append(filename)
    final_sample_filename_dict = {}
    for sample_id,filenamelist in sample_filename_dict.items():
        if len(filenamelist)>2: # if more than two files per sample, make sure to identify the pairs correctly
            shared_stems = []
            for i in filenamelist:
                longest_common_substrings = [longest_common_substring(i,j) for j in filenamelist]
                selected_longest_common_substrings = [i for i in longest_common_substrings if not any([i.endswith(ext) for ext in included_extenstions])]
                longest_common_substrings_length = np.array([len(i) for i in selected_longest_common_substrings])
                best_match = np.sort(longest_common_substrings_length)[-1] # pick the longest, because the match with itself is already filtered out in previous step
                shared_file_stem = selected_longest_common_substrings[np.where(longest_common_substrings_length==best_match)[0][0]]
                shared_stems.append(shared_file_stem)
            paired_filename_list = []
            for namestem in np.unique(shared_stems):
                paired_filename_list.append([os.path.join(work_dir,i) for i in filenamelist if i.startswith(namestem)])
        else:
            paired_filename_list = [os.path.join(work_dir,i) for i in filenamelist]
        final_sample_filename_dict.setdefault(sample_id,paired_filename_list)

    # Count total number of reads of a sample
    final_sample_list = []
    raw_read_counts = []
    for sample_id, fastq_files in final_sample_filename_dict.items():
        print("%s: Counting all reads (forward + reverse) belonging to this sample..."%sample_id)
        fastq_files = list(np.array(fastq_files).flatten())
        raw_total_read_count = sum([read_count(i) for i in fastq_files])
        print(raw_total_read_count)
        if raw_total_read_count >= read_threshold:
            final_sample_list.append(sample_id)
            raw_read_counts.append(raw_total_read_count)
        else:
            print('***Not enough reads found for sample %s. Excluding this sample from further processing. Adjust the --read_min flag to set lower minimum read count threshold.***' %sample_id )
            pass

    # Clean and trim with fastp
    #! Allow manual input of adapter sequences
    #! Allow merging of reads
    clean_reads_outdirs = {}
    for sample_id in final_sample_list:
        read_files = final_sample_filename_dict[sample_id]
        if len(np.array(read_files).flatten()) == 2:
            sample_outdir = quality_trim(read_files[0],read_files[1],sample_id,out_dir,args,pair_index=0)
        else:
            for i,read_pair in enumerate(read_files):
                sample_outdir = quality_trim(read_pair[0], read_pair[1], sample_id, out_dir, args, pair_index=i)
        clean_reads_outdirs.setdefault(sample_id,sample_outdir)

    # Count all remaining clean reads per sample
    clean_read_counts = []
    for sample_id in final_sample_list:
        sample_outdir = clean_reads_outdirs[sample_id]
        print("%s: Counting remaining clean reads (forward + reverse) belonging to this sample..."%sample_id)
        fastq_files = glob.glob(os.path.join(sample_outdir,"*.fastq.gz"))
        cleaned_total_read_count = sum([read_count(i) for i in fastq_files])
        print(cleaned_total_read_count)
        clean_read_counts.append(cleaned_total_read_count)

    stats_df = pd.DataFrame(np.vstack([final_sample_list,raw_read_counts,clean_read_counts]).T,columns = ['sample_id','raw_reads_total','cleaned_reads_total'])
    stats_df.to_csv(os.path.join(out_dir,'sample_stats.txt'),sep = '\t',index=False)

