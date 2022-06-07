#author: Tobias Andermann, tobias.andermann@bioenv.gu.se

'''
Assemble trimmed Illumina read files (fastq)
'''

import os
import sys
import re
import glob
import shutil
import argparse
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
import time
import multiprocessing.pool
from functools import partial

# Complete path function
class CompletePath(argparse.Action):
    """give the full path of an input file/folder"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def add_arguments(parser):
    parser.add_argument(
        '--input',
        required=True,
        action=CompletePath,
        default=None,
        help='Call the folder that contains the cleaned fastq read files. The fastq files should be organized in a separate subfolder for each sample (default output of secapr clean_reads function).'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CompletePath,
        default=None,
        help='The output directory where results will be saved'
    )
    parser.add_argument(
        '--kmer',
        type=str,
        help='Set the kmer value. Provide a list of kmers for Spades, e.g. "--kmer 21,33,55". Default is 21,33,55,77,99,127. Note that Spades only accepts uneven kmer values.'
    )
    parser.add_argument(
        '--contig_length',
        type=int,
        default=200,
        help='Set the minimum contig length for the assembly. Contigs that are shorter than this threshold will be discarded.'
    )
    parser.add_argument(
        '--max_memory',
        type=str,
        help='Set the maximum memory to be used during assembly in GB. This can be necessary when working with computing nodes with limited memory or to avoid over-allocation of computing resources on clusters which can in some cases cause your assembly to be stopped or interrupted.'
    )
    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help='For parallel processing you can set the number of cores you want to run the assembly on.'
    )
    parser.add_argument(
        '--instances',
        type=int,
        default=1,
        help='How many parallel assemblies to run at a time. This will multiply the cores and max_memory arguments (for each SPAdes run). For example, max_memory=4, cores=2 and instances=4 will use 8 threads and 16 GB.'
    )

def assembly_spades(sorted_fastq_files,n_library_numbers,output_folder,id_sample,kmer,args):
    print("De-novo assembly with spades of sample %s:" %id_sample)
    command = [
        "spades.py",
        "-k",
        kmer,
        "--only-assembler"
    ]
    for i in np.arange(n_library_numbers):
        command += [
        "--pe-1",
        i+1,
        sorted_fastq_files[i*3],
        "--pe-2",
        i+1,
        sorted_fastq_files[i*3+1],
        "--pe-s",
        i+1,
        sorted_fastq_files[i*3+2],
        ]
    command += ["-o", output_folder]
    if args.cores > 1:
        command+=["--threads", args.cores]
    if args.max_memory:
        command+=["--memory", args.max_memory]
    command = list(np.array(command).astype(str))
    print ("Building contigs........")
    with open(os.path.join(output_folder, "%s_spades_screen_out.txt" %id_sample), 'w') as log_err_file:
        p = subprocess.Popen(command, stdout=log_err_file)
        p.communicate()
        p.wait()
    print(("%s assembled. Statistics are printed into %s" %(id_sample, os.path.join(output_folder, "%s_spades_screen_out.txt" %id_sample))))
    # except:
    #     print(("Could not assemble %s" %id_sample))

def count_contigs(contig_file):
    """Return a count of contigs from a fasta file"""
    return sum([1 for line in open(contig_file, 'r').readlines() if line.startswith('>')])

def remove_short_contigs(contig_file,min_length):
    fasta_content = list(SeqIO.parse(open(contig_file),'fasta'))
    new_fasta_content = []
    for record in fasta_content:
        contig_length = len(str(record.seq))
        if contig_length < min_length:
            pass
        else:
            new_fasta_content.append(record)
    SeqIO.write(new_fasta_content,contig_file, 'fasta-2line')

def get_stats_spades(contig_file,sample_id,sample_contig_count_dict):
    #contig_count_cmd = subprocess.Popen(["tail", "-n", "2", "%s/%s.fa" %('/'.join(sample_output_folder.split('/')[:-2]),sample_id)], stdout=subprocess.PIPE)
    #contig_count_pre = contig_count_cmd.communicate()[0]
    contig_count = count_contigs(contig_file)
    #contig_count = contig_count_pre.split(' ')[0].replace('>','')
    sample_contig_count_dict.setdefault(sample_id,contig_count)
    stats_df=pd.DataFrame.from_dict(sample_contig_count_dict, orient='index').reset_index()
    stats_df.columns = ['sample_id', 'total_contig_count']
    print(('#'*50))
    print(stats_df)
    return(stats_df,contig_file)
    #contig_count, header, percent, sequence = contig_count_pre.split("\t")    

def process_subfolder(pool_args):
    sample_contig_count_dict = {}
    subfolder, args = pool_args
    # Set working directory
    out_dir = os.path.join(args.output,'stats')
    # Get all the other input variables
    min_length = args.contig_length

    if args.kmer:
        kmer = str(args.kmer)
    else:
        kmer = '21,33,55,77,99,127'
    if args.max_memory:
        max_memory = args.max_memory
    else:
        max_memory = None
    sample_id = os.path.basename(subfolder)
    # Loop through each sample-folder and find read-files
    sample_output_folder = os.path.join(out_dir, sample_id)
    if not os.path.exists(sample_output_folder):
        os.makedirs(sample_output_folder)
    fastq_files = glob.glob(os.path.join(subfolder,'*.fastq.gz'))
    sorted_fastq_files = ['']*len(fastq_files)
    n_library_numbers = int(len(fastq_files) / 3)
    for element in fastq_files:
        for i in np.arange(n_library_numbers):
            if element.endswith("_%i_clean-READ1.fastq.gz"%i):
                sorted_fastq_files[i*3] = element
            elif element.endswith("_%i_clean-READ2.fastq.gz"%i):
                sorted_fastq_files[i*3+1] = element
            elif element.endswith("_%i_clean-unpaired.fastq.gz"%i):
                sorted_fastq_files[i*3+2] = element
    print(('#' * 50))
    print(("Processing sample %s" % sample_id))
    start = time.time()
    assembly_spades(sorted_fastq_files, n_library_numbers, sample_output_folder, sample_id, kmer, args)
    contig_file = os.path.join(sample_output_folder, 'contigs.fasta')
    new_contig_file = '%s/../../%s.fa' % (sample_output_folder, sample_id)
    mv_contig = "cp %s %s" % (contig_file, new_contig_file)
    os.system(mv_contig)
    contig_count_df, contig_file = get_stats_spades(new_contig_file, sample_id,sample_contig_count_dict)
    remove_short_contigs(new_contig_file, min_length)
    end = time.time()
    print('Assembled contigs for sample %s in %i minutes' % (sample_id, int(np.round((end - start) / 60))))
    return(contig_count_df)

def main(args):
    input_folder = args.input
    out_folder = args.output
    out_dir = os.path.join(out_folder,'stats')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    instances = args.instances
    subfolder_list = [subfolder for subfolder, __, __ in os.walk(input_folder) if os.path.basename(subfolder) != os.path.basename(input_folder)]
    if instances > 1:
        print(("Running assemblies in parallel as %d instances" %instances))
        pool = multiprocessing.Pool(instances)
        pool_args = [[subfolder,args] for subfolder in subfolder_list]
        contig_count_df_list = list(pool.map(partial(process_subfolder), pool_args))
        pool.close()
    else:
        contig_count_df_list = [process_subfolder([subfolder,args]) for subfolder in subfolder_list]
    #print(contig_count_df_list)
    contig_count_df = pd.concat(contig_count_df_list)
    try:
        previous_stats_df = pd.read_csv(os.path.join(input_folder,'sample_stats.txt'),sep='\t')
        counter = 0
        for index,row in previous_stats_df.iterrows():
            sample_name = str(row['sample_id'])
            if sample_name in list(contig_count_df['sample_id']):
                new_info = contig_count_df[contig_count_df['sample_id']==sample_name]['total_contig_count']
                new_value = new_info.values[0]
                new_name = new_info.name
                headers = np.array(row.index)
                old_values = row.values
                new_index = np.append(headers,new_name)
                new_values = np.append(old_values,new_value)
                if counter == 0:
                    new_values_previous = new_values
                else:
                    new_values_previous = np.vstack([new_values_previous, new_values])
                counter += 1
        new_stats_df = pd.DataFrame(data=new_values_previous,columns=new_index)
    except:
        print('No previous stats file found, creating new stats file.')
        new_stats_df = contig_count_df
    new_stats_df.to_csv(os.path.join(args.output,'sample_stats.txt'),sep="\t",index=False)





# if assembler == "abyss":
#     assembly_abyss(forward, backward, single_f, single_b, sample_output_folder, sample_id, kmer, 1, args)
#     files = glob.glob(os.path.join(sample_output_folder, '*'))
#     links = [f for f in files if os.path.islink(f)]
#     for l in links:
#         if l.endswith("-contigs.fa"):
#             contig_file = os.path.realpath(l)
#             mv_contig = "mv %s %s/../../%s.fa" % (contig_file, sample_output_folder, sample_id)
#             os.system(mv_contig)
#     # mv_cmd1 = "mv %s/%s* %s" %(home_dir,sample_id,sample_output_folder)
#     # os.system(mv_cmd1)
#     # mv_cmd2 = "mv %s/coverage.hist %s" %(home_dir,sample_output_folder)
#     # os.system(mv_cmd2)
#     contig_count_df, contig_file = get_stats_abyss(sample_output_folder, sample_id,
#                                                    sample_contig_count_dict)
#     remove_short_contigs(contig_file, min_length)

#
# def assembly_abyss(forw,backw,singlef,singleb,output_folder,id_sample,kmer,cores,args):
#     print("\nWARNING: Abyss is very memory heavy (dependend on the size of your read files and the chosen kmer value) "
#           "and it may throw an error or run painstakingly slow because it's running out of memory. "
#           "For average sized read files amounting to approx. 1 GB per sample (forward and backward reads) at kmer 35 "
#           "Abyss will require around 6GB of memory (keep that in mind when parallelizing on multiple cores).\n")
#     print(("De-novo assembly with abyss of sample %s:" %id_sample))
#     try:
#         kmer = int(kmer)
#     except:
#         quit('\n\nError: Provided kmer value could not be formatted as integer. Please provide single numeric kmer value when choosing the Abyss assembler.')
#     command = [
#         "abyss-pe",
#         "--directory={}".format(output_folder),
#         "k={}".format(kmer),
#         "j={}".format(cores),
#         'name={}'.format(id_sample),
#         'in={} {}'.format(forw,backw)
#     ]
#     if args.single_reads:
#         command.append('se={} {}'.format(singlef,singleb))
#     try:
#         print ("Building contigs........")
#         with open(os.path.join(output_folder, "%s_abyss_screen_out.txt" %id_sample), 'w') as log_err_file:
#             p = subprocess.Popen(command, stdout=log_err_file)
#             p.communicate()
#             p.wait()
#         print(("%s assembled. Statistics are printed into %s" %(id_sample, os.path.join(output_folder, "%s_abyss_screen_out.txt" %id_sample))))
#     except:
#         print(("Could not assemble %s" %id_sample))
#
# def get_stats_abyss(sample_output_folder,sample_id,sample_contig_count_dict):
#     #contig_count_cmd = subprocess.Popen(["tail", "-n", "2", "%s/%s.fa" %('/'.join(sample_output_folder.split('/')[:-2]),sample_id)], stdout=subprocess.PIPE)
#     #contig_count_pre = contig_count_cmd.communicate()[0]
#     contig_file = "%s/%s.fa" %('/'.join(sample_output_folder.split('/')[:-2]),sample_id)
#     contig_count = count_contigs(contig_file)
#     #contig_count = contig_count_pre.split(' ')[0].replace('>','')
#     sample_contig_count_dict.setdefault(sample_id,contig_count)
#     stats_df=pd.DataFrame.from_dict(sample_contig_count_dict, orient='index').reset_index()
#     stats_df.columns = ['sample', 'total_contig_count']
#     print(('#'*50))
#     print(stats_df)
#     return(stats_df,contig_file)
#     #contig_count, header, percent, sequence = contig_count_pre.split("\t")
