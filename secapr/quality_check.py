# encoding: utf-8
#author: Tobias Andermann, tobias.andermann@bioenv.gu.se
"""
This script runs a fastqc test on all fastq samples in a user-provided folder and creates an overview plot,
"""

import matplotlib
matplotlib.use('Agg')

import os
import sys
import glob
import fnmatch
import shutil
import ConfigParser
from .utils import CompletePath
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import zipfile
import collections


def add_arguments(parser):
    parser.add_argument(
        '--input',
        required=True,
        action=CompletePath,
        default=None,
        help='The directory containing fastq files'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CompletePath,
        default=None,
        help='The output directory where quality-test results will be saved'
    )

def get_test_results(fastqc_log_content):
    test_results = [i for i in fastqc_log_content if i.startswith('>>')]
    test_names = [string.split('\t')[0].replace('>>','') for string in test_results if not string == '>>END_MODULE']
    test_results = [string.split('\t')[-1] for string in test_results if not string == '>>END_MODULE']
    return test_names,test_results

def plot_fastqc_results(fastqc_out_folder):
    zip_files = []
    for root, dirnames, filenames in os.walk(fastqc_out_folder):
        for filename in fnmatch.filter(filenames, '*.zip'):
            zip_files.append(os.path.join(root, filename))
    sample_test_results_dict = {}
    for file in zip_files:
        sample_name = file.split('/')[-1].replace('_fastqc.zip','')
        archive = zipfile.ZipFile(file,'r')
        target_file = [i for i in archive.namelist() if i.endswith('fastqc_data.txt')][0]
        fastqc_log = archive.read(target_file)
        fastqc_log_formatted = str(fastqc_log).replace('\\t','\t').split('\n')
        labels,results = get_test_results(fastqc_log_formatted)
        #print(results)
        num_results = [0 if i == 'pass' else i for i in results]
        num_results = [1 if i == 'warn' else i for i in num_results]
        num_results = [2 if i == 'fail' else i for i in num_results]
        sample_test_results_dict[sample_name] = num_results

    label_abbrevations = []
    for i in labels:
        split_string = i.split(' ')
        abbrevation = []
        for j in split_string:
            letter = j[0]
            abbrevation.append(letter)
        abbrevation = ''.join(abbrevation)
        label_abbrevations.append(abbrevation)
    # plot the sample overview
    ordered_dict = collections.OrderedDict(sorted(sample_test_results_dict.items()))
    samples = list(ordered_dict.keys())
    values = np.array(list(ordered_dict.values()))

    fig = plt.figure(figsize=(8,len(samples)))
    plt.imshow(values, interpolation='nearest', cmap=colors.ListedColormap(['green','yellow','red']))
    plt.yticks(range(values.shape[0]), samples)
    plt.xticks(range(values.shape[1]), label_abbrevations)
    plt.xlabel('FastQC test (abbrevated names)')
    plt.ylabel('Sample name')
    plt.title('FastQC results by sample')
    fig.savefig(os.path.join(fastqc_out_folder,'quality_summary_all_samples_1.pdf'), dpi = 500,transparent=True,bbox_inches='tight')

    # plot the test overview
    all_pass_counts = [list(col).count(0) for col in values.T]
    all_warn_counts = [list(col).count(1) for col in values.T]
    all_fail_counts = [list(col).count(2) for col in values.T]

    barWidth=0.3
    r2 = np.arange(len(all_pass_counts))
    r1 = [x - barWidth for x in r2]
    r3 = [x + barWidth for x in r2]

    fig = plt.figure(figsize=(8,1+len(samples)/8))
    plt.bar(r1, all_pass_counts, color='green', width=barWidth, edgecolor='black', label='pass')
    plt.bar(r2, all_warn_counts, color='yellow', width=barWidth, edgecolor='black', label='warn')
    plt.bar(r3, all_fail_counts, color='red', width=barWidth, edgecolor='black', label='fail')
    plt.xticks(range(values.shape[1]), label_abbrevations)
    for border in np.array(r3)+0.66*barWidth:
        plt.axvline(border,color='black',linestyle='--',alpha=0.5)
    #plt.yticks(range(len(samples)+1), range(len(samples)+1))
    plt.xlim(0-barWidth-0.75*barWidth,)
    plt.xlabel('FastQC test (abbrevated names)')
    plt.ylabel('number of samples')
    plt.title('FastQC results by test type')
    #plt.legend()
    fig.savefig(os.path.join(fastqc_out_folder,'quality_summary_all_samples_2.pdf'), dpi = 500,transparent=True,bbox_inches='tight')

def main(args):
    # Set working directory
    out_folder = args.output
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # Get list of all fastq-files
    input_folder = args.input
    matches = []
    for root, dirnames, filenames in os.walk(input_folder):
        for filename in fnmatch.filter(filenames, '*.fastq'):
            if not filename.endswith('-single.fastq'):
                matches.append(os.path.join(root, filename))
    if len(matches) == 0:
        print('No files with the ending .fastq found in input folder. Please check path and ensure that all readfiles are unzipped and have the filending ".fastq"')
        sys.exit()
    fastq_df = pd.DataFrame(index=np.arange(0,len(matches)), columns=['filepaths'])
    fastq_df['filepaths'] = matches
    fastq_list_path = os.path.join(out_folder,'fastq_file_list.txt')
    fastq_df.to_csv(fastq_list_path,index=False,header=False,sep='\t')

    # run FASTQC
    fastqc_cmd = [
        'fastqc -o %s -f fastq $(cat %s)' %(out_folder,fastq_list_path)
    ]
    with open(os.path.join(out_folder, "fastqc_screen_out.txt"), 'w') as log_err_file:
        p = subprocess.Popen(fastqc_cmd, stdout=log_err_file, stderr=log_err_file, shell=True)
        p.communicate()

    plot_fastqc_results(out_folder)
