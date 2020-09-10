#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 09:39:51 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import glob
import zipfile
import collections

def get_test_results(fastqc_log_content):
    test_results = [i for i in fastqc_log_content if i.startswith('>>')]
    test_names = [string.split('\t')[0].replace('>>','') for string in test_results if not string == '>>END_MODULE']
    test_results = [string.split('\t')[-1] for string in test_results if not string == '>>END_MODULE']
    return test_names,test_results


input_dir = '/Users/tobias/GitHub/seqcap_processor/data/raw/test_folder_quality_check/'
zip_files = glob.glob('%s*.zip'%input_dir)
sample_test_results_dict = {}
for file in zip_files:
    sample_name = file.split('/')[-1].replace('_fastqc.zip','')
    archive = zipfile.ZipFile(file,'r')
    target_file = [i for i in archive.namelist() if i.endswith('fastqc_data.txt')][0]
    fastqc_log = archive.read(target_file)
    fastqc_log_formatted = str(fastqc_log).replace('\\t','\t').split('\\n')
    labels,results = get_test_results(fastqc_log_formatted)
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


# plot the sampel overview
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
fig.savefig('/Users/tobias/Desktop/test1.pdf', dpi = 500,transparent=True)#bbox_inches='tight',

# plot the test overview
all_pass_counts = [list(col).count(0) for col in values.T]
all_warn_counts = [list(col).count(1) for col in values.T]
all_fail_counts = [list(col).count(2) for col in values.T]

barWidth=0.3
r2 = np.arange(len(all_pass_counts))
r1 = [x - barWidth for x in r2]
r3 = [x + barWidth for x in r2]

fig = plt.figure(figsize=(8,len(samples)))
plt.bar(r1, all_pass_counts, color='green', width=barWidth, edgecolor='black', label='pass')
plt.bar(r2, all_warn_counts, color='yellow', width=barWidth, edgecolor='black', label='warn')
plt.bar(r3, all_fail_counts, color='red', width=barWidth, edgecolor='black', label='fail')
plt.xticks(range(values.shape[1]), label_abbrevations)
for border in np.array(r3)+0.66*barWidth:
    plt.axvline(border,color='black',linestyle='--',alpha=0.5)
plt.yticks(range(len(samples)+1), range(len(samples)+1))
plt.xlim(0-barWidth-0.75*barWidth,)
plt.xlabel('FastQC test (abbrevated names)')
plt.ylabel('number of samples')
plt.title('FastQC results by test type')
plt.legend()
fig.savefig('/Users/tobias/Desktop/test.pdf', dpi = 500,transparent=True)#bbox_inches='tight',




