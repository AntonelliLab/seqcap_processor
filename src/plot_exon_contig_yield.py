#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 11:40:16 2017

@author: tobias
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read the input data
input_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/target_contigs/match_table.txt'
workdir = '/'.join(input_file.split('/')[:-1])
matrix = pd.read_csv(input_file,sep='\t',index_col=0)
data = np.matrix(matrix).T
y_labels = matrix.columns
x_labels = np.array(matrix.index)
num_x_labels = range(len(x_labels))

# Split dataset into thirds for better readability
third_data = np.split(data, 3,axis=1)
third_x_labels = np.split(np.matrix(x_labels), 3,axis=1)
third_num_x_labels = np.split(np.matrix(num_x_labels),3,axis=1)

# Plot the matrices
for i in range(len(third_data)):
    fig = plt.figure()
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_aspect(1)
    res = ax.imshow(np.array(third_data[i]), cmap='GnBu')
    height,width = third_data[i].shape
    #cb = fig.colorbar(res)
    plt.xlabel('exon index',fontsize=7)
    plt.ylabel('sample index',fontsize=7)
    xlabels = list(np.array(third_num_x_labels[i])[0])
    plt.xticks(np.arange(width)[::30],xlabels[::30],fontsize=8)
    plt.yticks(fontsize=8)
    #ax.tick_params(left='off',labelleft='off')
    fig.savefig(os.path.join(workdir,'contig_exon_matrix_%i.png'%i), dpi = 500)

# Write overview of exon indeces
key_to_exon_index = pd.DataFrame({'index':num_x_labels,'locus_name': x_labels})
key_to_exon_index.to_csv(os.path.join(workdir,'key_to_exon_index.txt'),index=False,sep='\t')

# Write overview of sample indeces
key_to_sample_index = pd.DataFrame({'index':range(len(y_labels)),'sample_ID': y_labels})
key_to_sample_index.to_csv(os.path.join(workdir,'key_to_sample_index.txt'),index=False,sep='\t')
