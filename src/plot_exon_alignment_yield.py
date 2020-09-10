#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 13:55:09 2017

@author: tobias
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Get the alignment files and make list of loci with alignments
alignment_folder = '/Users/tobias/GitHub/seqcap_processor/data/processed/alignments/contig_alignments'
alignment_files = glob.glob(os.path.join(alignment_folder, '*.fa*'))
list_of_loci_with_alignments = [re.sub('.fasta','',al.split('/')[-1]) for al in alignment_files]

# Get the list of all exon loci from the match table (x-axis)
match_table = '/Users/tobias/GitHub/seqcap_processor/data/processed/target_contigs/match_table.txt'
matrix = pd.read_csv(match_table,sep='\t',index_col=0)
x_labels = np.array(matrix.index)
#num_x_labels = range(len(x_labels))

# Split x-axis into thirds for better readability of plots
#third_x_labels = np.split(np.matrix(x_labels), 3,axis=1)

# Create 1-dimensional matrix and fill with info which loci have alignment data 
presence_absence_df = pd.DataFrame({'loci':x_labels,'presence':0})
for locus in list_of_loci_with_alignments:
    row_index = presence_absence_df[presence_absence_df.loci == locus].index
    presence_absence_df.loc[row_index,'presence'] = 1


# Plot the data
for i in range(len(third_x_labels)):
    fig = plt.figure()
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_aspect(1)
    res = ax.imshow(np.array(third_x_labels[i]), cmap='binary')
    height,width = third_data[i].shape
    #cb = fig.colorbar(res)
    plt.xlabel('exon index',fontsize=7)
    plt.ylabel('sample index',fontsize=7)
    xlabels = list(np.array(third_num_x_labels[i])[0])
    plt.xticks(np.arange(width)[::30],xlabels[::30],fontsize=8)
    plt.yticks(fontsize=8)
    #ax.tick_params(left='off',labelleft='off')