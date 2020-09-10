#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 11:40:16 2017

@author: tobias
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Get data for axes
contig_input_file = '/Users/tobias/Desktop/target_contigs/match_table.txt'
workdir = '/'.join(contig_input_file.split('/')[:-1])
contig_matrix = pd.read_csv(contig_input_file,sep='\t',index_col=0)

x_labels = np.array(contig_matrix.index)
num_x_labels = range(len(x_labels))


#______________________________Contig Data_____________________________________
# Read the contig data
data_1_contig_present = np.matrix(contig_matrix).T
data_1_y_labels = contig_matrix.columns
# replace substring in sample name
data_1_y_labels = np.core.defchararray.replace(np.array(data_1_y_labels,dtype=str), 'sample_', 'contigs ')


#_______________________________Contig Alignment Data__________________________
# Get the alignment files and make list of loci with alignments
alignment_folder = '/Users/tobias/Desktop/target_contigs/msa_alignments'
alignment_files = glob.glob(os.path.join(alignment_folder, '*.fa*'))
list_of_loci_with_alignments = [re.sub('.fasta','',al.split('/')[-1]) for al in alignment_files]
# Create 1-dimensional matrix and fill with info which loci have alignment data 
presence_absence_df = pd.DataFrame({'loci':x_labels,'presence':0})
for locus in list_of_loci_with_alignments:
    row_index = presence_absence_df[presence_absence_df.loci == locus].index
    presence_absence_df.loc[row_index,'presence'] = 1
    
data_2_contig_alignment = np.matrix(presence_absence_df.presence)
data_2_y_labels = np.array('contig alignment')


#_______________________________Reference-assembly Data__________________________
# Get the data as pandas dataframe
read_cov_file = '/Users/tobias/GitHub/target_contigs/remapped_reads/average_cov_per_locus.txt'
unsorted_read_cov_data = pd.read_csv(read_cov_file, sep = '\t',index_col=0)
# sort columns in df
temp_read_cov_data = unsorted_read_cov_data[sorted(unsorted_read_cov_data.columns)].sort_index()
# add row of 0's for all missing loci
loci_in_df = list(temp_read_cov_data.index)
for locus in list(x_labels):
    if locus not in loci_in_df:
        temp_read_cov_data.loc[locus] = [0.0]*len(temp_read_cov_data.columns)
# sort by index again
read_cov_data = temp_read_cov_data.sort_index()

# turn df into matrix
data_3_read_cov = np.matrix(read_cov_data).T
# lets use the same labels as for the contig data
data_3_y_labels = np.core.defchararray.replace(data_1_y_labels, 'contigs ', 'coverage ')



#___________________________Combine all Data___________________________________
combined_data = np.vstack([data_1_contig_present, data_2_contig_alignment,data_3_read_cov])
tmp_combined_y_labels =  np.append(data_1_y_labels,data_2_y_labels)
combined_y_labels = np.append(tmp_combined_y_labels,data_3_y_labels)

height,width = combined_data.shape


#_______________________________Plot Combined Data_____________________________
fig = plt.figure(figsize=(20,8))
#fig.subplots_adjust(top=1, bottom=0.0, left=0.2, right=0.99)
for i,m in enumerate(combined_data):
    ax = plt.subplot(height, 1, i+1)
    ax.tick_params(left='off',bottom='off',labelleft='off')
    # Only plot x-axis for last row
    if not i == height-1:
        ax.xaxis.set_major_formatter(plt.NullFormatter())
    #plt.axis("off")
    if combined_y_labels[i] == 'contig alignment':
        plt.imshow(combined_data[i], aspect='auto', cmap='binary', origin='lower')
    elif 'contigs' in combined_y_labels[i]:
        plt.imshow(combined_data[i], aspect='auto', cmap='GnBu', origin='lower')
    else:
        plt.imshow(combined_data[i], aspect='auto', cmap='hot_r', origin='lower',clim=(0.0, 10))
    pos = list(ax.get_position().bounds)
    fig.text(pos[0] - 0.01, pos[1], combined_y_labels[i], horizontalalignment='right')
plt.xlabel('exon index')
#plt.colorbar()
plt.show()
fig.savefig(os.path.join(workdir,'exon_yield_all_datatypes.png'), dpi = 500)



#________________________________Plot Legend___________________________________
# Make a figure and axes with dimensions as desired.
fig = plt.figure(figsize=(1, 8))
# the values stand for [x0,x1,width,height] --> all in relation to total size as given by 'figsize='
ax1 = fig.add_axes([0.1,0.05,.4,.9])
# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.hot_r
norm = mpl.colors.Normalize(vmin=0, vmax=10)
# plot a basic continuous colorbar with ticks and labels.
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
cb1.set_label('Read coverage')
#plt.show()
fig.savefig(os.path.join(workdir,'legend.png'), dpi = 500)
