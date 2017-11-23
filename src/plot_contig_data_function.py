#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 18:43:28 2017

@author: tobias
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Get data for axes
#contig_input_file = '../../data/processed/target_contigs/match_table.txt'
#alignment_folder = '/Users/tobias/GitHub/seqcap_processor/data/processed/alignments/contig_alignments'


def plot_contig_yield(contig_input_file):
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
    
    
    #___________________________Plotting settings___________________________________
    height,width = data_1_contig_present.shape
    fig = plt.figure(figsize=(20,8))
    #fig.subplots_adjust(top=1, bottom=0.0, left=0.2, right=0.99)
    for i,m in enumerate(data_1_contig_present):
        ax = plt.subplot(height, 1, i+1)
        ax.tick_params(left='off',bottom='off',labelleft='off')
        # Only plot x-axis for last row
        if not i == height-1:
            ax.xaxis.set_major_formatter(plt.NullFormatter())
        #plt.axis("off")
        if data_1_y_labels[i] == 'contig alignment':
            plt.imshow(data_1_contig_present[i], aspect='auto', cmap='binary', origin='lower')
        else:
            plt.imshow(data_1_contig_present[i], aspect='auto', cmap='GnBu', origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], data_1_y_labels[i], horizontalalignment='right')
    plt.xlabel('exon index')
    #plt.colorbar()
    plt.show()
    
    
def plot_contigs_and_alignments_yield(contig_input_file,alignment_folder):
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
    alignment_files = glob.glob(os.path.join(alignment_folder, '*.fa*'))
    list_of_loci_with_alignments = [re.sub('.fasta','',al.split('/')[-1]) for al in alignment_files]
    # Create 1-dimensional matrix and fill with info which loci have alignment data 
    presence_absence_df = pd.DataFrame({'loci':x_labels,'presence':0})
    for locus in list_of_loci_with_alignments:
        row_index = presence_absence_df[presence_absence_df.loci == locus].index
        presence_absence_df.loc[row_index,'presence'] = 1
        
    data_2_contig_alignment = np.matrix(presence_absence_df.presence)
    data_2_y_labels = np.array('contig alignment')
    
    
    #_________________________Combine contig and alignment data_____________________
    contig_data_subset = np.vstack([data_1_contig_present, data_2_contig_alignment])
    y_labels_contig_data =  np.append(data_1_y_labels,data_2_y_labels)
    
    
    #___________________________Plotting settings___________________________________
    height,width = contig_data_subset.shape
    fig = plt.figure(figsize=(20,8))
    #fig.subplots_adjust(top=1, bottom=0.0, left=0.2, right=0.99)
    for i,m in enumerate(contig_data_subset):
        ax = plt.subplot(height, 1, i+1)
        ax.tick_params(left='off',bottom='off',labelleft='off')
        # Only plot x-axis for last row
        if not i == height-1:
            ax.xaxis.set_major_formatter(plt.NullFormatter())
        #plt.axis("off")
        if y_labels_contig_data[i] == 'contig alignment':
            plt.imshow(contig_data_subset[i], aspect='auto', cmap='binary', origin='lower')
        else:
            plt.imshow(contig_data_subset[i], aspect='auto', cmap='GnBu', origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], y_labels_contig_data[i], horizontalalignment='right')
    plt.xlabel('exon index')
    #plt.colorbar()
    plt.show()