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
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')

#contig_input_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/target_contigs/match_table.txt'
#alignment_folder = '/Users/tobias/GitHub/seqcap_processor/data/processed/alignments/contig_alignments'
#read_cov_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/remapped_reads/average_cov_per_locus.txt'
#read_cov_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/selected_loci_50/overview_selected_loci.txt'

##plot_contigs_alignments_read_cov(contig_input_file,alignment_folder,read_cov_file,number_of_rows=2)
#selected_loci = plot_contigs_alignments_read_cov(contig_input_file,alignment_folder,read_cov_file,reduce=True,norm_value=10)
#selected_loci.savefig(os.path.join('/Users/tobias/GitHub/seqcap_processor/data/processed/selected_loci_50/','contig_exon_coverage_matrix_reduced.png'), dpi = 500)
#
#output_folder = '/Users/tobias/GitHub/seqcap_processor/data/processed/'
#align = general_scale_bar(2,tick_labels=['No','Yes'],x0=.1,x1=.25,plot_height=.5,plot_width=.3,font_size = 26,color1='white',color2=(0.0, 0.26666666666666666, 0.10588235294117647),height=4,width=3,plot_label='Alignment present')
#align.savefig(os.path.join(output_folder,'legend_presence_absence_alignments_green.png'), dpi = 500)
#
#contig = general_scale_bar(2,tick_labels=['No','Yes'],x0=.1,x1=.25,plot_height=.5,plot_width=.3,font_size = 26,color1='white',color2=(0.031372549019607843, 0.25098039215686274, 0.50588235294117645),height=4,width=3,plot_label='Contig present')
#contig.savefig(os.path.join(output_folder,'legend_presence_absence_contig_blue.png'), dpi = 500)
#
#
#legend = plot_heatmap_legend(0,10,font_size=26)
#legend.savefig(os.path.join('/Users/tobias/GitHub/seqcap_processor/data/processed/','legend_read_coverage.png'), dpi = 500)
#

def plot_contig_yield_linux_cluster(contig_input_file,outdir):
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
    # print a text file with the loci indeces and the corresponding loci names
    new_locus_list = x_labels
    locus_index_overview = pd.DataFrame({'loci':new_locus_list})
    locus_index_overview.to_csv(os.path.join(workdir,'locus_index_overview.txt'),sep='\t',header=False)
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
    fig.savefig(os.path.join(outdir,'contig_yield_overview.png'),bbox_inches='tight', dpi = 500)


def plot_contigs_and_alignments_yield_linux_cluster(contig_input_file,alignment_folder,outdir):
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
    alignment_files = glob.glob(os.path.join(alignment_folder, '*.fasta'))
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
    # print a text file with the loci indeces and the corresponding loci names
    new_locus_list = x_labels
    locus_index_overview = pd.DataFrame({'loci':new_locus_list})
    locus_index_overview.to_csv(os.path.join(workdir,'locus_index_overview.txt'),sep='\t',header=False)
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
    fig.savefig(os.path.join(outdir,'contig_yield_and_msas_overview.png'),bbox_inches='tight', dpi = 500)


def plot_contigs_alignments_read_cov_linux_cluster(contig_input_file,alignment_folder,read_cov_file,outdir,number_of_rows=False,font_size=12,reduce=False,string_to_remove_from_sample_names='sample_',norm_value=False):
    mpl.rcParams.update({'font.size': font_size})
    workdir = '/'.join(read_cov_file.split('/')[:-1])
    contig_matrix = pd.read_csv(contig_input_file,sep='\t',index_col=0)
    x_labels = np.array(contig_matrix.index)
    num_x_labels = range(len(x_labels))
    #______________________________1. Contig Data_____________________________________
    # Read the contig data
    data_1_contig_present = np.matrix(contig_matrix).T
    data_1_y_labels = contig_matrix.columns
    # replace substring in sample name
    data_1_y_labels = np.core.defchararray.replace(np.array(data_1_y_labels,dtype=str), string_to_remove_from_sample_names, 'contigs ')
    #_______________________________2. Contig Alignment Data__________________________
    # Get the alignment files and make list of loci with alignments
    alignment_files = glob.glob(os.path.join(alignment_folder, '*.fasta'))
    list_of_loci_with_alignments = [re.sub('.fasta','',al.split('/')[-1]) for al in alignment_files]
    # Create 1-dimensional matrix and fill with info which loci have alignment data 
    presence_absence_df = pd.DataFrame({'loci':x_labels,'presence':0})
    for locus in list_of_loci_with_alignments:
        row_index = presence_absence_df[presence_absence_df.loci == locus].index
        presence_absence_df.loc[row_index,'presence'] = 1
    data_2_contig_alignment = np.matrix(presence_absence_df.presence)
    data_2_y_labels = np.array('contig alignment')
    #_______________________________3. Reference-assembly Data__________________________
    # Get the data as pandas dataframe
    unsorted_read_cov_data = pd.read_csv(read_cov_file, sep = '\t',index_col=0)
    locus_selection=False
    if 'sum_per_locus' in unsorted_read_cov_data.columns:
        unsorted_read_cov_data = unsorted_read_cov_data.iloc[:,:-1]
        locus_selection=True
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
    data_3_y_labels = np.core.defchararray.replace(data_1_y_labels, 'contigs', 'coverage')
    #___________________________Combine all Data___________________________________
    combined_data = np.vstack([data_1_contig_present, data_2_contig_alignment,data_3_read_cov])
    tmp_combined_y_labels =  np.append(data_1_y_labels,data_2_y_labels)
    combined_y_labels = np.append(tmp_combined_y_labels,data_3_y_labels)
    height,width = combined_data.shape
    # Define the range for the heatmap, above the maximum everything will be colored in the highest color
    norm=None
    if norm_value:  
        norm = mpl.colors.Normalize(vmin=0, vmax=norm_value)
    if locus_selection and reduce:
        print('Reducing final matrix to selected loci.')
        # Reduce the data matrix to only those columns for which we have extracted values in the selected loci
        boolean = combined_data[-1]>0
        reduced_data = np.matrix([np.matrix.tolist(ind[boolean])[0] for ind in combined_data])
        # also print a text file with the loci indeces and the corresponding loci names
        new_locus_list = x_labels[np.array(boolean)[0]]
        locus_index_overview = pd.DataFrame({'loci':new_locus_list})
        locus_index_overview.to_csv(os.path.join(workdir,'locus_index_overview.txt'),sep='\t',header=False)
        # continue with plotting
        fig = plt.figure(figsize=(13.5,8))
        for i,m in enumerate(reduced_data):
            ax = plt.subplot(height, 1, i+1)
            ax.tick_params(left='off',bottom='off',labelleft='off')
            # Only plot x-axis for last row
            if not i == height-1:
                ax.xaxis.set_major_formatter(plt.NullFormatter())
            #plt.axis("off")
            if combined_y_labels[i] == 'contig alignment':
                plt.imshow(reduced_data[i], aspect='auto', cmap='Greens_r', origin='lower')
            elif 'contigs' in combined_y_labels[i]:
                plt.imshow(reduced_data[i], aspect='auto', cmap='GnBu', origin='lower')
            else:
                plt.imshow(reduced_data[i], aspect='auto', cmap='hot_r',norm=norm, origin='lower')#,clim=(0.0, 10))
            pos = list(ax.get_position().bounds)
            fig.text(pos[0] - 0.01, pos[1], combined_y_labels[i], horizontalalignment='right')
        plt.xlabel('exon index')
    else:
        # print a text file with the loci indeces and the corresponding loci names
        new_locus_list = x_labels
        locus_index_overview = pd.DataFrame({'loci':new_locus_list})
        locus_index_overview.to_csv(os.path.join(workdir,'locus_index_overview.txt'),sep='\t',header=False)
        if not number_of_rows:
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
                    plt.imshow(combined_data[i], aspect='auto', cmap='Greens', origin='lower')
                elif 'contigs' in combined_y_labels[i]:
                    plt.imshow(combined_data[i], aspect='auto', cmap='GnBu', origin='lower')
                else:
                    plt.imshow(combined_data[i], aspect='auto', cmap='hot_r',norm=norm, origin='lower')#,clim=(0.0, 10))
                pos = list(ax.get_position().bounds)
                fig.text(pos[0] - 0.01, pos[1], combined_y_labels[i], horizontalalignment='right')
            plt.xlabel('exon index')
            #plt.colorbar()
        elif number_of_rows:
            images = []
            #_______________________________Plot Split Data_____________________________
            # Split dataset for better readability
            columns_per_row = int(combined_data.shape[1]/number_of_rows)
            remainder = combined_data.shape[1]%number_of_rows
            # start and end of first row, adding the remainder to the first row
            b = 0
            n = columns_per_row+remainder
            subset_dict = {}
            # iterate through chunks
            for i in range(number_of_rows):
                data_chunk = combined_data[:,b:n]
                subset_dict.setdefault('%i,%i' %(b,n),data_chunk)
                b=n
                n+=columns_per_row
            for j in subset_dict.keys():
                split_data = subset_dict[j]
                data_range = j.split(',')
                num_x_labels = np.arange(int(data_range[0]),int(data_range[-1]))
                fig = plt.figure(figsize=(20,8))
                #fig.subplots_adjust(top=1, bottom=0.0, left=0.2, right=0.99)
                for i,m in enumerate(split_data):
                    ax = plt.subplot(height, 1, i+1)
                    ax.tick_params(left='off',bottom='off',labelleft='off')
                    # Only plot x-axis for last row
                    if not i == height-1:
                        ax.xaxis.set_major_formatter(plt.NullFormatter())
                    #plt.axis("off")
                    if combined_y_labels[i] == 'contig alignment':
                        plt.imshow(split_data[i], aspect='auto', cmap='Greens', origin='lower')
                    elif 'contigs' in combined_y_labels[i]:
                        plt.imshow(split_data[i], aspect='auto', cmap='GnBu', origin='lower')
                    else:
                        plt.imshow(split_data[i], aspect='auto', cmap='hot_r',norm=norm, origin='lower')#,clim=(0.0, 10))
                    pos = list(ax.get_position().bounds)
                    fig.text(pos[0] - 0.01, pos[1], combined_y_labels[i], horizontalalignment='right')
                # make sure to have 10 ticks on the x-axis (ensured by dividing the total lenght by 9 and using the resulting value as stepsize)
                tick_step_size = split_data[i].shape[1]/9
                # get the desired indeces of the x-values that shall carry ticks on x-axis
                xi = np.arange(0, split_data[i].shape[1],int(tick_step_size))
                # get the corresponding x-values from the num_x_labels variable (dicitonary keys)
                x = np.arange(num_x_labels[0], num_x_labels[-1], int(tick_step_size))
                plt.xticks(xi,x)
                plt.xlabel('exon index')
                fig.savefig(os.path.join(workdir,'contig_exon_coverage_matrix_%s.png'%j), dpi = 500)
                images.append(fig)
    if not number_of_rows:
        fig.savefig(os.path.join(outdir,'contig_yield_and_msas_and_readcov_overview.png'),bbox_inches='tight', dpi = 500)
    else:
        images.savefig(os.path.join(outdir,'contig_yield_and_msas_and_readcov_overview.png'),bbox_inches='tight', dpi = 500)


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
    # print a text file with the loci indeces and the corresponding loci names
    new_locus_list = x_labels
    locus_index_overview = pd.DataFrame({'loci':new_locus_list})
    locus_index_overview.to_csv(os.path.join(workdir,'locus_index_overview.txt'),sep='\t',header=False)
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
    return fig
    

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
    alignment_files = glob.glob(os.path.join(alignment_folder, '*.fasta'))
    list_of_loci_with_alignments = np.array([re.sub('.fasta','',al.split('/')[-1]) for al in alignment_files]).astype(int)
    # Create 1-dimensional matrix and fill with info which loci have alignment data 
    presence_absence_df = pd.DataFrame({'loci':x_labels,'presence':0})
    for locus in list_of_loci_with_alignments:
        row_index = presence_absence_df[presence_absence_df.loci == locus].index
        presence_absence_df.loc[row_index,'presence'] = 1
    data_2_contig_alignment = np.array(presence_absence_df.presence)
    data_2_y_labels = np.array('contig alignment')
    #_________________________Combine contig and alignment data_____________________
    contig_data_subset = np.vstack([data_1_contig_present, data_2_contig_alignment])
    y_labels_contig_data =  np.append(data_1_y_labels,data_2_y_labels)
    # print a text file with the loci indeces and the corresponding loci names
    new_locus_list = x_labels
    locus_index_overview = pd.DataFrame({'loci':new_locus_list})
    locus_index_overview.to_csv(os.path.join(workdir,'locus_index_overview.txt'),sep='\t',header=False)
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
    return fig

#contig_input_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/target_contigs/match_table.txt'
#alignment_folder = '/Users/tobias/GitHub/seqcap_processor/data/processed/alignments/contig_alignments'
#read_cov_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/selected_loci/overview_selected_loci.txt'
#norm_value=10

def plot_contigs_alignments_read_cov(contig_input_file,alignment_folder,read_cov_file,number_of_rows=False,font_size=12,reduce=False,string_to_remove_from_sample_names='sample_',norm_value=False):
    mpl.rcParams.update({'font.size': font_size})
    workdir = '/'.join(read_cov_file.split('/')[:-1])
    contig_matrix = pd.read_csv(contig_input_file,sep='\t',index_col=0)
    x_labels = np.array(contig_matrix.index)
    num_x_labels = range(len(x_labels))
    #______________________________1. Contig Data_____________________________________
    # Read the contig data
    data_1_contig_present = np.array(contig_matrix).T
    data_1_y_labels = contig_matrix.columns
    # replace substring in sample name
    data_1_y_labels = np.array([i.replace(string_to_remove_from_sample_names,'')+'_contigs' for i in data_1_y_labels])
    #_______________________________2. Contig Alignment Data__________________________
    # Get the alignment files and make list of loci with alignments
    alignment_files = glob.glob(os.path.join(alignment_folder, '*.fasta'))
    list_of_loci_with_alignments = np.array([re.sub('.fasta','',al.split('/')[-1]) for al in alignment_files]).astype(int)
    # Create 1-dimensional matrix and fill with info which loci have alignment data 
    presence_absence_df = pd.DataFrame({'loci':x_labels,'presence':0})
    for locus in list_of_loci_with_alignments:
        row_index = presence_absence_df[presence_absence_df.loci == locus].index
        presence_absence_df.loc[row_index,'presence'] = 1
    data_2_contig_alignment = np.array(presence_absence_df.presence)
    data_2_y_labels = np.array('contig alignment')
    #_______________________________3. Reference-assembly Data__________________________
    # Get the data as pandas dataframe
    unsorted_read_cov_data = pd.read_csv(read_cov_file, sep = '\t',index_col=0)
    locus_selection=False
    if 'sum_per_locus' in unsorted_read_cov_data.columns:
        unsorted_read_cov_data = unsorted_read_cov_data.iloc[:,:-1]
        locus_selection=True
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
    data_3_y_labels = np.core.defchararray.replace(data_1_y_labels, 'contigs', 'coverage')
    #___________________________Combine all Data___________________________________
    combined_data = np.vstack([data_1_contig_present, data_2_contig_alignment,data_3_read_cov])
    tmp_combined_y_labels =  np.append(data_1_y_labels,data_2_y_labels)
    combined_y_labels = np.append(tmp_combined_y_labels,data_3_y_labels)
    height,width = combined_data.shape
    # Define the range for the heatmap, above the maximum everything will be colored in the highest color
    norm=None
    if norm_value:  
        norm = mpl.colors.Normalize(vmin=0, vmax=norm_value)
    if locus_selection and reduce:
        print('Reducing final matrix to selected loci.')
        # Reduce the data matrix to only those columns for which we have extracted values in the selected loci
        boolean = combined_data[-1]>0
        reduced_data = np.matrix([np.matrix.tolist(ind[boolean])[0] for ind in combined_data])
        # also print a text file with the loci indeces and the corresponding loci names
        new_locus_list = x_labels[np.array(boolean)[0]]
        locus_index_overview = pd.DataFrame({'loci':new_locus_list})
        locus_index_overview.to_csv(os.path.join(workdir,'locus_index_overview.txt'),sep='\t',header=False)
        # continue with plotting
        fig = plt.figure(figsize=(13.5,8))
        for i,m in enumerate(reduced_data):
            ax = plt.subplot(height, 1, i+1)
            ax.tick_params(left='off',bottom='off',labelleft='off')
            # Only plot x-axis for last row
            if not i == height-1:
                ax.xaxis.set_major_formatter(plt.NullFormatter())
            #plt.axis("off")
            if combined_y_labels[i] == 'contig alignment':
                plt.imshow(reduced_data[i], aspect='auto', cmap='Greens_r', origin='lower')
            elif 'contigs' in combined_y_labels[i]:
                plt.imshow(reduced_data[i], aspect='auto', cmap='GnBu', origin='lower')
            else:
                plt.imshow(reduced_data[i], aspect='auto', cmap='hot_r',norm=norm, origin='lower')#,clim=(0.0, 10))
            pos = list(ax.get_position().bounds)
            fig.text(pos[0] - 0.01, pos[1], combined_y_labels[i], horizontalalignment='right')
        plt.xlabel('exon index')
    else:
        # print a text file with the loci indeces and the corresponding loci names
        new_locus_list = x_labels
        locus_index_overview = pd.DataFrame({'loci':new_locus_list})
        locus_index_overview.to_csv(os.path.join(workdir,'locus_index_overview.txt'),sep='\t',header=False)
        if not number_of_rows:
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
                    plt.imshow(combined_data[i], aspect='auto', cmap='Greens', origin='lower')
                elif 'contigs' in combined_y_labels[i]:
                    plt.imshow(combined_data[i], aspect='auto', cmap='GnBu', origin='lower')
                else:
                    plt.imshow(combined_data[i], aspect='auto', cmap='hot_r',norm=norm, origin='lower')#,clim=(0.0, 10))
                pos = list(ax.get_position().bounds)
                fig.text(pos[0] - 0.01, pos[1], combined_y_labels[i], horizontalalignment='right')
            plt.xlabel('exon index')
            #plt.colorbar()
        elif number_of_rows:
            images = []
            #_______________________________Plot Split Data_____________________________
            # Split dataset for better readability
            columns_per_row = int(combined_data.shape[1]/number_of_rows)
            remainder = combined_data.shape[1]%number_of_rows
            # start and end of first row, adding the remainder to the first row
            b = 0
            n = columns_per_row+remainder
            subset_dict = {}
            # iterate through chunks
            for i in range(number_of_rows):
                data_chunk = combined_data[:,b:n]
                subset_dict.setdefault('%i,%i' %(b,n),data_chunk)
                b=n
                n+=columns_per_row
            for j in subset_dict.keys():
                split_data = subset_dict[j]
                data_range = j.split(',')
                num_x_labels = np.arange(int(data_range[0]),int(data_range[-1]))
                fig = plt.figure(figsize=(20,8))
                #fig.subplots_adjust(top=1, bottom=0.0, left=0.2, right=0.99)
                for i,m in enumerate(split_data):
                    ax = plt.subplot(height, 1, i+1)
                    ax.tick_params(left='off',bottom='off',labelleft='off')
                    # Only plot x-axis for last row
                    if not i == height-1:
                        ax.xaxis.set_major_formatter(plt.NullFormatter())
                    #plt.axis("off")
                    if combined_y_labels[i] == 'contig alignment':
                        plt.imshow(split_data[i], aspect='auto', cmap='Greens', origin='lower')
                    elif 'contigs' in combined_y_labels[i]:
                        plt.imshow(split_data[i], aspect='auto', cmap='GnBu', origin='lower')
                    else:
                        plt.imshow(split_data[i], aspect='auto', cmap='hot_r',norm=norm, origin='lower')#,clim=(0.0, 10))
                    pos = list(ax.get_position().bounds)
                    fig.text(pos[0] - 0.01, pos[1], combined_y_labels[i], horizontalalignment='right')
                # make sure to have 10 ticks on the x-axis (ensured by dividing the total lenght by 9 and using the resulting value as stepsize)
                tick_step_size = split_data[i].shape[1]/9
                # get the desired indeces of the x-values that shall carry ticks on x-axis
                xi = np.arange(0, split_data[i].shape[1],int(tick_step_size))
                # get the corresponding x-values from the num_x_labels variable (dicitonary keys)
                x = np.arange(num_x_labels[0], num_x_labels[-1], int(tick_step_size))
                plt.xticks(xi,x)
                plt.xlabel('exon index')
                fig.savefig(os.path.join(workdir,'contig_exon_coverage_matrix_%s.png'%j), dpi = 500)
                images.append(fig)
    if not number_of_rows:
        return fig
    else:
        return images

def plot_heatmap_legend(min_value,max_value,font_size=26,width=3,height=8):
    #________________________________Plot Legend___________________________________
    # Make a figure and axes with dimensions as desired.
    mpl.rcParams.update({'font.size': font_size})
    fig = plt.figure(figsize=(width, height))
    # the values stand for [x0,x1,width,height] --> all in relation to total size as given by 'figsize='
    ax1 = fig.add_axes([0.1,0.05,.3,.9])
    # Set the colormap and norm to correspond to the data for which
    # the colorbar will be used.
    cmap = mpl.cm.hot_r
    norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
    # plot a basic continuous colorbar with ticks and labels.
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                    norm=norm,
                                    orientation='vertical')
    cb1.set_label('Read coverage')
    cb1.ax.set_yticklabels(['0','2','4','6','8','>10'])
    #plt.show()
    return fig


def plot_presence_absence_legend(font_size=22,width=2,height=2):
    mpl.rcParams.update({'font.size': font_size})
    fig = plt.figure(figsize=(width, height))
    # the values stand for [x0,x1,width,height] --> all in relation to total size as given by 'figsize='
    ax1 = fig.add_axes([0.1,0.1,.4,.8])
    # Set the colormap and norm to correspond to the data for which
    cmap = mpl.cm.binary
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # define black and white manually, make beginning to center black and center to end white
    center = int(len(cmaplist)/2)
    cmaplist[0:center + 1] = [(1.0,1.0,1.0,1.0)] * ((center + 1) - 0)
    cmaplist[center+1:] = [(0.0,0.0,0.0,1.0)]*(len(cmaplist)+1-center)
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist)
    bounds = np.linspace(0,1,3)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    # plot a basic continuous colorbar with ticks and labels.
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                    norm=norm,
                                    orientation='vertical',spacing='proportional', ticks=[0.25,0.75], boundaries=bounds)
    cb1.ax.set_yticklabels(['Absence','Presence'])
    #cb1.set_label('Presence/Absence')
    #fig.savefig(os.path.join(output_folder,'legend_presence_absence_alignments.png'), dpi = 500)

def general_scale_bar(n_categories,tick_labels='auto',font_size=12,x0=0.1,x1=0.1,plot_height=0.8,plot_width=0.4,colormap=False,color1='yellow',color2='red',height=2,width=2,plot_label=False,norm_max=1):
#general_scale_bar(2,tick_labels=['Absence','Presence'],color1='white',color2='darkgreen',height=2,width=2)
    if tick_labels == 'auto':
        labels = np.arange(0,n_categories)
    else:
        labels=tick_labels
    mpl.rcParams.update({'font.size': font_size})
    fig = plt.figure(figsize=(width, height))
    # the values stand for [x0,x1,width,height] --> all in relation to total size as given by 'figsize='
    ax = fig.add_axes([x0,x1,plot_width,plot_height])
    # define a colormap by providing the beginning and end color and number of categories between these colors
    cmap = mpl.colors.LinearSegmentedColormap.from_list("my_colormap", (color1, color2), N=n_categories, gamma=1.0)
    bounds = np.linspace(0,1,n_categories+1)
    # the ticks should be in the cneter between any two bounds, so move all 1/2 unit up
    ticks = bounds[:-1]+(1/n_categories/2)
    norm = mpl.colors.Normalize(0, norm_max)
    if colormap:
        cb = mpl.colorbar.ColorbarBase(ax, cmap=colormap,norm=norm,
                                       orientation='vertical',spacing='proportional', ticks=ticks, boundaries=bounds)
    else:
        cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,norm=norm,
                                       orientation='vertical',spacing='proportional', ticks=ticks, boundaries=bounds)
    cb.ax.set_yticklabels(labels)
    if plot_label:
        cb.set_label(plot_label)
    return fig


def slice_data_columns_per_row(combined_data,columns_per_row):
    number_of_rows = int(combined_data.shape[1]/columns_per_row)
    b = 0
    n = columns_per_row
    subset_dict = {}
    # iterate through chunks
    for i in range(number_of_rows):
        data_chunk = combined_data[:,b:n]
        subset_dict.setdefault('%i,%i' %(b,n),data_chunk)
        b=n
        n+=columns_per_row
    # take care of the remainder, if there is some left
    if b<combined_data.shape[1]:
        n=combined_data.shape[1]
        data_chunk = combined_data[:,b:n]
        subset_dict.setdefault('%i,%i' %(b,n),data_chunk)


def plot_multiple_plots_in_one_ipython_cell():
    import glob
    import os
    import numpy as np
    import scipy.misc
    import matplotlib.pyplot as plt
    #%matplotlib inline

    # get all subplot files:
    work_dir = "../../data/processed/remapped_reads/"
    images = []
    for file in glob.glob(os.path.join(work_dir,"*.png")):
        images.append(file)

    # convert images into data objects
    img_data = []
    for img in images:
        img_data.append(scipy.misc.imread(img))

    # plot all images in one window
    plt.figure(figsize=(200,100))
    columns = 1
    for i, image in enumerate(np.array(img_data)):
        plt.subplot(len(images) / columns + 1, columns, i + 1)
        plt.axis('off')
        plt.imshow(image)