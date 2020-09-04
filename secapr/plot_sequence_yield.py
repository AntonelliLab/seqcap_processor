# encoding: utf-8
'''
Plot overview of extracted sequences
'''


import logging
import argparse
from secapr.helpers import is_dir, is_file, FullPaths
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')

#import pdb
log = logging.getLogger(__name__)

def add_arguments(parser):
    parser.add_argument(
        '--extracted_contigs',
        required=True,
        type=is_dir,
        action=FullPaths,
        help="The directory containing the extracted target contigs (output from find_target_contigs function)."
    )
    parser.add_argument(
        '--alignments',
        required=False,
        type=is_dir,
        action=FullPaths,
        help="The directory containing the contig alignments. Provide this path if you want to add a line to the plot showing for which loci alignments could be created."
    )
    parser.add_argument(
        '--read_cov',
        required=False,
        type=is_dir,
        action=FullPaths,
        help="The directory containing the reference assembly results. Provide this path if you want to display the read coverage for each locus and sample."
    )
    parser.add_argument(
        '--coverage_norm',
        default=10,
        type=int,
        help="Here you can adjust the color scale of the read-coverage plot. This value will define the maximum of the color scale, e.g. if set to '10' read-coverage of 10 and above will be colored black, while everything below (0-10) will be stretched across the color spectrum from yellow, red to black."
    )
    parser.add_argument(
        '--output',
        required=True,
        action=FullPaths,
        help="The directory in which to store the plots."
    )


def plot_contig_yield(contig_input_file,outdir,alignment_folder=False,read_cov_file=False,norm_value = None):
    workdir = '/'.join(contig_input_file.split('/')[:-1])
    contig_matrix = pd.read_csv(contig_input_file,sep='\t',index_col=0)
    x_labels = np.array(contig_matrix.index)
    num_x_labels = list(range(len(x_labels)))
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
    combined_data = data_1_contig_present
    combined_y_labels =  data_1_y_labels

    #_______________________________Contig Alignment Data__________________________
    # Get the alignment files and make list of loci with alignments
    if alignment_folder:
        alignment_files = glob.glob(os.path.join(alignment_folder, '*.fasta'))
        list_of_loci_with_alignments = [re.sub('.fasta','',al.split('/')[-1]) for al in alignment_files]
        # Create 1-dimensional matrix and fill with info which loci have alignment data 
        presence_absence_df = pd.DataFrame({'loci':x_labels,'presence':0})
        for locus in list_of_loci_with_alignments:
            row_index = presence_absence_df[presence_absence_df.loci == int(locus)].index
            presence_absence_df.loc[row_index,'presence'] = 1
        data_2_contig_alignment = np.matrix(presence_absence_df.presence)
        data_2_y_labels = np.array('Contig alignment')

        #___________________________Combine all Data___________________________________
        combined_data = np.vstack([data_1_contig_present, data_2_contig_alignment])
        combined_y_labels =  np.append(data_1_y_labels,data_2_y_labels)

    #_______________________________Reference-assembly Data__________________________
    # Get the data as pandas dataframe
    if read_cov_file:
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
        data_3_y_labels = np.core.defchararray.replace(data_1_y_labels, 'contigs', 'coverage')

        #___________________________Combine all Data___________________________________
        combined_data = np.vstack([data_1_contig_present, data_2_contig_alignment,data_3_read_cov])
        tmp_combined_y_labels =  np.append(data_1_y_labels,data_2_y_labels)
        combined_y_labels = np.append(tmp_combined_y_labels,data_3_y_labels)

    #___________________________Plotting settings___________________________________
    height,width = combined_data.shape
    norm=None
    if norm_value:  
        norm = mpl.colors.Normalize(vmin=0, vmax=norm_value)
    switch = 'off'

    fig, axes2d = plt.subplots(nrows=height, ncols=1,sharex=True, sharey=True,figsize=(width/40,height/2))
    for i, ax in enumerate(axes2d):
        ax = plt.subplot(height, 1, i+1)
        ax.tick_params(left=False,bottom=False,labelleft=True)
        # Only plot x-axis for last row
        if not i == height-1:
            ax.xaxis.set_major_formatter(plt.NullFormatter())
        #plt.axis("off")
        if combined_y_labels[i] == 'Contig alignment':
            ax.imshow(combined_data[i], aspect='auto', cmap='Greens', origin='lower')
            switch = 'on'
        else:
            if switch == 'off':
                ax.imshow(combined_data[i], aspect='auto', cmap='GnBu', origin='lower') 
            else :
                ax.imshow(combined_data[i], aspect='auto', cmap='hot_r',norm=norm, origin='lower')#,clim=(0.0, 10))
        plt.yticks([0],[combined_y_labels[i]])
    ax.set_xlabel('Exon index')

    if read_cov_file:
        fig.add_subplot(211, frameon=False)
        # hide tick and tick label of the big axis
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.ylabel("Contig present (yes/no)",labelpad=50)
        fig.add_subplot(212, frameon=False)
        # hide tick and tick label of the big axis
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.ylabel("Read coverage (# of reads)",labelpad=50)

    else:
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axis
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.ylabel("Contig present (yes/no)",labelpad=50)

    if not read_cov_file:
        if not alignment_folder:
            fig.savefig(os.path.join(outdir,'contig_yield_overview.png'),bbox_inches='tight', dpi = 500)
        else:
            fig.savefig(os.path.join(outdir,'contig_alignment_yield_overview.png'),bbox_inches='tight', dpi = 500)
    else:
        fig.savefig(os.path.join(outdir,'contig_alignment_read_cov_yield_overview.png'),bbox_inches='tight', dpi = 500)

def main(args):
    target_contig_folder = args.extracted_contigs
    extracted_contigs = os.path.join(target_contig_folder,'match_table.txt')
    outdir = args.output
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print('Plotting overview of extracted sequences. This may take several minutes...')
    if not args.alignments:
        plot_contig_yield(extracted_contigs,outdir)
    else:
        if not args.read_cov:
            plot_contig_yield(extracted_contigs,outdir,args.alignments)
        else:
            read_cov_info = os.path.join(args.read_cov,'average_cov_per_locus.txt') 
            plot_contig_yield(extracted_contigs,outdir,args.alignments,read_cov_info,norm_value=args.coverage_norm)



