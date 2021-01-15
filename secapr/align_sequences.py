#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 22:04:36 2021

@author: Tobias Andermann (tobiasandermann88@gmail.com)
"""

import os, glob, sys, shutil
import subprocess
import numpy as np
from secapr.helpers import FullPaths, CreateDir, is_file
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline, MuscleCommandline
import multiprocessing

def add_arguments(parser):
    parser.add_argument(
        "--sequences",
        required=True,
        action=FullPaths,
        type=is_file,
        help='The fasta file containing individual sequences for several samples and loci'
    )
    parser.add_argument(
        "--outdir",
        required=True,
        action=CreateDir,
        help='The directory in which to store the resulting alignments.'
    )
    parser.add_argument(
        "--aligner",
        choices=["muscle", "mafft"],
        default="mafft",
        help='The alignment engine to use.'
    )
    parser.add_argument(
        "--exclude_ambiguous",
        action="store_true",
        default=False,
        help='Don\'t allow reads in alignments containing N-bases.'
    )
    parser.add_argument(
        "--gap_opening_penalty",
        type=float,
        default=3.0,
        help='Set gap opening penalty for aligner.'
    )
    parser.add_argument(
        "--gap_extension_penalty",
        type=float,
        default=2.0,
        help='Set gap extension penalty for aligner.'
    )
    parser.add_argument(
        "--min_seqs_per_locus",
        type=int,
        default=3,
        help='Minimum number of sequences required for building alignment.'
    )
    parser.add_argument(
        "--no_trim",
        action="store_true",
        default=False,
        help='Suppress trimming of alignments. By default secapr uses trimal to trim gappy positions from alignments.'
    )
    parser.add_argument(
        "--trimal_setting",
        choices=["manual", "gappyout", "strict", "strictplus"],
        default="gappyout",
        help='Use one of trimal automated scenarios. These will overwrite all other trimming flags (see below). See trimal tutorial for more info about settings.'
    )
    parser.add_argument(
        "--window_size",
        type=int,
        default=5,
        help='Sliding window size for trimming.'
    )
    parser.add_argument(
        "--seq_proportion",
        type=float,
        default=0.7,
        help='The proportion of sequences required. All alignment columns with fewer sequences will be deleted (0-1).'
    )
    parser.add_argument(
        "--conserve_alignment_percentage",
        type=int,
        default=0,
        help='This setting will ensure to conserve the specified percentage of the alignment when trimming (0-100).'
    )
    parser.add_argument(
        "--min_length",
        type=int,
        default=0,
        help='The minimum length of alignments to keep.'
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help='Process alignments in parallel using --cores for alignment. This is the number of PHYSICAL CPUs.'
    )

def align_seqs(pool_input):
    counter,total,sequence_collection,aligner,gap_opening_penalty,gap_extension_penalty,no_trim,trimal_setting,window_size,seq_proportion,conserve_alignment_percentage,min_length,outdir = pool_input
    filename = os.path.basename(sequence_collection).replace('sequence_collection_locus_','')
    if aligner == 'mafft':
        cline = MafftCommandline(input=sequence_collection,adjustdirection=True,maxiterate=1000,op=gap_opening_penalty,ep=gap_extension_penalty)
    elif aligner == 'muscle':
        cline = MuscleCommandline(input=sequence_collection,maxiters=1000,gapopen=gap_opening_penalty,gapextend=gap_extension_penalty)
    stdout, stderr = cline()
    alignment_out = os.path.join(outdir,filename)
    sys.stdout.write('\rAligning sequence collections %i/%i '%(int(counter+1),total))
    sys.stdout.flush()
    with open(alignment_out, "w") as handle:
        handle.write(stdout)
    
    if not no_trim:
        # trim alignments with trimal
        if trimal_setting != 'manual':
            cmd = [
                "trimal",
                "-in",
                alignment_out,
                "-out",
                alignment_out,
                '-%s'%trimal_setting
                   ]
        else:
            cmd = [
                "trimal",
                "-in",
                alignment_out,
                "-out",
                alignment_out,
                '-w',
                str(window_size),
                '-gt',
                str(seq_proportion),
                '-cons',
                str(conserve_alignment_percentage)
                   ]                
        # run trimal command
        proc = subprocess.Popen(cmd,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE
            )
        stderr,stdout = proc.communicate()
    if min_length:
        align = AlignIO.read(alignment_out, "fasta")
        al_length = len(align[0])
        if al_length < min_length:
            # delete file if smaller than minlength
            os.remove(alignment_out)
            #too_short_alignments.append(filename.replace('.fasta',''))
            return(filename.replace('.fasta','')) # Return locus name in case alignment is too short


def main(args):
    sequences = args.sequences
    outdir = args.outdir
    aligner = args.aligner
    exclude_ambiguous = args.exclude_ambiguous
    gap_opening_penalty = args.gap_opening_penalty
    gap_extension_penalty = args.gap_extension_penalty
    min_seqs_per_locus = args.min_seqs_per_locus    
    no_trim = args.no_trim
    trimal_setting = args.trimal_setting
    window_size = args.window_size
    seq_proportion = args.seq_proportion
    conserve_alignment_percentage = args.conserve_alignment_percentage
    min_length = args.min_length
    cores = args.cores
    

    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # create locus dictionary
    locus_dict = {}
    for record in SeqIO.parse(sequences, "fasta"):
        locus = record.description.split("|")[1]
        locus_dict.setdefault(locus,[])
        if 'N' in record.seq and exclude_ambiguous:
            print('Found ambiguous bases in sequence %s. Dropping this sequence. Disable --exclude_ambiguous if you want to keeo these sequences.'%record.description)
            pass
        else:
            locus_dict[locus].append(record)
    
    # check which loci have to be dropped due to too few loci
    locus_names = list(locus_dict.keys())
    sequences_per_locus = np.array([len(locus_dict[i]) for i in locus_names])
    exclude_these_loci = list(np.array(locus_names)[sequences_per_locus<min_seqs_per_locus])
    
    # write sequence collections to file
    tmp_dir = os.path.join(outdir,'tmp')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    for locus in locus_names:
        if locus not in exclude_these_loci:            
            SeqIO.write(locus_dict[locus], os.path.join(tmp_dir,'sequence_collection_locus_%s.fasta'%locus),'fasta')
    
    # align sequence collections
    #too_short_alignments = []
    seq_collections = glob.glob(tmp_dir+'/*.fasta')[:40] #XXXXXXXXXXXXXXXXXXX REMOVE [:20]
    
    
    pool_input = [[i,len(seq_collections),seq_collection,aligner,gap_opening_penalty,gap_extension_penalty,no_trim,trimal_setting,window_size,seq_proportion,conserve_alignment_percentage,min_length,outdir] for i,seq_collection in enumerate(seq_collections)]
    pool = multiprocessing.Pool(cores)
    all_output = np.array(pool.map(align_seqs, pool_input))
    pool.close()
               
    too_short_alignments = all_output[~(all_output == None)].astype(str) # get the loci if there are any
    
    shutil.rmtree(tmp_dir)
    print('\n\nExcluded these loci due to too few taxa present (< %i):\n'%min_seqs_per_locus, sorted(exclude_these_loci))
    if len(too_short_alignments) > 0:
        print('\n\nExcluded these loci due to too short alignments (< %i):\n'%min_length,sorted(too_short_alignments))
    
    












