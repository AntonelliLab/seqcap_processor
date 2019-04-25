#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:46:57 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import os
import shutil
import numpy as np
import pandas as pd
from numpy import genfromtxt
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
import glob
import sys
#help(MafftCommandline)

def fix_line_wrap(alignment_file):
    file_in = open(alignment_file)
    final = {}    
    for line in file_in:
        line = line.strip()
        if line.startswith(">"):
            id = line
            final[id] = " "
        else:
            final[id] += line    
    file_out = open(alignment_file, "w")    
    for key, val in final.items():
        file_out.write(key)
        file_out.write("\n")
        file_out.write(val)
        file_out.write("\n")

# the following paths need to be provided
contig_folder = '/Users/tobias/GitHub/seqcap_processor/data/processed/contigs'
ref_file = '/Users/tobias/GitHub/seqcap_processor/data/raw/palm_reference_sequences.fasta'
root_dir = '/Users/tobias/GitHub/seqcap_processor/data/processed/target_contigs_paralogs_info'
outdir = '/Users/tobias/Desktop/test'

if os.path.exists(outdir):
    shutil.rmtree(outdir)
os.makedirs(outdir)
ref_index_info = os.path.join(root_dir,'reference_fasta_header_info.txt')

subdirs = list(os.walk(root_dir))[0][1]
for sample in subdirs:
    print('\nProcessing sample %s'%sample)
    # get the paths for the sample
    sample_path = os.path.join(root_dir,sample)
    contig_orientation = os.path.join(sample_path,'contig_orientation.txt')
    para_info = os.path.join(sample_path,'info_paralogous_loci.txt')
    contig_file = os.path.join(contig_folder,'%s.fa'%sample)
    # read the data
    ref_index_df = pd.read_csv(ref_index_info,sep='\t',header=-1)
    keys = ref_index_df[0].values
    values = ref_index_df[1].values
    id_ref_dict = dict(zip(keys,values))
    ref_seqs = list(SeqIO.parse(ref_file, "fasta"))
    contig_seqs = list(SeqIO.parse(contig_file, "fasta"))
    contig_orientation_df = pd.read_csv(contig_orientation,sep='\t')
    para_data = genfromtxt(para_info, delimiter='\t')
    print('%i paralogous loci found.'%len(para_data))
    sample_out_dir = os.path.join(outdir,sample)
    if not os.path.exists(sample_out_dir):
        os.makedirs(sample_out_dir)
    sample_sequence_outdir = os.path.join(sample_out_dir,'paralog_seq_collections')
    if not os.path.exists(sample_sequence_outdir):
        os.makedirs(sample_sequence_outdir)
    # print ref and contig sequences for each paralogous locus into a separate sequence collection
    for counter,i in enumerate(para_data):
        records = []
        reference = int(i[0])
        reference_id = id_ref_dict[reference]
        ref_seq = [ref for ref in ref_seqs if reference_id==ref.id][0]
        records.append(ref_seq)
        contig_list_tmp = i[1:]
        contig_list = np.unique(contig_list_tmp[~np.isnan(contig_list_tmp)].astype(int).astype(str))
        for contig_id in contig_list:
            contig_seq = [contig for contig in contig_seqs if contig_id == contig.id][0]
            orientation = contig_orientation_df[contig_orientation_df.contig_id == int(contig_id)].orientation.values[0]
            if orientation == '+':
                pass
            else:
                contig_seq.seq = contig_seq.seq.reverse_complement()
            records.append(contig_seq)
        sys.stdout.write('\rPrinting sequence collections %i/%i '%(int(counter+1),len(para_data)))
        SeqIO.write(records, os.path.join(sample_sequence_outdir,'paralog_contigs_collection_locus_%i.fasta'%reference),'fasta')
    print('\rPrinting sequence collections %i/%i '%(int(counter+1),len(para_data)))
    # align the sequences and print as fasta alignment file
    sample_alignment_outdir = os.path.join(sample_out_dir,'paralog_alignments')
    if not os.path.exists(sample_alignment_outdir):
        os.makedirs(sample_alignment_outdir)
    seq_colls = glob.glob(os.path.join(sample_sequence_outdir,'*'))
    for counter, sequence_collection in enumerate(seq_colls):
        filename = sequence_collection.split('/')[-1].replace('paralog_contigs_collection_','alignment_')
        cline = MafftCommandline(input=sequence_collection,op=6.,ep=1.)
        stdout, stderr = cline()
        alignment_out = os.path.join(sample_alignment_outdir,filename)
        sys.stdout.write('\rAligning sequence collections %i/%i '%(int(counter+1),len(para_data)))
        with open(alignment_out, "w") as handle:
            handle.write(stdout)
        fix_line_wrap(alignment_out)




