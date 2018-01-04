#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 16:51:25 2017

@author: tobias
"""
import os
import re
import random
import subprocess
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna

def read_fasta(fasta):
    name, seq = None, []
    for line in fasta:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def get_seq_dict(dictionary, fasta_path):
    new_dict = {}
    for locus in dictionary.keys():
        new_dict.setdefault(locus,[])
    alignment = SeqIO.parse(open(fasta_path),'fasta')
    for seq in alignment:
        ref_loc = str(seq.id).split('_')[0]
        new_dict[ref_loc].append(seq)
        seq.name=''
        seq.description=''
    return new_dict


def create_reference_fasta(out_dir,path_to_alignments):
    # Create a list of fasta files from the input directory
    file_list = [fn for fn in os.listdir(path_to_alignments) if fn.endswith(".fasta")]
    reference_list = []
    for fasta_alignment in file_list:
        sequence_name = re.sub(".fasta","",fasta_alignment)
        orig_aln = os.path.join(path_to_alignments,fasta_alignment)
        sep_reference = "%s/%s" %(out_dir,fasta_alignment)
        reference_list.append(sep_reference)
        cons_cmd = "cons -sequence %s -outseq %s -name %s -plurality 0.1 -setcase 0.1" %(orig_aln,sep_reference,sequence_name)
        os.system(cons_cmd)
    reference = os.path.join(out_dir,"joined_fasta_library.fasta")
    join_fastas = "cat %s/*.fasta > %s" %(out_dir,reference)
    os.system(join_fastas)
    return reference



fasta_file = '/Users/tobias/Desktop/cos2.fasta'

locus_bait_dict = {}
with open(fasta_file) as f:
    for name, seq in read_fasta(f):
        locus_name = re.sub('>','',name.split('_')[0])
        locus_bait_dict.setdefault(locus_name,[])
        locus_bait_dict[locus_name].append(seq)

locus_fasta_dict = get_seq_dict(locus_bait_dict,fasta_file)

            

out_path = '/Users/tobias/Desktop/merging_probes/sequence_files'
for locus in locus_fasta_dict:
    filename = '%s_sequences.fasta' %locus
    with open(os.path.join(out_path,filename), "w") as out_file:
        seq_list = locus_fasta_dict[locus]
        index = 0
        for sequence in seq_list:
            sequence.id = '%s_%i' %(locus,index)
            sequence.name=''
            sequence.description=''                
            index += 1
            out_file.write(sequence.format('fasta'))

# align the sequence fasta files
aln_path = '/Users/tobias/Desktop/merging_probes/alignments'
for fasta in os.listdir(out_path):
    fasta_file = os.path.join(out_path,fasta)
    new_file_name = re.sub('_sequences.fasta','_sequence_alignment.fasta',fasta)
    aln = os.path.join(aln_path,new_file_name)
    aln_stdout = open(aln, 'w')
    # run MAFFT on the temp file
    cmd = ["mafft","--maxiterate", "1000", fasta_file]
    # just pass all ENV params
    proc = subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=aln_stdout)
    stderr = proc.communicate()
    aln_stdout.close()


output = '/Users/tobias/Desktop/merging_probes/new_reference'
create_reference_fasta(output,aln_path)

