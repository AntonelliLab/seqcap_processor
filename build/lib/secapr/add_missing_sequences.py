# encoding: utf-8
#author: Tobias Hofmann, tobias.andermann@bioenv.gu.se
"""
This script will add dummy sequences '?' for missing taxa in each alignments, making sure that all alignments in the input folder contain the same taxa (as required for e.g. *BEAST).",
"""

import os
import sys
import glob
import shutil
import configparser
import pickle
from .utils import CompletePath

def add_arguments(parser):
    parser.add_argument(
        '--input',
        required=True,
        action=CompletePath,
        default=None,
        help='The directory containing fasta alignments'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CompletePath,
        default=None,
        help='The output directory where results will be safed'
    )


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


def main(args):
    # Set working directory
    work_dir = args.input
    out_dir = args.output
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Create a dictionary with the name-pattern as key and all file-names sharing that name-pattern
    fasta_dict = {}
    for fasta in os.listdir(work_dir):
        if fasta.endswith(".fasta") or fasta.endswith(".fa"):
            fasta_dict.setdefault("all",[]).append(fasta)



    # Get the list of taxa names (headers) for each locus, key is out-file, values are in-files
    for key, value in fasta_dict.items():
        # Creates a list of all headers that are present in the concatenated alignments, accounting for differences in the taxon composition of each alignment
        list_headers=[]
        # Each k is a separate fasta input file belonging to the same locus ()to be joined)
        for k in sorted(value):
            with open("%s/%s" %(work_dir,k)) as f:
                for name, seq in read_fasta(f):
                    if not name in list_headers:
                        list_headers.append(name)

        # "value" is a list of all fasta files to be concatenated
        # Find the missing taxa in each fasta input file and simulate a sequence of correct length (only "n")
        for k in sorted(value):
            taxa_names_single = []
            present_seq = []
            length_alignment = 0
            with open("%s/%s" %(work_dir,k)) as f:
                for name, seq in read_fasta(f):
                    taxa_names_single.append(name)
                    present_seq.append((name,seq))
                    length_alignment = len(seq)
            # Make a list of all missing taxa in each fasta input file
            missing_taxa = []
            for header in list_headers:
                if header not in taxa_names_single:
                    missing_taxa.append(header)
            simulated_seq = []
            for mistax in missing_taxa:
                fake_string = "n" * length_alignment
                simulated_seq.append((mistax,fake_string))
            all_seq = sorted(simulated_seq+present_seq)
            out_fasta = open(os.path.join(out_dir, k), 'w')
            for seqname, sequence in all_seq:
                out_fasta.write(seqname+"\n")
                out_fasta.write(sequence+"\n")
            out_fasta.close()
    try:        
        pickle_in = os.path.join(args.input,'.secapr_files/sequence_origin.pickle')
        with open(pickle_in, 'rb') as handle:
            sequence_origin = pickle.load(handle)
        pickle_path = os.path.join(args.output,'.secapr_files/sequence_origin.pickle')
        with open(pickle_path, 'wb') as handle:
            pickle.dump(sequence_origin, handle, protocol=pickle.HIGHEST_PROTOCOL)
    except:
        pass