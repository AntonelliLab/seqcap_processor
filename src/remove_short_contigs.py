#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 15:21:03 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



contig_file = '/Users/tobias/Desktop/1063.fa'
#contig_file_new = '/Users/tobias/GitHub/seqcap_processor/data/processed/contigs/1063_removed_short_contigs.fa'
min_length = 200
fasta =  open(contig_file,'r')
fasta_content = list(fasta)
counter = 0
indeces_to_keep = []
for i,line in enumerate(fasta_content):
    if not line.startswith('>'):
        contig_length = len(line.replace('\n',''))
        if contig_length < min_length:
            pass
        else:
            # line number of header
            indeces_to_keep.append(i-1)
            # line number of sequence
            indeces_to_keep.append(i)

new_fasta_content = list(np.array(fasta_content)[indeces_to_keep])
new_fasta = open(contig_file,'w')
for line in new_fasta_content:
    new_fasta.write(line)
new_fasta.close()


