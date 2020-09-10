#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:46:32 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

old_contigs = '/Users/tobias/Desktop/extracted_target_contigs_all_samples_old.fasta'
new_contigs = '/Users/tobias/Desktop/extracted_target_contigs_all_samples_new.fasta'

from Bio import AlignIO
from Bio import SeqIO

contig_lengths_old = []
for record in SeqIO.parse(old_contigs, "fasta"):
    contig_lengths_old.append(len(record.seq))
contig_lengths_new = []
for record in SeqIO.parse(new_contigs, "fasta"):
    contig_lengths_new.append(len(record.seq))

len(contig_lengths_old)
len(contig_lengths_new)

np.mean(contig_lengths_old)
np.mean(contig_lengths_new)

plt.hist(contig_lengths_old,100)
plt.title('Old contigs')

plt.hist(contig_lengths_new,100)
plt.title('New contigs')