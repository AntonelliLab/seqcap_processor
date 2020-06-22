#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:08:28 2020

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import sys,os,glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)
np.random.seed(1234)

from Bio import SeqIO
file = '/Users/tobias/GitHub/seqcap_processor/data/test/target_sequences/extracted_target_contigs_all_samples.fasta'
contigs_file_content = SeqIO.parse(open(file),'fasta')
seqlengths = []
for i in contigs_file_content:
    seqlength = len(str(i.seq))
    seqlengths.append(seqlength)


np.mean(np.array(seqlengths))
