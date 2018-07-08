#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 10:47:27 2018

@author: tobias
"""

import numpy as np
import pandas as pd
from functools import reduce

read_cov_overview = '/Users/tobias/GitHub/seqcap_processor/data/processed/remapped_reads/average_cov_per_locus.txt'
read_cov_overview_df = pd.read_csv(read_cov_overview,sep='\t')

loci_names = read_cov_overview_df['locus'].values
good_exons = []
sample_exon_count_dict = {}
for sample in read_cov_overview_df.columns:
    if not sample =='locus':
        values = read_cov_overview_df[sample].values
        num_loci_high_coverage = len(values[values>3])
        sample_exon_count_dict.setdefault(sample,num_loci_high_coverage)
        good_loci_names = loci_names[values>3]
        good_exons.append(good_loci_names)

loci_present_in_all_samples = reduce(np.intersect1d, (good_exons))


exons_per_sample = list(sample_exon_count_dict.values())
np.mean(exons_per_sample)
np.std(exons_per_sample)
