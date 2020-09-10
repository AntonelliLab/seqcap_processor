#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:16:59 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



contig_headers = ['109189', '156474', '490477', '183356', '244827', '249589', '376762', '473401', '506353', '552484']
contig_headers = ['549799', '593734']
contig_headers = ['252306', '271392', '399821', '180012', '383994', '497380', '573398']
lastz_data = pd.read_csv('/Users/tobias/Desktop/target_contigs_test/1061/1061.lastz',sep='\t')
contig_header_values = np.array([i.split(' ')[0].replace('>','') for i in lastz_data.name1.values if i.split(' ')[0].replace('>','') in contig_headers]).astype(int)
contig_length_values = np.array([i.split(' ')[1] for i in lastz_data.name1.values if i.split(' ')[0].replace('>','') in contig_headers]).astype(int)
longest_contig = contig_header_values[list(contig_length_values).index(np.max(contig_length_values))]
