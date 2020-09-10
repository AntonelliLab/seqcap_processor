#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 13:48:44 2017

@author: tobias
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Get the data as pandas dataframe
log_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/selected_loci/average_cov_per_locus.txt'
data_input = pd.read_csv(log_file, sep = '\t')
