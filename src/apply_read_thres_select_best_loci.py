#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 10:27:32 2017

@author: tobias
"""

import pandas as pd

# Define a threshold of read-depth (anything below will be considered as a locus not present in the data)
threshold = 3

# Load the input file with all read-coverage values
cov_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/selected_loci/average_cov_per_locus.txt'
coverage_df = pd.read_csv(cov_file,sep='\t')

# Return boolean for every field, depending on if its greater than the threshold
thres_test = coverage_df.ix[:,1:]>threshold

# Extract only those rows for which all fields returned 'True' and store in new df
selected_rows = pd.DataFrame([])
for line in thres_test.iterrows():
    line = line[1]
    if line.all():
        selected_rows = selected_rows.append(line)

# Store all indices of the selected data (selected_rows) in a list
indeces = list(selected_rows.index.get_values())

# Use indices to extract rows from oriignal df and create new one from it
loci_passing_test = coverage_df.iloc[indeces,:].copy()
list_of_good_loci = list(loci_passing_test.locus)

# Calculate the read-depth sum across all samples for each locus and store as new column in df 
loci_passing_test['sum'] = loci_passing_test.ix[:,1:].sum(axis=1)

# Sort the df by the 'sum' column to have the best covered loci on top
loci_passing_test.sort_values('sum', axis=0, ascending=False, inplace=True)

# Write the sorted df to a csv file
loci_passing_test.to_csv('/Users/tobias/GitHub/seqcap_processor/data/processed/selected_loci/loci_passing_threshold_%s' %str(threshold), sep = '\t', index = False)

