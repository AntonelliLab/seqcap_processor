#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:33:17 2018

@author: tobias
"""


df_1 = '/Users/tobias/GitHub/seqcap_processor/data/processed/cleaned_trimmed_reads_test/sample_overview.txt'
df_a = pd.read_csv(df_1,sep='\t')
        
        
match_table = '/Users/tobias/GitHub/seqcap_processor/data/processed/target_contigs/match_table.txt'
table = pd.read_csv(table,sep='\t',index_col=0)      


sample_count_dict = {}
for column in table.columns:
	sample_count_dict.setdefault(column.replace('sample_',''),sum(table[column]))

type(new_value)



test_str = '/Users/tobias/GitHub/seqcap_processor/data/processed/contigs_test/stats/1063.fa'
'/'.join(test_str.split('/')[:-2])

import os
import subprocess
file = contig_file
lines = 2
    tail_out = subprocess.Popen(['tail','-n', lines,file],stdout=subprocess.PIPE)
    tail_out.communicate[0]
    lines = stdout.readlines(); stdout.close()
    return lines[:,-offset]

tail(contig_file,3)

contig_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/contigs/sample_1061.fa'


stats_file = '/Users/tobias/GitHub/seqcap_processor/data/processed/cleaned_trimmed_reads/1085_clean/1085_stats.txt'
def get_read_count_from_stats_file(stats_file):
    F = open(stats_file,'r') 
    for line in F:
        if line.startswith('Input'):
            reads_before = line.split(' ')[3]
            reads_after = line.split(' ')[6]
    return(reads_before,reads_after)
