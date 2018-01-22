#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:33:17 2018

@author: tobias
"""



df_1 = '/Users/tobias/GitHub/seqcap_processor/data/processed/cleaned_trimmed_reads/sample_overview.txt'
df_2 = '/Users/tobias/GitHub/seqcap_processor/data/processed/contigs_test/sample_stats.txt'
df_a = pd.read_csv(df_1,sep='\t')
df_b = pd.read_csv(df_2,sep='\t')

new = df_a.iloc[[7,9],:].copy().reset_index()
new['contig_count'] = df_b.total_contig_count
counter = 0
for index,row in df_a.iterrows():
    sample_name = row['sample']
    if sample_name in list(df_b['sample']):
        print(row)
        new_info = df_b[df_b['sample']==sample_name]['total_contig_count']
        new_value = new_info.values[0]
        new_name = new_info.name
        headers = np.array(row.index)
        old_values = row.values
        new_index = np.append(headers,new_name)
        new_values = np.append(old_values,new_value)
        if counter == 0:
            new_values_previous = new_values
        else:
            new_values_previous = np.stack((new_values_previous, new_values), axis=0)
        counter += 1


pd.DataFrame(data=new_values,index=False,columns=new_index)
        
        
        new_values.
        
        
        
        new_item = pd.Series({new_name:new_value})
        new_row = row.append(new_item)


      row.index
        if counter ==0:
            previous_row = new_row
        else:
            previous_row = pd.concat([previous_row, new_row], axis=1)
        counter +=1
pd.DataFrame(previous_row)
previous_row.values
previous_row.index

pd.DataFrame(data=data[1:,1:],index=data[1:,0],columns=data[0,1:])


        print(new_row)
        
        
        
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
