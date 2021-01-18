#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 19:01:42 2017

@author: tobias
"""

# this program estimates some basic population genetic parameters from allele alignments

import os
import itertools
import numpy as np
import random
import matplotlib.pyplot as plt
from Bio import AlignIO


input_dir = '/Users/tobias/Desktop/abc_modeling_course_tjarno_2017/my_data/topaza-uce-allele-alignments'
aln_format = 'fasta'
min_num_seqs = 20
output_dir = '/Users/tobias/Desktop/abc_modeling_course_tjarno_2017/my_data/'
outgroup = 'Florisuga'


def read_aln_in_folder(input_dir, aln_format):
    file_list = []
    for file in os.listdir(input_dir):
        if file.endswith(".%s" %aln_format):
            file_list.append(os.path.join(input_dir, file))
    return file_list



def get_variable_positions(aln_list,aln_format, min_num_seqs):

    # this dict will only contain the valid columns withut ambiguities
    valid_columns_dict = {}
    # this dict will contains all variable positions
    variable_columns_dict = {}  

    for file in aln_list:
        alignment = AlignIO.read(open(file), aln_format)
        if len(alignment) < min_num_seqs:
            pass
        else:
            # get the name of the alignment
            name = file.split('/')[-1]
            # iterate through columns of each alignment
            for r in range(0,alignment.get_alignment_length()):
                not_real = False
                column = alignment[:,r]
                states = set(column)
                if not 'N' in states:
                    valid_columns_dict.setdefault(name,[])
                    valid_columns_dict[name].append(r)
                    valid_states = []
                    gap_openings = []
                    for i in range(0,len(states)):
                        state = list(states)[i]
                        # only these states count as true bases
                        if state in ['A','C','T','G']:
                            valid_states.append(state)
                        elif state == '-':
                            valid_states.append(state)
                            if '-' in alignment[:,r-1]:
                                not_real = True
                        
                    if len(valid_states) > 1:
                        variable_columns_dict.setdefault(name,[])
                        if not_real:
                            pass
                        else:
                            variable_columns_dict[name].append(r)
                    
    return valid_columns_dict,variable_columns_dict
                

def get_possible_combinations_for_n_sequences(n):
    combo_list = []
    for subset in itertools.combinations(range(0,n), 2):
        if subset not in combo_list:
            combo_list.append(subset)
    return combo_list


def tajimas_estimator_per_locus(input_dir,aln_format,snp_columns,valid_columns,min_num_seqs):
    locus_tajima_dict = {}
    locus_snp_count_dict = {}
    for locus in snp_columns: #iterating over all loci
        alignment = AlignIO.read(open(os.path.join(input_dir,locus)), aln_format)
        d_count = 0    
        # iterate through columns
        snp_pos_freq_dict = {}
        for position in snp_columns[locus]:
            column = alignment[:,position]
            snp_pos_freq_dict.setdefault(position,[[x,column.count(x)] for x in set(column)])
            if len([[x,column.count(x)] for x in set(column)]) <= 2:
                # p = occurrence of one of the two alleles, doesn't matter which one, since it's always simmetrical
                p = [[x,column.count(x)] for x in set(column)][0][1]
                n =  min_num_seqs
                # this formula calculates the pairwise differences between all possible pairs (formula derived by myself by seeing what fits the actually counted pairwise differnces for different n's and p's)
                delta = n*p-p**2
                d_count += delta
            # deal with positions with more than 2 variants by counting 'manually'
            else:
                for pair in possible_combos:
                    if alignment[pair[0],position] != alignment[pair[1],position]:
                        d_count += 1
        locus_snp_count_dict.setdefault(locus,[])
        locus_snp_count_dict[locus].append(snp_pos_freq_dict)
        tajimas_estimator = d_count/len(possible_combos)
        locus_len = len(valid_columns[locus])
        corrected_tajimas_estimator = tajimas_estimator/locus_len
        locus_tajima_dict.setdefault(locus,corrected_tajimas_estimator)
    return locus_tajima_dict,locus_snp_count_dict


def plot_expected_heterozyosity(locus_tajima_dict,output_dir):
    theta_estimates = np.array(list(locus_tajima_dict.values()))
    def expected_heterozygosity(theta):
        y = theta/(theta+1)
        return y
    # create a series of numbers in order to plot the function
    t = np.arange(0.0001, 1000.0, 0.001)
    f=plt.figure()
    plt.plot(theta_estimates, expected_heterozygosity(theta_estimates), 'ro',label="calculated tajimas-theta-estimator",alpha=0.2)
    plt.axis([0.00009, 1005.0, -0.1, 1.1])
    plt.xscale('log')
    plt.plot(t, expected_heterozygosity(t), label="expected heterozygosity function")
    plt.legend(loc='upper left')
    plt.xlabel('theta')
    plt.ylabel('expected heterozygosity')
    plt.title('Expected heterozygosity for UCE loci of Topaza (n = %s)'%len(theta_estimates))
    plt.show()
    f.savefig(os.path.join(output_dir,"expected_heterozygosity_uces.pdf"), bbox_inches='tight')


def get_fsf_stats(input_dir,aln_format,outgroup,locus_snp_count_dict):
    # this function counts how often every new mutation is shared among all sampled alleles
    locus_allele_counts_dict = {}
    for locus in locus_snp_count_dict:
        locus_allele_counts_dict.setdefault(locus,[])
        alignment = AlignIO.read(open(os.path.join(input_dir,locus)), aln_format)
        for position in locus_snp_count_dict[locus]:
            position_list = position.keys()
            # iterate through the target columns of the alignment (containing SNPs)
            for pos in position_list:
                ancestral = ''
                outgroup_pos = []
                ingroup_pos = []
                for sequence in alignment:
                    # make sure to catch the outgroup state to estimate the ancestral state
                    if outgroup in sequence.id:
                        outgroup_pos.append(sequence.seq[int(pos)])
                    else:
                        ingroup_pos.append(sequence.seq[int(pos)])
                if len(set(outgroup_pos))>1:
                    secure_random = random.SystemRandom()
                    ancestral_proposal = secure_random.choice(outgroup_pos) 
                else:
                    ancestral_proposal = outgroup_pos[0]             
                if ancestral_proposal in ingroup_pos:
                    ancestral = ancestral_proposal
                else:
                    count = [[x,ingroup_pos.count(x)] for x in set(ingroup_pos)]
                    # sort the count list by the number of occurrences
                    count.sort(key=lambda x: x[1])
                    # get the most frequent occurrence
                    ancestral = count[-1][0]
                base_count_dict = {}
                [base_count_dict.setdefault(x,ingroup_pos.count(x)) for x in set(ingroup_pos)]
                for base in ingroup_pos:
                    if base != ancestral:
                        locus_allele_counts_dict[locus].append(base_count_dict[base])
    return locus_allele_counts_dict


alignment_list = read_aln_in_folder(input_dir,aln_format)
valid_columns, snp_columns = get_variable_positions(alignment_list,aln_format,min_num_seqs)
possible_combos = get_possible_combinations_for_n_sequences(min_num_seqs)
locus_tajima_dict,locus_snp_count_dict = tajimas_estimator_per_locus(input_dir,aln_format,snp_columns,valid_columns,min_num_seqs)
plot_expected_heterozyosity(locus_tajima_dict,output_dir)
locus_allele_counts_dict = get_fsf_stats(input_dir,aln_format,outgroup,locus_snp_count_dict)

                
            
            
# REMOVE OUTGROUP!!!!!!!!!                
                
     
# get the site frequency spectrum stats:
    # determine the ancestral state at each polymorphism:
        # 1. if outgroup allele also present in ingroup, take it as ancestral state
        # 2. if outgroup allele not present in ingroup, choose the most common ingroup allele as ancestral
        # 3. if all allele copies are equally common, choose one at random as ancestral
    # for all derived mutations, count by how many sequences those are shared




# find good way to partition data by population

# what to do with columns containing N's?? (we are excluding them for now)



