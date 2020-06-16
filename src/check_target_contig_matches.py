#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 20:32:09 2020

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import sys,os,glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

np.set_printoptions(suppress=True)
np.random.seed(1234)


def contigs_matching_exons(lastz_df):
    # make a dictionary with all contig names that match a exon locus
    exon_contig_dict = {}
    contig_exon_dict = {}
    contig_orientation_dict = {}
    contig_multi_exon_dict = {}
    for row in lastz_df.iterrows():
        locus = row[1].name2
        locus_name = re.sub('\_p[0-9]* \|.*', '', locus)
        locus_name = re.sub('^>', '', locus_name)
        contig_header = row[1].name1
        #print(contig_header)
        #contig_name = re.sub('^\>([0-9]*) .*', '\\1', contig_header)
        contig_name = re.sub('^\>([^\s]*) .*', '\\1', contig_header)
        contig_name = re.sub('^>', '', contig_name)
        #print(contig_name)
        exon_contig_dict.setdefault(locus_name,[])
        exon_contig_dict[locus_name].append(contig_name)
        contig_exon_dict.setdefault(contig_name,[])
        contig_exon_dict[contig_name].append(locus_name)
        orientation = row[1].strand2
        contig_orientation_dict.setdefault(contig_name,orientation)
    
    orientation_df = pd.DataFrame(index=np.arange(0,len(contig_orientation_dict)), columns=['contig_id','orientation'])
    orientation_df['contig_id'] = list(contig_orientation_dict.keys())
    orientation_df['orientation'] = list(contig_orientation_dict.values())
    for contig in list(contig_exon_dict.keys()):
        if len(contig_exon_dict[contig]) > 1:
            contig_multi_exon_dict.setdefault(contig,contig_exon_dict[contig])
    return exon_contig_dict, contig_exon_dict, contig_orientation_dict, contig_multi_exon_dict, orientation_df


def find_duplicates(exon_contig_dict,contig_exon_dict):
    # get exons that have multiple contigs matching them = paralogs
    invalid_exon_loci = []
    exons_with_multiple_hits = []
    for exon in list(exon_contig_dict.keys()):
        if len(exon_contig_dict[exon]) > 1:
            exons_with_multiple_hits.append(exon)
            invalid_exon_loci.append(exon)
    # get contigs matching multiple exons = long contigs
    contigs_matching_multiple_exons = []
    for contig in list(contig_exon_dict.keys()):
        if len(contig_exon_dict[contig]) > 1:
            contigs_matching_multiple_exons.append(contig)
            for exon in contig_exon_dict[contig]:
                invalid_exon_loci.append(exon)
    return sorted(list(np.unique(np.array((invalid_exon_loci))))), exons_with_multiple_hits, contigs_matching_multiple_exons

def find_longest_contig(contig_name,lastz_df):
    contig_lengths = np.array([int(i.split('length_')[1].split('_')[0]) for i in contig_name])
    longest_contig = contig_name[np.where(np.max(contig_lengths))[0][0]]
    return longest_contig


def get_list_of_valid_exons_and_contigs(exon_contig_dict,loci_with_issues,possible_paralogous,contigs_matching_multiple_exon_dict,keep_duplicates_boolean,keep_paralogs_boolean,outdir,lastz_df):
    # summarize all exons that should be excluded form further processing (duplicates)

    # remove all duplicates
    invalid_exons_unique = list(set(loci_with_issues))
    
    if keep_duplicates_boolean:
        # keep all exons that are only affected by possible duplicate issues (i.e. long exons)
        exons_affected_by_dupls = [contigs_matching_multiple_exon_dict[i] for i in contigs_matching_multiple_exon_dict.keys()]
        exons_affected_by_dupls = sorted(list(set([item for sublist in exons_affected_by_dupls for item in sublist])))
        keep_these_exons = list(set(exons_affected_by_dupls)-set(possible_paralogous))
        invalid_exons_unique = list(set(invalid_exons_unique)-set(keep_these_exons))
        
    if keep_paralogs_boolean:
        # keep the potential paralogs
        keep_these_exons = possible_paralogous
        invalid_exons_unique = list(set(invalid_exons_unique)-set(keep_these_exons))
        print('Warning: Found %i paralogous loci. The longest matching contig for each paralogous locus will be kept, due to the use of the --keep_paralogs flag. It is not recommendable to use paralogous loci for phylogenetic inference!'%len(invalid_exons_temp))

    print('Removing',len(invalid_exons_unique),'exon loci with potential issues.')

    # print the info to file
    # paralog info file
    paralogous_exons = {}
    for exon in possible_paralogous:
        paralogous_exons.setdefault(exon,exon_contig_dict[exon])
    paralog_info = pd.DataFrame.from_dict(paralogous_exons, orient='index')
    paralog_info.to_csv(os.path.join(outdir,'info_paralogous_loci.txt'),header=False,sep="\t")
    # dupl info file
    duplicate_info = pd.DataFrame.from_dict(contigs_matching_multiple_exon_dict, orient='index')
    duplicate_info.to_csv(os.path.join(outdir,'info_exons_spanning_multiple_loci.txt'),header=False,sep="\t")
    
    # get list of valid contig names
    valid_contig_names = []
    for exon in exon_contig_dict:
        if exon not in invalid_exons_unique:
            contig_name = exon_contig_dict[exon]
            contig_name = find_longest_contig(contig_name,lastz_df)
            valid_contig_names.append(str(contig_name).replace('>',''))
    return valid_contig_names


def extract_target_contigs(sample_id,contig_sequences,valid_contig_names,contig_exon_dict,contig_orientation_dict,subfolder):
    printed_contigs_counter = 0
    # define the output file where extracted contigs will be stored
    global_match_output_name = 'extracted_target_contigs_all_samples.fasta'
    global_match_output_file = os.path.join('/'.join(subfolder.split('/')[:-1]),global_match_output_name)
    sample_match_output_name = 'extracted_target_contigs_%s.fasta'%sample_id
    sample_match_output_file = os.path.join(subfolder,sample_match_output_name)
    # extract valid contigs form contig file and print to fasta file with exon-names+ sample_id as headers
    with open(global_match_output_file, "a") as out_file:
        with open(sample_match_output_file, "w") as sample_file:
            for fasta in contig_sequences:
                if str(fasta.id) in valid_contig_names:
                    orientation = contig_orientation_dict[fasta.id]
                    if orientation == '-':
                        seq = fasta.seq.reverse_complement()
                    else:
                        seq = fasta.seq
                    # get the corresponding exon locus name from the dictionary
                    if len(contig_exon_dict[fasta.id])>1:
                        for matching_locus in contig_exon_dict[fasta.id]:
                            header = '%s_%s |%s' %(matching_locus,sample_id,matching_locus)
                            new_fasta = SeqRecord(seq, id=header, name='', description='')
                            out_file.write(new_fasta.format('fasta'))
                            sample_file.write(new_fasta.format('fasta'))
                            printed_contigs_counter += 1                         
                    else:        
                        header = '%s_%s |%s' %(contig_exon_dict[fasta.id][0],sample_id,contig_exon_dict[fasta.id][0])
                        new_fasta = SeqRecord(seq, id=header, name='', description='')
                        out_file.write(new_fasta.format('fasta'))
                        sample_file.write(new_fasta.format('fasta'))
                        printed_contigs_counter += 1 
    return printed_contigs_counter



contig_file = '/Users/tobias/GitHub/seqcap_processor/data/test/spades_output/first_try/contigs.fasta'
reference_fasta = '/Users/tobias/GitHub/seqcap_processor_old/data/processed/target_contigs/formatted_reference_library.fasta'
subfolder = '/Users/tobias/GitHub/seqcap_processor/data/test/target_contigs/spades_first_try/'
output_file = os.path.join(subfolder,'lastz_out.lastz')
min_coverage = 80
min_identity = 80
keep_duplicates = True
keep_paralogs = False
sample_id = '1061'

with open(output_file, 'w') as lastz_out_file:
    lastz_command = [
        'lastz',
        '%s[multiple,nameparse=full]'%contig_file,
        '%s[nameparse=full]'%reference_fasta,
        '--strand=both',
        '--seed=12of19',
        '--transition',
        '--nogfextend',
        '--nochain',
        '--gap=400,30',
        '--xdrop=910',
        '--ydrop=8370',
        '--hspthresh=3000',
        '--gappedthresh=3000',
        '--noentropy',
        '--coverage=%i'%min_coverage,
        '--identity=%i'%min_identity,
        '--ambiguous=iupac',
        '--format=general:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity'
    ]
    run_lastz = subprocess.Popen(lastz_command, stdout=lastz_out_file, stderr=None)
    run_lastz.communicate()
    
lastz_df = pd.read_csv(output_file,sep='\t')
# store the data in dictionaries for convenience
exon_contig_dict, contig_exon_dict, contig_orientation_dict, contig_multi_exon_dict, orientation_df = contigs_matching_exons(lastz_df)

# mark duplicate loci
loci_with_issues, possible_paralogous, contigs_covering_several_loci = find_duplicates(exon_contig_dict,contig_exon_dict)
# remove duplicate loci from the list of targeted loci and contigs
target_contigs = get_list_of_valid_exons_and_contigs(exon_contig_dict,loci_with_issues,possible_paralogous,contig_multi_exon_dict,keep_duplicates,keep_paralogs,subfolder,lastz_df)
# load the actual contig sequences
contig_sequences = SeqIO.parse(open(contig_file),'fasta')
# write those contigs that match the reference library to the file
extracted_contig_counter = extract_target_contigs(sample_id,contig_sequences,target_contigs,contig_exon_dict,contig_orientation_dict,subfolder)
print('Extracted %i contigs matching reference exons\n' %extracted_contig_counter)








