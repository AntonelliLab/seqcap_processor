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

def contigs_matching_exons(blast_df):
    # make a dictionary with all contig names that match a exon locus
    exon_contig_dict = {}
    contig_exon_dict = {}
    contig_orientation_dict = {}
    contig_multi_exon_dict = {}
    for row in blast_df.iterrows():
        locus = str(row[1].qseqid)
        locus_name = re.sub('\_p[0-9]* \|.*', '', locus)
        locus_name = re.sub('^>', '', locus_name)
        contig_header = str(row[1].sseqid)
        #print(contig_header)
        #contig_name = re.sub('^\>([0-9]*) .*', '\\1', contig_header)
        contig_name = re.sub('^\>([^\s]*) .*', '\\1', contig_header)
        contig_name = re.sub('^>', '', contig_name)
        #print(contig_name)
        exon_contig_dict.setdefault(locus_name,[])
        exon_contig_dict[locus_name].append(contig_name)
        contig_exon_dict.setdefault(contig_name,[])
        contig_exon_dict[contig_name].append(locus_name)
        orientation = row[1].sstrand
        contig_orientation_dict.setdefault(contig_name,orientation)
    
    # remove double listings of loci/contigs
    for i in exon_contig_dict.keys():
        exon_contig_dict[i] = list(set(exon_contig_dict[i]))
    for i in contig_exon_dict.keys():
        contig_exon_dict[i] = list(sorted(np.unique(np.array(contig_exon_dict[i]).astype(int)).astype(str)))
        
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


def find_longest_contig(contig_names,blast_df,contig_file):
    try:
        contigs_file_content = SeqIO.parse(open(contig_file),'fasta')
        contig_info = [(i.id,int(i.description.split(' ')[1])) for i in contigs_file_content if i.id in contig_names]
        longest_contig = sorted(contig_info, key=lambda tup: tup[1],reverse=True)[0][0]
    except:
        contig_lengths = np.array([int(i.split('length_')[1].split('_')[0]) for i in contig_names])
        longest_contig = contig_names[np.where(np.max(contig_lengths))[0][0]]        
    return longest_contig

def get_list_of_valid_exons_and_contigs(exon_contig_dict,loci_with_issues,possible_paralogous,contigs_matching_multiple_exon_dict,keep_duplicates_boolean,keep_paralogs_boolean,outdir,blast_df,contig_file):
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
        print('Warning: Found %i paralogous loci. The longest matching contig for each paralogous locus will be kept, due to the use of the --keep_paralogs flag. It is not recommendable to use paralogous loci for phylogenetic inference!'%len(possible_paralogous))
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
            contig_names = exon_contig_dict[exon]
            contig_names = find_longest_contig(contig_names,blast_df,contig_file)
            valid_contig_names.append(str(contig_names).replace('>',''))
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
                    if orientation == 'minus':
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


#contig_file = '/Users/tobias/GitHub/seqcap_processor/data/test/contigs/1061.fa'
contig_file = '/Users/tobias/GitHub/seqcap_processor/data/test/spades_output/atram_default/contigs.fasta'
#contig_file = '/Users/tobias/GitHub/seqcap_processor/data/test/spades_output/first_try/contigs.fasta'
reference_fasta = '/Users/tobias/GitHub/seqcap_processor_old/data/processed/target_contigs/formatted_reference_library.fasta'
subfolder = '/Users/tobias/GitHub/seqcap_processor/data/test/target_contigs/spades_atram_default/'
keep_duplicates = True
keep_paralogs = False
sample_id = '1061'
seed_length = 11
min_similarity = 0.9
min_length_fraction = 0.9
target_length = 50


# create blast_db
makeblastdb_logfile = os.path.join(subfolder,'makeblastdb.log')
with open(makeblastdb_logfile, 'w') as makeblastdb_out_file:
    makeblastdb_command = ['makeblastdb', '-in', contig_file, '-dbtype', 'nucl']
    run_makeblastdb = subprocess.Popen(makeblastdb_command, stdout=makeblastdb_out_file, stderr=None)
    run_makeblastdb.communicate()

# run blast
blast_cmd = [
    'blastn',
    '-db',
    contig_file,
    '-query',
    reference_fasta,
    '-word_size',
    str(seed_length),
    '-gapopen',
    '100',
    '-gapextend',
    '20',
    '-strand',
    'both',
    '-outfmt',
    '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand'
]

blast_out = os.path.join(subfolder,'%s_all_blast_hits.txt'%sample_id)
blast_err = os.path.join(subfolder,'%s_blast_screen.txt'%sample_id)
#print(blast_cmd)
with open(blast_out, 'w') as out, open(blast_err, 'w') as err:
    blast = subprocess.Popen(blast_cmd, stderr = err, stdout=out)
    blast.wait()


print('Filtering best matches, using min_similarity and min_length_fraction values ...')
selected_matches_out = os.path.join(subfolder,'%s_selected_blast_hits.txt'%sample_id)
blast_hits = pd.read_csv(blast_out,sep='\t',header=None)
blast_hits.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','sstrand']
selected_matches_sim = blast_hits[blast_hits.pident>=min_similarity*100]
selected_matches = selected_matches_sim[selected_matches_sim.length>=min_length_fraction*target_length]
selected_matches.to_csv(selected_matches_out,sep='\t',index=False)

# store the data in dictionaries for convenience
exon_contig_dict, contig_exon_dict, contig_orientation_dict, contig_multi_exon_dict, orientation_df = contigs_matching_exons(selected_matches)
# mark duplicate loci
loci_with_issues, possible_paralogous, contigs_covering_several_loci = find_duplicates(exon_contig_dict,contig_exon_dict)
# remove duplicate loci from the list of targeted loci and contigs
target_contigs = get_list_of_valid_exons_and_contigs(exon_contig_dict,loci_with_issues,possible_paralogous,contig_multi_exon_dict,keep_duplicates,keep_paralogs,subfolder,selected_matches,contig_file)
# load the actual contig sequences
contig_sequences = SeqIO.parse(open(contig_file),'fasta')
# write those contigs that match the reference library to the file
extracted_contig_counter = extract_target_contigs(sample_id,contig_sequences,target_contigs,contig_exon_dict,contig_orientation_dict,subfolder)
print('Extracted %i contigs matching reference exons\n' %extracted_contig_counter)








