# encoding: utf-8
# author: Tobias Andermann, tobias.andermann@bioenv.gu.se
'''
Extract the contigs that match the reference database
'''


import re
import os
import sys
import glob
import logging
import argparse
import subprocess
import numpy as np
import pandas as pd
#import pylab as plt
from secapr.helpers import is_dir, is_file, FullPaths
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#import pdb
log = logging.getLogger(__name__)

def add_arguments(parser):
    parser.add_argument(
        '--contigs',
        required=True,
        type=is_dir,
        action=FullPaths,
        help="The directory containing the assembled contigs in fasta format."
    )
    parser.add_argument(
        '--reference',
        required=True,
        type=is_file,
        action=FullPaths,
        help="The fasta-file containign the reference sequences (probe-order-file)."
    )
    parser.add_argument(
        '--output',
        required=True,
        action=FullPaths,
        help="The directory in which to store the extracted target contigs and lastz results."
    )
    parser.add_argument(
        '--target_length',
        default=50,
        type=int,
        help="The required length of the matching sequence stretch between contigs and target sequences. This does not have to be a perfect match but can be adjusted with the --min_identity flag [default=50]."
    )
    parser.add_argument(
        '--min_identity',
        default=90,
        type=int,
        help="The minimum percent identity required for a match [default=90]."
    )
    parser.add_argument(
        '--seed_length',
        default=11,
        type=int,
        help="Length of initial seed sequence for finding BLAST matches. The seed has to be a perfect match between a given contig and a reference locus (default=11)."
    )
    parser.add_argument(
        "--remove_multilocus_contigs",
        action='store_true',
        default=False,
        help="Drop those contigs that span across multiple exons.",
    )
    parser.add_argument(
        "--keep_paralogs",
        action='store_true',
        default=False,
        help="Keep loci with signs of paralogy (multiple contig matches). The longest contig will be selected for these loci.",
    )
    parser.add_argument(
        '--disable_stats',
        action='store_true',
        default=False,
        help='Disable generation of stats (can be necessary because previous stats files can\'t be found if files are used that were not previously processed with SECAPR)'
    )



def contig_count(contig):
    """Return a count of contigs from a fasta file"""
    return sum([1 for line in open(contig, 'r').readlines() if line.startswith('>')])


def new_get_probe_name(header, regex):
    match = re.search(regex, header)
    #print match
    return match.groups()[0]


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
    duplicate_info.to_csv(os.path.join(outdir,'info_contigs_spanning_multiple_loci.txt'),header=False,sep="\t")
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
                            out_file.write(new_fasta.format('fasta-2line'))
                            sample_file.write(new_fasta.format('fasta-2line'))
                            printed_contigs_counter += 1                         
                    else:        
                        header = '%s_%s |%s' %(contig_exon_dict[fasta.id][0],sample_id,contig_exon_dict[fasta.id][0])
                        new_fasta = SeqRecord(seq, id=header, name='', description='')
                        out_file.write(new_fasta.format('fasta-2line'))
                        sample_file.write(new_fasta.format('fasta-2line'))
                        printed_contigs_counter += 1 
    return printed_contigs_counter



def main(args):
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    # Write new reference file with numbers as fasta-headers
    reference_fasta = open(args.reference,'r')
    new_fasta = os.path.join(args.output,'formatted_reference_library.fasta')
    new_reference_fasta = open(new_fasta,'w')
    if args.remove_multilocus_contigs:
        keep_duplicates = False
    else:
        keep_duplicates = True
    counter_sequence_dict = {}
    counter = 0
    for line in reference_fasta:
        if line.startswith('>'):
            old_header = line.replace('>','').strip()
            counter_sequence_dict.setdefault(counter,old_header)
            new_header = '>%i\n' %(counter)
            counter += 1
            new_reference_fasta.write(new_header)
        else:
            new_reference_fasta.write(line)
    new_reference_fasta.close()
    # write the translation dictionary between new numerical identifiers and previous fasta headers to file
    header_info_file = os.path.join(args.output,'reference_fasta_header_info.txt')
    header_info = pd.DataFrame.from_dict(counter_sequence_dict, orient='index')
    header_info.to_csv(header_info_file,sep='\t',header=False,index=True)
    # get the fasta headers from the new formatted reference file
    exons = [seq.id for seq in SeqIO.parse(open(new_fasta, 'r'), 'fasta')]
    sorted_exon_list = list(exons)
    # Get the paths to the contig fasta files for all samples
    fasta_files = glob.glob(os.path.join(args.contigs, '*.fa*'))
    sample_ids = [os.path.basename(fasta).split('.')[0] for fasta in fasta_files]
    # Create a dataframe filled with 0's
    contig_match_df = pd.DataFrame(index=sorted_exon_list,columns=sample_ids)
    for locus in sorted_exon_list:
        contig_match_df.loc[locus] = [0]*len(sample_ids)
    # Print some log screen output
    log.info("Processing contig data")
    log.info("{}".format("-" * 65))
    # Start processing by iterating through contig files (=samples)
    for contig_file in sorted(fasta_files):
        # Get the name of the sample
        sample_id = os.path.basename(contig_file).split('.')[0]#.replace('-', "_")
        # Make subfolder for each sample
        subfolder = os.path.join(args.output,sample_id)
        if not os.path.exists(subfolder):
            os.makedirs(subfolder)        
        # Print some stats to screen
        total_count_of_contig = contig_count(contig_file)
        print('%s:\n'%sample_id,'Total contigs: %i\nSearching for contigs with matches in reference database.'%total_count_of_contig)
        # Blast the the contigs against the reference file
        # create blast_db
        reference_lib = os.path.join(subfolder,os.path.basename(contig_file))
        mv_contig = 'cp %s %s'%(contig_file,reference_lib)
        os.system(mv_contig)
        makeblastdb_logfile = os.path.join(subfolder,'makeblastdb.log')
        with open(makeblastdb_logfile, 'w') as makeblastdb_out_file:
            makeblastdb_command = ['makeblastdb', '-in', reference_lib, '-dbtype', 'nucl']
            run_makeblastdb = subprocess.Popen(makeblastdb_command, stdout=makeblastdb_out_file, stderr=None)
            run_makeblastdb.communicate()
        # run blast
        blast_cmd = [
            'blastn',
            '-db',
            reference_lib,
            '-query',
            new_fasta,
            '-word_size',
            str(args.seed_length),
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
        selected_matches_sim = blast_hits[blast_hits.pident>=args.min_identity]
        selected_matches = selected_matches_sim[selected_matches_sim.length>=args.target_length]
        selected_matches.to_csv(selected_matches_out,sep='\t',index=False)

        # store the data in dictionaries for convenience
        exon_contig_dict, contig_exon_dict, contig_orientation_dict, contig_multi_exon_dict, orientation_df = contigs_matching_exons(selected_matches)
        # mark duplicate loci
        loci_with_issues, possible_paralogous, contigs_covering_several_loci = find_duplicates(exon_contig_dict,contig_exon_dict)
        # remove duplicate loci from the list of targeted loci and contigs
        target_contigs = get_list_of_valid_exons_and_contigs(exon_contig_dict,loci_with_issues,possible_paralogous,contig_multi_exon_dict,keep_duplicates,args.keep_paralogs,subfolder,selected_matches,contig_file)
        # load the actual contig sequences
        contig_sequences = SeqIO.parse(open(contig_file),'fasta')
        # write those contigs that match the reference library to the file
        extracted_contig_counter = extract_target_contigs(sample_id,contig_sequences,target_contigs,contig_exon_dict,contig_orientation_dict,subfolder)
        # Fill the extracted target contig into the dataframe
        for contig in target_contigs:
            for exon in contig_exon_dict[contig]:
                contig_match_df.loc[exon,sample_id] = 1        
        print('Extracted %i contigs matching reference exons\n' %extracted_contig_counter)
        log.info("{}".format("-" * 65))

    contig_match_df.to_csv(os.path.join(args.output,'match_table.txt'),sep='\t',index=True,encoding='utf-8')
    # Print summary stats
    table = pd.read_csv(os.path.join(args.output,'match_table.txt'), delimiter = '\t',index_col=0)
    with open(os.path.join(args.output,'summary_stats.txt'), "w") as out_file:
        out_file.write('Total number of samples: %i\nTotal number of targeted exons: %i\n\n'%(len(table.columns),len(table)))
        complete_loci_counter = 0
        for locus in table.iterrows():
            if sum(locus[1]) == len(locus[1]):
                complete_loci_counter += 1
        out_file.write('%i exons are shared by all samples.\n\n'%complete_loci_counter)
        count_list = []
        sample_count_dict = {}
        for column in table.columns:
            sample_count_dict.setdefault(column,sum(table[column]))
            count_list.append(sum(table[column]))
            out_file.write('%s: %i extracted contigs\n'%(column,sum(table[column])))
        out_file.write('mean: %f stdev: %f'%(np.mean(count_list),np.std(count_list)))
    try:
        stats_df = pd.read_csv(os.path.join(args.contigs,'sample_stats.txt'),sep='\t')
        new_df = stats_df.copy()
        new_df['target_contigs'] = [0]*len(new_df)
        for key in list(sample_count_dict.keys()):
            index_row = new_df[new_df['sample'].apply(str) == str(key)].index 
            new_df.iloc[index_row,-1] = int(sample_count_dict[key])
        new_df.to_csv(os.path.join(args.output,'sample_stats.txt'),sep='\t',index=False)
    except:
        print('No stats-file found in contig-folder. Instead printing stats to screen:')
        for key in list(sample_count_dict.keys()):
            print('Sample %s: %i contigs extracted.'%(key,int(sample_count_dict[key])))



#    # Make plot of match table
#    # Get the data from the df
#    sample_labels = contig_match_df.columns
#    locus_labels = np.array(contig_match_df.index)
#    data = np.matrix(contig_match_df).T
#    # Define the figure and plot to png file
#    fig, ax = plt.subplots()
#    mat = ax.imshow(data, cmap='GnBu', interpolation='nearest')
#    plt.xticks(range(data.shape[1])[::20], locus_labels[::20],fontsize=3)
#    plt.yticks([0],fontsize=0)
#    #plt.yticks(range(data.shape[0])[::3], sample_labels[::3],fontsize=3)
#    plt.xticks(rotation=90)
#    plt.xlabel('exon',fontsize=3)
#    plt.ylabel('sample',fontsize=3)
#    fig.savefig(os.path.join(args.output,'contig_exon_matrix.png'), dpi = 500)
