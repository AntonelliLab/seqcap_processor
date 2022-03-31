# encoding: utf-8
'''
Create new reference library and map raw reads against the library (reference-based assembly)
'''


#author: Tobias Andermann, tobias.andermann@bioenv.gu.se

import os
import sys
import re
import csv
import glob
import shutil
import argparse
import configparser
import subprocess
import numpy as np
import pickle
import pandas as pd
from Bio import SeqIO

from secapr.utils import CompletePath
from secapr.helpers import CreateDir


# Get arguments
def add_arguments(parser):
    parser.add_argument(
        '--reads',
        required=True,
        action=CompletePath,
        default=None,
        help='Call the folder that contains the trimmed reads, organized in a separate subfolder for each sample. The name of the subfolder has to start with the sample name, delimited with an underscore [_] (default output of clean_reads function).'
    )
    parser.add_argument(
        '--reference_type',
        choices=["alignment-consensus", "sample-specific", "user-ref-lib"],
        default="user-ref-lib",
        help='Please choose which type of reference you want to map the samples to. "alignment-consensus" will create a consensus sequence for each alignment file which will be used as a reference for all samples. This is recommendable when all samples are rather closely related to each other. "sample-specific" will extract the sample specific sequences from an alignment and use these as a separate reference for each individual sample. "user-ref-lib" enables to input one single fasta file created by the user which will be used as a reference library for all samples.'
    )
    parser.add_argument(
        '--reference',
        required=True,
        action=CompletePath,
        default=None,
        help='When choosing "alignment-consensus" or "sample-specific" as reference_type, this flag calls the folder containing the alignment files for your target loci (fasta-format). In case of "user-ref-lib" as reference_type, this flag calls one single fasta file that contains a user-prepared reference library which will be applied to all samples.'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CreateDir,
        default=None,
        help='The output directory where results will be safed.'
    )
    parser.add_argument(
        '--keep_duplicates',
        action='store_true',
        default=False,
        help='Use this flag if you do not want to discard all duplicate reads with Picard.'
    )
    parser.add_argument(
        '--min_coverage',
        type=int,
        default=4,
        help='Set the minimum read coverage. Only positions that are covered by this number of reads will be called in the consensus sequence, otherwise the program will add an ambiguity at this position.'
    )
    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help='Number of computational cores for parallelization of computation.'
    )
    parser.add_argument(
        '--k',
        type=int,
        default=50,
        help='If the part of the read that sufficiently matches the reference is shorter than this threshold, it will be discarded (minSeedLen).'
    )
    parser.add_argument(
        '--w',
        type=int,
        default=21,
        help='Avoid introducing gaps in reads that are longer than this threshold.'
    )
    parser.add_argument(
        '--d',
        type=int,
        default=100,
        help='Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score.'
    )
    parser.add_argument(
        '--r',
        type=float,
        default=1.5,
        help='Trigger re-seeding for a MEM longer than minSeedLen*FLOAT.'
    )
    parser.add_argument(
        '--c',
        type=int,
        default=10000,
        help='Discard a match if it has more than INT occurence in the genome'
    )
    parser.add_argument(
        '--a',
        type=int,
        default=1,
        help='Matching score. Acts as a factor enhancing any match (higher value makes it less conservative = allows reads that have fewer matches, since every match is scored higher).'
    )
    parser.add_argument(
        '--b',
        type=int,
        default=4,
        help='Mismatch penalty. The accepted mismatch rate per read on length k is approximately: {.75 * exp[-log(4) * B/A]}'
    )
    parser.add_argument(
        '--o',
        type=int,
        default=10,
        help='Gap opening penalty'
    )    
    parser.add_argument(
        '--e',
        type=int,
        default=5,
        help='Gap extension penalty'
    )    
    parser.add_argument(
        '--l',
        type=int,
        default=4,
        help='Clipping penalty. During extension, the algorithm keeps track of the best score reaching the end of query. If this score is larger than the best extension score minus the clipping penalty, clipping will not be applied.'
    )    
    parser.add_argument(
        '--u',
        type=int,
        default=2,
        help='Penalty for an unpaired read pair. The lower the value, the more unpaired reads will be allowed in the mapping.'
    )


def create_reference_fasta(reference_folder,alignments):
    # Create a list of fasta files from the input directory
    file_list = [fn for fn in os.listdir(alignments) if fn.endswith(".fasta")]
    reference_list = []
    print('Creating consensus sequences from input alignments...')
    temp_single_refs = []
    for fasta_alignment in file_list:
        sequence_name = re.sub(".fasta","",fasta_alignment)
        orig_aln = os.path.join(alignments,fasta_alignment)
        sep_reference = "%s/%s" %(reference_folder,fasta_alignment)
        reference_list.append(sep_reference)
        cons_cmd = ["cons", "-sequence", orig_aln, "-outseq", sep_reference, "-name", sequence_name, "-plurality", "0.1", "-setcase", "0.1"]
        proc = subprocess.Popen(cons_cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        stderr,stdout = proc.communicate()
        temp_single_refs.append(sep_reference)
    print('Done.')
    reference = os.path.join(reference_folder,"joined_fasta_library.fasta")
    join_fastas = "cat %s/*.fasta > %s" %(reference_folder,reference)
    os.system(join_fastas)
    for i in temp_single_refs:
        os.remove(i)
    return reference


def create_sample_reference_fasta(reference_folder,sample_id,alignments):
    print(("Creating reference library for %s ..." %sample_id))
#    get the sequence header with the correct fasta id and extract sequence
#    store these sequences in separate fasta file for each locus at out_dir/reference_seqs/sample_id
#    header of sequence remains the locus name
#    remove all "-" and "?" in sequence (not in header!!)
    sample_reference_folder = os.path.join(reference_folder,sample_id)
    if not os.path.exists(sample_reference_folder):
        os.makedirs(sample_reference_folder)
    file_list = [fn for fn in os.listdir(alignments) if fn.endswith(".fasta")]
    for fasta_alignment in file_list:
        locus_id = fasta_alignment.replace(".fasta", "")
        sample_reference_fasta = os.path.join(sample_reference_folder,fasta_alignment)
        fasta_sequences = SeqIO.parse(open("%s/%s" %(alignments,fasta_alignment)),'fasta')
        outfile = open(sample_reference_fasta, 'w')
        for fasta in fasta_sequences:
            if fasta.id == sample_id:
                sequence = re.sub('[-,?]','',str(fasta.seq))
                sequence = sequence.upper()
                outfile.write(">%s\n%s\n" %(locus_id,sequence))
        outfile.close()
    reference = os.path.join(sample_reference_folder,"joined_fasta_library.fasta")
    join_fastas = "cat %s/*.fasta > %s" %(sample_reference_folder,reference)
    os.system(join_fastas)
    return reference


def mapping_bwa(forward,backward,reference,sample_id,sample_output_folder, args, log):
    #Indexing
    command1 = ["bwa","index",reference]
    #print(command1)
    bwa_out = os.path.join(log, "bwa_screen_out.txt")
    try:
        with open(bwa_out, 'w') as logfile:
            sp1 = subprocess.Popen(command1, shell=False, stderr = subprocess.STDOUT, stdout=logfile)
            sp1.wait()
    except:
        print(("Running bwa (%s) caused an error." %bwa))
        sys.exit()

    #Mapping
    command2 = ["bwa","mem","-t",str(args.cores),"-k",str(args.k),"-w",str(args.w),"-d",str(args.d),"-r",str(args.r),"-c",str(args.c),"-A",str(args.a),"-B",str(args.b),"-O",str(args.o),"-E",str(args.e),"-L",str(args.l),"-U",str(args.u),"-M",reference,forward,backward]
    """
    Copied from bwa manual (http://bio-bwa.sourceforge.net/bwa.shtml#3):
        -k INT Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [19]
        -w INT Band width. Essentially, gaps longer than INT will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. [100]
        -d INT Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. [100]
        -r FLOAT Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. [1.5]
        -c INT Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [10000] 
        -A INT Matching score. [1]
        -B INT Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [4] 
        -O INT Gap open penalty. [6]
        -E INT Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [1] 
        -L INT Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. [5]
        -U INT Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. [9] 
        -M Mark shorter split hits as secondary (for Picard compatibility). 
    """
    sam_name = os.path.join(sample_output_folder,'%s.sam'%sample_id)
    print ("Mapping...")
    #print(command2)
    with open(sam_name, 'w') as out, open(bwa_out, 'a') as err:
        sp2 = subprocess.Popen(command2, stderr = err, stdout=out)
        sp2.wait()

    #Converting to bam-format with samtools
    print ("Converting to bam...")
    raw_bam = os.path.join(sample_output_folder,"%s_raw.bam" %sample_id)
    command3 = ["samtools","view","-b","-o",raw_bam,"-S",sam_name]
    #print(command3)
    sp3 = subprocess.Popen(command3,stderr=subprocess.PIPE)
    sp3.wait()
    sorted_bam = "%s/%s.sorted.bam" %(sample_output_folder,sample_id)
    command4 = ["samtools", "sort", raw_bam, "-o", sorted_bam]
    #print(command4)
    sp4 = subprocess.Popen(command4)
    sp4.wait()

    #Indexing bam files
    print ("Indexing bam...")
    command5 = ["samtools","index",sorted_bam]
    #print(command5)
    sp5 = subprocess.Popen(command5)
    sp5.wait()
    
    #Remove some big and unnecessary intermediate files
    os.remove(sam_name)
    os.remove(raw_bam)

    return sorted_bam


def mapping_clc(forward,backward,reference,sample_id,sample_output_folder):
    print ("Mapping...")
    cas = "%s/%s.cas" %(sample_output_folder,sample_id)
    command1 = "%s -o %s -d %s -q -p fb ss 100 1000 -i %s %s -l %d -s %d --cpus %d" %(clc_mapper,cas,reference,forward,backward,length,similarity,args.cores)
    os.system(command1)

    print ("Converting to bam...")
    bam = "%s/%s.bam" %(sample_output_folder,sample_id)
    command2 = "%s -a %s -o %s -f 33 -u" %(clc_cas_to_sam,cas,bam)
    os.system(command2)

    print ("Sorting bam...")
    sorted_bam = "%s/%s.sorted.bam" %(sample_output_folder,sample_id)
    #sorted = "%s/%s.sorted" %(sample_output_folder,sample_id)
    command3 = "samtools sort -o %s %s" %(sorted.bam,bam)
    os.system(command3)

    print ("Indexing bam...")
    command4 = "samtools index %s" %(sorted_bam)
    os.system(command4)

    print ("Removing obsolete files...")
    command5 = "rm %s %s" %(cas,bam)
    os.system(command5)

    return sorted_bam


def clean_with_samtools(sample_output_folder,sample_id,sorted_bam,log):
	
    samtools_out = "%s/%s_no_dupls_sorted.bam" %(sample_output_folder,sample_id)
    dupl_log = "%s/%s_dupls.log" %(log,sample_id)
    run_samtools_rmdup = [
        "samtools",
        "rmdup",
        sorted_bam,
        samtools_out
    ]
    print ("Removing duplicate reads with samtools...")
    try:
        with open(os.path.join(log, "samtools_screen_out.txt"), 'w') as log_err_file:
            pi = subprocess.Popen(run_samtools_rmdup, stderr=log_err_file)
            pi.communicate()
    except OSError:
        print('Not enough reads mapped to reference in order to run samtools rmdup. Try using the "--keep_duplicates" flag in order to avoid the use of duplicate cleaning.') 
        quit()
    #print ("Duplicates successfully removed.")
    # Cleaning up a bit
    has_duplicates = "%s/including_duplicate_reads" %sample_output_folder
    if not os.path.exists(has_duplicates):
        os.makedirs(has_duplicates)
    mv_duplicates_1 = 'mv %s/*.bam %s' %(sample_output_folder,has_duplicates)
    mv_duplicates_2 = 'mv %s/*.bam.bai %s' %(sample_output_folder,has_duplicates)
    mv_final = 'mv %s/*_no_dupls_sorted.bam %s' %(has_duplicates,sample_output_folder)
    os.system(mv_duplicates_1)
    os.system(mv_duplicates_2)
    os.system(mv_final)

    print ("Indexing duplicate-free bam...")
    index_cleaned_bam = "samtools index %s" %(samtools_out)
    os.system(index_cleaned_bam)
    dupl_bam_name = sorted_bam.split('/')[-1]
    dupl_out = os.path.join(has_duplicates,dupl_bam_name)
    return samtools_out, dupl_out


def clean_with_picard(sample_output_folder,sample_id,sorted_bam,log):

    picard_out = "%s/%s_no_dupls_sorted.bam" %(sample_output_folder,sample_id)
    dupl_log = "%s/%s_dupls.log" %(log,sample_id)
    run_picard = [
        "picard",
        "MarkDuplicates",
        "I=%s" %sorted_bam,
        "O=%s" %picard_out,
        "M=%s" %dupl_log,
        "REMOVE_DUPLICATES=true",
        "VALIDATION_STRINGENCY=LENIENT"
    ]
    print ("Removing duplicate reads with Picard...")
    try:
        with open(os.path.join(log, "picard_screen_out.txt"), 'w') as log_err_file:
            pi = subprocess.Popen(run_picard, stderr=log_err_file)
            pi.communicate()
    except OSError:
        print('Not enough reads mapped to reference in order to run Picard. Try using the "--keep_duplicates" flag in order to avoid the use of Picard.') 
        quit()
    print ("Duplicates successfully removed.")
    # Cleaning up a bit
    has_duplicates = "%s/including_duplicate_reads" %sample_output_folder
    if not os.path.exists(has_duplicates):
        os.makedirs(has_duplicates)
    mv_duplicates_1 = 'mv %s/*.bam %s' %(sample_output_folder,has_duplicates)
    mv_duplicates_2 = 'mv %s/*.bam.bai %s' %(sample_output_folder,has_duplicates)
    mv_final = 'mv %s/*_no_dupls_sorted.bam %s' %(has_duplicates,sample_output_folder)
    os.system(mv_duplicates_1)
    os.system(mv_duplicates_2)
    os.system(mv_final)

    print ("Indexing Picard-cleaned bam...")
    index_picard_bam = "samtools index %s" %(picard_out)
    os.system(index_picard_bam)
    dupl_bam_name = sorted_bam.split('/')[-1]
    dupl_out = os.path.join(has_duplicates,dupl_bam_name)
    return picard_out, dupl_out


def bcf_cons(pileup,out_fasta,cov):
    with open(pileup) as f:
        content = [x.strip('\n') for x in f.readlines()]

    loci = []
    seq_dict = {}
    for line in content:
        if '#' not in line:
            # Split the tab delimited lines into their segments
            element = line.split('\t')
            # Create a list of all loci names
            seq_name = element[0]
            if seq_name not in loci:
                loci.append(seq_name)
            # By default call every position a uncertainty
            basecall = "N"

            # Turn all lower case values in upper case
            sample = element[4].upper()
            # make a directory with all different basecalls and count their occurences
            calls = dict((letter,sample.count(letter)) for letter in set(sample))

            # The basecall in the reference
            reference = element[2]
            # These characters signal a match with the reference
            match_ref = "." ","
            # List of base characters
            bases = "A" "G" "C" "T"

            # find out how many agreements with reference are among the basecalls. These are all . and , basecalls listed in match_ref.
            # reset the counter before every round (every line in file aka position in sequence)
            list_matches = 0
            for key,value in list(calls.items()):
                if key in match_ref:
                    list_matches += value
            if list_matches >= cov:
                basecall = reference

            # find if there are any well supported SNPs and make the most prominent call
            for key in sorted(calls, key=calls.get, reverse=True):
                if key in bases:
                    if int(calls[key]) >= cov:
                        if int(calls[key]) >= list_matches:
                            basecall = key
                            break
            # add the final basecall to the dictionary and to the respective key if it already exists, otherwise create new key
            seq_dict.setdefault(element[0],[])
            seq_dict[element[0]].append(basecall)
    # Join all basecalls for each key (=locus) into one sequence and deposit in new dictionary
    concat_basecalls = {}
    for key, value in list(seq_dict.items()):
        concat_basecalls[key] = "".join(value)
    #print concat_basecalls

    with open(out_fasta, "w") as f:
        for k, v in list(concat_basecalls.items()):
            f.write(">" + k+ "\n")
            f.write(v+ "\n")
    return out_fasta


def bam_consensus(reference,bam_file,name_base,out_dir,min_cov,cons_method = 'custom'):
    print ("Generating a consensus sequence from bam-file...")
    # Creating consensus sequences from bam-files

    # This is the custom written consensus building function
    if cons_method == 'custom':
        mpileup = [
            "samtools",
            "mpileup",
            "-A",
            "-f",
            reference,
            bam_file
        ]
        mpileup_file = os.path.join(out_dir, "%s.bcf" %name_base)
        with open(mpileup_file, 'w') as mpileupfile:
            mp = subprocess.Popen(mpileup, stdout=mpileupfile, stderr=subprocess.PIPE)
            mp.communicate()
            mp.wait()
        fasta_file = os.path.join(out_dir,"%s_temp.fasta" %name_base)
        fasta_file = bcf_cons(mpileup_file,fasta_file,min_cov)
    
    # This is the bcftools version which sometimes acts a bit funny
    else:
        # 1. mpilup__________
        mpileup = [
            "samtools",
            "mpileup",
            "-u",
            "-f",
            reference,
            bam_file
        ]
        mpileup_file = os.path.join(out_dir, "%s.mpileup" %name_base)
        with open(mpileup_file, 'w') as mpileupfile:
            mp = subprocess.Popen(mpileup, stdout=mpileupfile, stderr=subprocess.PIPE)
            mp.communicate()
            mp.wait()

        # 2. bcftools__________
        vcf_file = os.path.join(out_dir, "%s.vcf" %name_base)
        #bcf_cmd = [
        #    "bcftools",
        #    "view",
        #    "-c",
        #    "-g",
        #    mpileup_file
        #]
        bcf_cmd = [
            "bcftools",
            "call",
            "-c",
            mpileup_file
        ]
        with open(vcf_file, 'w') as vcffile:
            vcf = subprocess.Popen(bcf_cmd, stdout=vcffile, stderr=subprocess.PIPE)
            vcf.communicate()
            vcf.wait()

        # 3. vcfutils__________
        fq_file = os.path.join(out_dir, "%s.fq" %name_base)
        vcfutils_cmd = [
            "vcfutils.pl",
            "vcf2fq",
            "-d",
            str(min_cov),
            vcf_file
        ]
        with open(fq_file, 'w') as fqfile:
            fq = subprocess.Popen(vcfutils_cmd, stdout=fqfile)
            fq.communicate()
            fq.wait()
        
        # 4. seqtk__________
        # Converting fq into fasta files
        fasta_file = os.path.join(out_dir,"%s_temp.fasta" %name_base)
        make_fasta = [
            "seqtk",
            "seq",
            "-a",
            fq_file,
        ]
        with open(fasta_file, 'w') as fastafile:
            fasta = subprocess.Popen(make_fasta, stdout=fastafile)
            fasta.communicate()
            fasta.wait()

    # Independently of the chosen consensus building method, the following brings the resulting fasta file into the correct format for further processing
    # Create a new fasta file for final fasta printing
    final_fasta_file = os.path.join(out_dir,"%s.fasta" %name_base)

    if "allele_0" in name_base:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        if '_no_dupls' in name_base:
            sample_id = name_base.split("_no_dupls_allele_0")[0]
        elif '_selected_loci' in name_base:
            sample_id = name_base.split("_selected_loci")[0]
        else:
            sample_id = name_base.split("_allele_0")[0]
        with open(final_fasta_file, "w") as out_file:
            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                name = re.sub('_consensus_sequence','',name)
                name = re.sub('_\(modified\)','',name)
                #name = re.sub(r'(.*)',r'\1_%s_0 |\1' %sample_id ,name)
                sequence = re.sub('[a,c,t,g,y,w,r,k,s,m,n,Y,W,R,K,S,M]','N',sequence)
                out_file.write(">%s_%s_0 |%s\n%s\n" %(name,sample_id,name,sequence))
    
    elif "allele_1" in name_base:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        if '_no_dupls' in name_base:
            sample_id = name_base.split("_no_dupls_allele_1")[0]
        elif '_selected_loci' in name_base:
            sample_id = name_base.split("_selected_loci")[0]
        else:
            sample_id = name_base.split("_allele_1")[0]
        with open(final_fasta_file, "w") as out_file:
            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                name = re.sub('_consensus_sequence','',name)
                name = re.sub('_\(modified\)','',name)
                #name = re.sub(r'(.*)',r'\1_%s_1 |\1' %sample_id ,name)
                sequence = re.sub('[a,c,t,g,y,w,r,k,s,m,n,Y,W,R,K,S,M]','N',sequence)
                out_file.write(">%s_%s_1 |%s\n%s\n" %(name,sample_id,name,sequence))

    else:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        if '_no_dupls' in name_base:
            sample_id = name_base.split("_no_dupls_bam_consensus")[0]
        elif '_selected_loci' in name_base:
            sample_id = name_base.split("_selected_loci")[0]        
        else:
            sample_id = name_base.split("_bam_consensus")[0]
        with open(final_fasta_file, "w") as out_file:
            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                name = re.sub('_consensus_sequence','',name)
                name = re.sub('_\(modified\)','',name)
                #name = re.sub(r'(.*)',r'\1_%s |\1' %sample_id ,name)
                sequence = re.sub('[a,c,t,g,y,w,r,k,s,m,n,Y,W,R,K,S,M]','N',sequence)
                out_file.write(">%s_%s |%s\n%s\n" %(name,sample_id,name,sequence))
        tmp_folder = "%s/tmp" %out_dir
        if not os.path.exists(tmp_folder):
            os.makedirs(tmp_folder)

        if cons_method == 'custom':
            bcf = glob.glob('%s/*.bcf' %out_dir)[0]
            shutil.move(bcf,tmp_folder)
        else:
            mpileup = glob.glob('%s/*.mpileup' %out_dir)[0]
            vcf = glob.glob('%s/*.vcf' %out_dir)[0]
            shutil.move(vcf,tmp_folder)
            shutil.move(mpileup,tmp_folder)

    os.remove(fasta_file)
    if not cons_method == 'custom':
        os.remove(fq_file)

    return final_fasta_file


def join_fastas(out_dir,sample_out_list):
    allele_fastas = []
    unphased_fastas = []
    allele = False
    for folder in sample_out_list:
        for file in os.listdir(folder):
            if file.endswith('.fasta'):
                fasta = os.path.join(folder,file)
                if "allele_0" in fasta:
                    allele_fastas.append(fasta)
                    allele = True
                elif "allele_1" in fasta:
                    allele_fastas.append(fasta)
                    allele = True
                elif fasta.endswith('bam_consensus.fasta'):
                    unphased_fastas.append(fasta)

    if allele == True:
        joined_allele_fastas = "%s/joined_allele_fastas.fasta" %out_dir
        with open(joined_allele_fastas, 'w') as outfile:
            for fname in allele_fastas:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
    if not allele:
        joined_unphased_fastas = "%s/joined_unphased_fastas.fasta" %out_dir
        with open(joined_unphased_fastas, 'w') as outfile:
            for fname in unphased_fastas:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)    


def get_bam_path_dict(input_dir):
    type_input = ''
    subdirs = os.listdir(input_dir)
    sample_bam_dict = {}
    for subd in subdirs:
        if subd.endswith('_remapped'):
            type_input = 'unphased'
            sample_id = subd.split('_remapped')[0]
            bam = '%s/%s/%s*sorted.bam' %(input_dir,subd,sample_id)
            target_files = glob.glob(bam)
            unphased_bam = target_files[0]
            if len(target_files) > 1:
                print('Found multiple files matching the search, but there should only be one bam file containing the re-mapped reads. Please remove any non-relevant bam-files from target directory and run this function again.')
                print(target_files)
                exit()
            else:
                sample_bam_dict.setdefault(subd,[])
                sample_bam_dict[subd].append(unphased_bam)
        elif subd.endswith('_phased'):
            type_input = 'phased'
            sample_id = subd.split('_phased')[0]
            bam = '%s/%s/phased_bam_files/%s*sorted_allele_[0,1].bam' %(input_dir,subd,sample_id)
            target_files = glob.glob(bam)
            allele_0_bam = target_files[0]
            allele_1_bam = target_files[1]
            if len(target_files) > 2:
                print('Found multiple files matching the search, but there should only be one bam file per allele containing the re-mapped reads. Please remove any non-relevant bam-files from target directory and run this function again.')
                print(target_files)
                exit()
            else:
                sample_bam_dict.setdefault(subd,[])
                sample_bam_dict[subd].append(allele_0_bam)
                sample_bam_dict[subd].append(allele_1_bam)
    return sample_bam_dict, type_input


def get_bam_read_cov(bam,output_folder):
    bam_name = bam.split("/")[-1]
    sample_base = bam_name.split(".bam")[0]
    if sample_base.endswith('_no_dupls_sorted'):
        sample_base = sample_base.split("_no_dupls_sorted")[0]
    elif sample_base.endswith('_sorted'):
        sample_base = sample_base.split("_sorted")[0]
    else:
        sample_base = sample_base.split("_")[0]
    sample_base = sample_base.split(".")[0]
    print(('Reading read-depth info for %s.........' %sample_base))
    sample_dir = '/'.join(bam.split("/")[:-1])
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    get_read_depth = ["samtools", "depth", bam]
    read_depth_file = os.path.join(sample_dir,"%s_read_depth_per_position.txt" %sample_base)
    with open(read_depth_file, 'w') as logfile:
        sp1 = subprocess.Popen(get_read_depth, shell=False, stderr = subprocess.STDOUT, stdout=logfile)
        sp1.wait()
    return sample_dir,read_depth_file

def get_bam_read_count(bam_path):
    read_count_cmd = ['samtools',
                      'view',
                      '-F',
                      '0x904',
                      '-c',
                      bam_path
                      ]
    fasta = subprocess.Popen(read_count_cmd,stdout=subprocess.PIPE)
    read_count = fasta.stdout.read()
    read_count = int(read_count.strip().decode('ascii'))
    return(read_count)


def get_complete_loci_list(subfolder_file_dict):
    print('Generating locus database.........')
    locus_list = []
    for subfolder in subfolder_file_dict:
        read_depth_file = subfolder_file_dict[subfolder]
        with open(read_depth_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            reader = list(reader)
            for row in reader:
                locus_name = row[0]
                if locus_name not in locus_list:
                    locus_list.append(locus_name)
    return sorted(locus_list)


def summarize_read_depth_files(subfolder,read_depth_file,complete_locus_list,locus_dict_all_samples,sample_list,reference):
    sample_id = read_depth_file.split('/')[-1].split('_read_depth_per_position.txt')[0]
    sample_list.append(sample_id)
    print(('Calculating coverage for all loci from bam files for %s.........' %sample_id))
    reference_library = SeqIO.parse(reference, "fasta")
    locus_original_length_dict = {}
    for reference in reference_library:
        reference_locus = reference.name
        reference_sequence = reference.seq
        locus_length = len(str(reference_sequence))
        locus_original_length_dict.setdefault(reference_locus,locus_length)
    sample_loci_dict = {}
    with open(read_depth_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        reader = list(reader)
        for row in reader:
            locus = row[0]
            locus_list = locus.split("_")[:3]
            locus_name = "_".join(locus_list)
            position = row[1]
            coverage = int(row[2])
            sample_loci_dict.setdefault(locus_name,[])
            sample_loci_dict[locus_name].append(coverage)
    for locus_name in complete_locus_list:
        if locus_name in sample_loci_dict:
            avg_read_depth = float(sum(sample_loci_dict[locus_name]))/float(locus_original_length_dict[locus_name])
            locus_dict_all_samples.setdefault(locus_name,[])
            locus_dict_all_samples[locus_name].append(avg_read_depth)
        else:
            avg_read_depth = 0.0
            locus_dict_all_samples.setdefault(locus_name,[])
            locus_dict_all_samples[locus_name].append(avg_read_depth)
    return locus_dict_all_samples 


def main(args):
    print('\n')
    mapper = 'bwa'
    # Set working directory
    out_dir = args.output
    #if not os.path.exists(out_dir):
    #    os.makedirs(out_dir)
    #else:
    #    raise IOError("The directory {} already exists.  Please check and remove by hand.".format(out_dir))
    # Get other input variables
    alignments = args.reference
    reads = args.reads
    min_cov = args.min_coverage
    reference = ''
    reference_folder = "%s/reference_seqs" %out_dir
    sample_out_list = []
    total_read_count_dict = {}
    if not os.path.exists(reference_folder):
        os.makedirs(reference_folder)
    if args.reference_type == "user-ref-lib":
        reference = args.reference
        manage_reference = "cp %s %s" %(reference,reference_folder)
        os.system(manage_reference)
    elif args.reference_type == "alignment-consensus":
        reference = create_reference_fasta(reference_folder,alignments)
    for subfolder in [ name for name in os.listdir(reads) if os.path.isdir(os.path.join(reads, name)) ]:
        subfolder_path = os.path.join(reads,subfolder)
        sample_id = subfolder
        print(('\n' + "#" * 50))
        print(("Processing sample %s" % sample_id))
        if args.reference_type == "sample-specific":
            reference = create_sample_reference_fasta(reference_folder,sample_id,alignments)
        # Safe the sample specific reference as a pickle file for downstream processing
        sample_output_folder = "%s/%s_remapped" %(out_dir,sample_id)
        sample_out_list.append(sample_output_folder)
        if not os.path.exists(sample_output_folder):
            os.makedirs(sample_output_folder)
        tmp_folder = "%s/tmp" %sample_output_folder
        if not os.path.exists(tmp_folder):
            os.makedirs(tmp_folder)
        pickle_path = os.path.join(tmp_folder,'%s_reference.txt' %sample_id)
        np.savetxt(pickle_path,np.array([reference]),fmt='%s')
        # with open(pickle_path, 'wb') as handle:
        #     pickle.dump(reference, handle, protocol=pickle.HIGHEST_PROTOCOL)


        # Loop through each sample-folder and find read-files
        read1_files = sorted(glob.glob(os.path.join(subfolder_path,'_clean-READ1.fastq.gz')))
        read2_files = sorted(glob.glob(os.path.join(subfolder_path, '_clean-READ2.fastq.gz')))
        for i,forward in enumerate(read1_files):
            backward = read2_files[i]


            # samtools merge joined.sam tpella5.sam tpella9.sam


            sorted_bam = ""
            log = os.path.join(sample_output_folder,'log')
            if not os.path.exists(log):
                os.makedirs(log)
            if mapper == "bwa":
                sorted_bam = mapping_bwa(forward,backward,reference,sample_id,sample_output_folder,args,log)
            if not args.keep_duplicates:
                sorted_bam, dupl_bam = clean_with_samtools(sample_output_folder,sample_id,sorted_bam,log)
            name_stem = '%s_bam_consensus' %sample_id
            bam_consensus_file = bam_consensus(reference,sorted_bam,name_stem,sample_output_folder,min_cov)
            if not args.keep_duplicates:
                dupl_output_folder = ('/').join(dupl_bam.split('/')[:-1])
                dupl_name_stem = '%s_with_duplicates_bam_consensus' %sample_id
                bam_consensus_with_duplicates = bam_consensus(reference,dupl_bam,dupl_name_stem,dupl_output_folder,min_cov)
    join_fastas(out_dir,sample_out_list)
    # create file with read-coverage overview
    print(('\n'+"#" * 50))
    sample_bam_dict, input_type = get_bam_path_dict(out_dir)
    # currently only available for unphased data
    if input_type == 'unphased':
        subfolder_list = []
        subfolder_file_dict = {}
        reference_file_dict = {}
        # iterating through samples
        for key in sample_bam_dict:
            if key.endswith('_remapped'):
                sample = key.split('_remapped')[0]
                path = os.path.join(out_dir,key)
                path2 = os.path.join(path,'tmp')
                reference_pickle = os.path.join(path2,'%s_reference.txt' %sample)
                reference_file_dict.setdefault(sample,reference_pickle)        
                bam = sample_bam_dict[key][0]
                total_read_count = get_bam_read_count(bam)
                total_read_count_dict.setdefault(sample,total_read_count)
                sample_dir, read_depth_file = get_bam_read_cov(bam,out_dir)
                subfolder_file_dict.setdefault(sample_dir,read_depth_file)
        locus_list = get_complete_loci_list(subfolder_file_dict)
        locus_dict_all_samples = {}
        sample_id_list = []
        for subfolder in subfolder_file_dict:
            sample = subfolder.split('/')[-1].split('_remapped')[0]
            reference = ''
            reference_pickle = reference_file_dict[sample]
            reference = str(np.loadtxt(reference_pickle,dtype=str))
            # with open(reference_pickle, 'rb') as handle:
            #     reference = pickle.load(handle)
            read_depth_file = subfolder_file_dict[subfolder]
            locus_dict_all_samples = summarize_read_depth_files(subfolder,read_depth_file,locus_list,locus_dict_all_samples,sample_id_list,reference)
        output_dict = {}
        # Create a separate list for each column in the final csv file
        final_locus_list = ['locus']
        for locus in locus_dict_all_samples:
            final_locus_list.append(locus)
        output_dict.setdefault('column1',final_locus_list)
        for sample in sample_id_list:
            index = sample_id_list.index(sample)
            key_int = 2+index
            key_name = 'column%i' %key_int
            output_dict.setdefault(key_name,[sample])
            for locus in locus_dict_all_samples:
                read_depth = locus_dict_all_samples[locus][index]
                output_dict[key_name].append(read_depth)
        final_data_list = []
        for column in sorted(output_dict):
            data = output_dict[column]
            final_data_list.append(data)
        output = open("%s/average_cov_per_locus.txt" %out_dir, "w")
        outlog=csv.writer(output, delimiter='\t')
        transformed_data = list(zip(*final_data_list))
        for row in transformed_data:
            outlog.writerow(row)
# The following lines don't work because of some strange issue with reading the text file with pandas:
#        txt_path = os.path.join(out_dir,"average_cov_per_locus.txt")
#        final_overview = pd.read_csv(txt_path, sep='\t')
#        no_of_loci = len(final_overview)
#        new_dict = {}
#        for column in final_overview:
#            if not column in ['sum_per_locus','locus']:
#                average = sum(final_overview[column])/no_of_loci
#                new_dict.setdefault(column,average)
#        avg_read_cov_across_all_loci = pd.DataFrame.from_dict(new_dict, orient='index', dtype=None)
#        avg_read_cov_across_all_loci.columns = ['average']
#        avg_read_cov_across_all_loci.sort_values(by='average',ascending=False).to_csv('%s/average_read_coverage_across_all_loci_per_sample.txt' %out_dir, sep = '\t', index = True,header=False)
        #print(total_read_count_dict)
        try:
            pickle_in = os.path.join(args.reference,'.secapr_files/sequence_origin.pickle')
            with open(pickle_in, 'rb') as handle:
                sequence_origin = pickle.load(handle)
            stats_file_path = '/'.join(sequence_origin.split('/')[:-1])
            stats_file = os.path.join(stats_file_path,'sample_stats.txt')
            stats_df = pd.read_csv(stats_file,sep='\t')
            new_df = stats_df.copy()
            new_df['reads_on_target'] = [0]*len(new_df)
            for key in list(total_read_count_dict.keys()):
                index_row = new_df[new_df['sample_id'] == int(key)].index
                new_df.iloc[index_row,-1] = int(total_read_count_dict[key])
            new_df.to_csv(os.path.join(args.output,'sample_stats.txt'),sep='\t',index=False)
        except:
            alternative_stats_file = 'NONE'
            print(('INFO: Stats file for samples could not be found (not available when using --reference_type user-ref-lib). Stats about read_coverage are instead written to file %s' %alternative_stats_file))
