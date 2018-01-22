#author: Tobias Hofmann, tobiashofmann@gmx.net

'''
Assemble trimmed Illumina read files (fastq)
'''

import os
import sys
import re
import glob
import shutil
import argparse
import commands
import subprocess
import pandas as pd
import numpy as np

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Input %%%

# Complete path function
class CompletePath(argparse.Action):
    """give the full path of an input file/folder"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def add_arguments(parser):
    parser.add_argument(
        '--input',
        required=True,
        action=CompletePath,
        default=None,
        help='Call the folder that contains the trimmed reads, organized in a separate subfolder for each sample. The name of the subfolder has to start with the sample name, delimited with an underscore [_]'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CompletePath,
        default=None,
        help='The output directory where results will be saved'
    )
#    parser.add_argument(
#        '--assembler',
#        choices=["abyss","trinity"],
#        default="abyss",
#        help="""The assembler to use (default = abyss)."""
#    )
    parser.add_argument(
        '--kmer',
        type=int,
        default=35,
        help='Set the kmer value'
    )
    parser.add_argument(
        '--contig_length',
        type=int,
        default=200,
        help='Set the minimum contig length for the assembly. Contigs that are shorter than this threshold will be discarded. [Only available for Trinity assembler]'
    )
    parser.add_argument(
        '--single_reads',
        action='store_true',
        default=False,
        help='Use this flag if you additionally want to use single reads for the assembly'
    )
    parser.add_argument(
        '--disable_stats',
        action='store_true',
        default=False,
        help='Use this flag if you want to disabel generation of stats (can be necessary because previous stats files can\'t be found if reads are used that were not previously processed with SECAPR) '
    )
    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help='For parallel processing you can set the number of cores you want to run the assembly on.'
    )

def main(args):
    # Set working directory
    out_folder = args.output
    out_dir = "%s/stats" %out_folder
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Get all the other input variables
    input_folder = args.input
    min_length = args.contig_length
    #trinity = args.trinity
    cores = args.cores
    #abyss = args.abyss
    kmer = args.kmer
    #assembler = args.assembler
    assembler = 'abyss'
    home_dir = os.getcwd()
    sample_contig_count_dict = {}
    if cores > 1:
        print("Running %s parallel on %d cores" %(assembler, cores))
    for subfolder, dirs, files in os.walk(input_folder):
        subfolder_path_elements = re.split("%s/" %input_folder, subfolder)
        if subfolder_path_elements[-1] != input_folder:
            sample_folder = subfolder_path_elements[-1]
            sample_id = re.split("_clean", sample_folder)[0]
            # Loop through each sample-folder and find read-files
            sample_output_folder = "%s/%s" %(out_dir, sample_id)
            if not os.path.exists(sample_output_folder):
                os.makedirs(sample_output_folder)
            for misc1, misc2, fastq in os.walk(subfolder):
                forward = ""
                backward = ""
                single_f = ""
                single_b = ""
                for element in fastq:
                    if sample_id in element and element.endswith("READ1.fastq"):
                        forward = "%s/%s" %(subfolder,element)
                    if sample_id in element and element.endswith("READ2.fastq"):
                        backward = "%s/%s" %(subfolder,element)
                    if sample_id in element and element.endswith("READ1-single.fastq"):
                        single_f = "%s/%s" %(subfolder,element)
                    if sample_id in element and element.endswith("READ2-single.fastq"):
                        single_b = "%s/%s" %(subfolder,element)
                if forward != "" and backward != "":
                    print ('#' * 50)
                    print ("Processing sample %s" %sample_id)
                    if assembler == "trinity":
                        assembly_trinity(forward,backward,sample_output_folder,sample_id,cores,min_length)
                        get_stats(sample_output_folder,sample_id)
                        cleanup_trinity_assembly_folder(sample_output_folder,sample_id)
                        print ("#" * 50)
                        mv_cmd = "mv %s/Trinity.fasta %s/%s.fasta" %(sample_output_folder,out_folder,sample_id)
                        os.system(mv_cmd)
                    elif assembler == "abyss":
                        assembly_abyss(forward,backward,single_f,single_b,sample_output_folder,sample_id,kmer,cores,args)
                        files = glob.glob(os.path.join(home_dir,'*'))
                        links = [f for f in files if os.path.islink(f)]
                        for l in links:
                            if l.endswith("-contigs.fa"):
                                contig_file = os.path.realpath(l)
                                mv_contig = "mv %s %s/../../%s.fa" %(contig_file,sample_output_folder,sample_id)
                                os.system(mv_contig)
                        mv_cmd1 = "mv %s/%s* %s" %(home_dir,sample_id,sample_output_folder)
                        os.system(mv_cmd1)
                        mv_cmd2 = "mv %s/coverage.hist %s" %(home_dir,sample_output_folder)
                        os.system(mv_cmd2)
                        contig_count_df = get_stats_abyss(sample_output_folder,sample_id,sample_contig_count_dict)
                else:
                    print ("Error: Read-files for sample %s could not be found.Please check if fastq file names end with 'READ1.fastq' and 'READ2.fastq' respectively." %sample_id)
                    raise SystemExit
    previous_stats_df = pd.read_csv(os.path.join(input_folder,'sample_overview.txt'),sep='\t')
    counter = 0
    for index,row in previous_stats_df.iterrows():
        sample_name = str(row['sample'])
        if sample_name in list(contig_count_df['sample']):
            new_info = contig_count_df[contig_count_df['sample']==sample_name]['total_contig_count']
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
    new_stats_df = pd.DataFrame(data=new_values_previous,columns=new_index)
    new_stats_df.to_csv(os.path.join(out_folder,'sample_stats.txt'),sep="\t",index=False)


def assembly_trinity(forw,backw,output_folder,id_sample,cores,min_length):
    print ("De-novo assembly with Trinity of sample %s:" %id_sample)
    command = [
        "Trinity",
        "--seqType",
        "fq",
        "--left",
        forw,
        "--right",
        backw,
         "--CPU",
        str(cores),
        "--min_contig_length",
        str(min_length),
        #"--JM",
        #"20G",
        "--output",
        output_folder
    ]
    try:
        print ("Building contigs........")
        with open(os.path.join(output_folder, "%s_trinity_screen_out.txt" %id_sample), 'w') as log_err_file:
            p = subprocess.Popen(command, stdout=log_err_file)
            p.communicate()
        print ("%s assembled. Trinity-stats are printed into %s" %(id_sample, os.path.join(output_folder, "%s_trinity_screen_out.txt" %sample_id)))
    except:
        print ("Could not assemble %s" %id_sample)

def assembly_abyss(forw,backw,singlef,singleb,output_folder,id_sample,kmer,cores,args):
    print ("De-novo assembly with abyss of sample %s:" %id_sample)
    command = [
        "abyss-pe",
        "k={}".format(kmer),
        "j={}".format(cores),
        'name={}'.format(id_sample),
        'in={} {}'.format(forw,backw)
    ]
    if args.single_reads:
        command.append('se={} {}'.format(singlef,singleb))
    try:
        print ("Building contigs........")
        with open(os.path.join(output_folder, "%s_abyss_screen_out.txt" %id_sample), 'w') as log_err_file:
            p = subprocess.Popen(command, stdout=log_err_file)
            p.communicate()
        print ("%s assembled. Statistics are printed into %s" %(id_sample, os.path.join(output_folder, "%s_abyss_screen_out.txt" %id_sample)))
    except:
        print ("Could not assemble %s" %id_sample)

def get_stats(sample_output_folder,sample_id):
    print ("Extracting statistics for", sample_id)
    # Read counts
    read_count_cmd = subprocess.Popen(["cat", "%s/both.fa.read_count" %sample_output_folder], stdout=subprocess.PIPE)
    read_count = read_count_cmd.communicate()[0]
    # Assembled read counts
    assembled_reads_cmd = subprocess.Popen(["wc", "-l", "%s/chrysalis/readsToComponents.out.sort" %sample_output_folder], stdout=subprocess.PIPE)
    assembled_reads = assembled_reads_cmd.communicate()[0]
    assembled_reads_count, file = assembled_reads.split(" ")
    # Contig count
    unimportant = ""
    contig_count_cmd = subprocess.Popen(["tail", "-n", "1", "%s/chrysalis/readsToComponents.out.sort" %sample_output_folder], stdout=subprocess.PIPE)
    contig_count_pre = contig_count_cmd.communicate()[0]
    print (contig_count_pre)
    contig_count, header, percent, sequence = contig_count_pre.split("\t")
    with open(os.path.join(sample_output_folder, "%s_stats.txt" %sample_id), 'w') as stat_file:
        stat_file.write("Statistics for sample %s\n" %sample_id)
        stat_file.write("Read-count in trimmed fastq read-files : %s" %read_count)
        stat_file.write("Reads assembled into contigs : %s\n" %assembled_reads_count)
        stat_file.write("Assembled contigs : %s\n" %contig_count)

def get_stats_abyss(sample_output_folder,sample_id,sample_contig_count_dict):
    contig_count_cmd = subprocess.Popen(["tail", "-n", "2", "%s/%s.fa" %('/'.join(sample_output_folder.split('/')[:-2]),sample_id)], stdout=subprocess.PIPE)
    contig_count_pre = contig_count_cmd.communicate()[0]
    contig_count = contig_count_pre.split(' ')[0].replace('>','')
    sample_contig_count_dict.setdefault(sample_id,contig_count)
    stats_df=pd.DataFrame.from_dict(sample_contig_count_dict, orient='index').reset_index()
    stats_df.columns = ['sample', 'total_contig_count']
    print('#'*50)
    print(stats_df)
    return(stats_df)
    #contig_count, header, percent, sequence = contig_count_pre.split("\t")    

def cleanup_trinity_assembly_folder(sample_output_folder, sample_id):
# This function is copied (and slightly modified) from phyluce, written by Brant Faircloth
    print ("Removing unnecessary files from the Trinity folder for %s" %sample_id)
    files = glob.glob(os.path.join(sample_output_folder, '*'))
    # check the names to make sure we're not deleting something improperly
    names = [os.path.basename(f) for f in files]
    try:
        assert "Trinity.fasta" in names
        assert "%s_trinity_screen_out.txt" %sample_id in names
    except:
        raise IOError("Neither Trinity.fasta nor %s_trinity_screen_out.txt were found in output." %sample_id)
    for file in files:
        if not os.path.basename(file) in ("Trinity.fasta", "%s_trinity_screen_out.txt" %sample_id, "%s_stats.txt" %sample_id):
            if os.path.isfile(file) or os.path.islink(file):
                os.remove(file)
            elif os.path.isdir(file):
                shutil.rmtree(file)
