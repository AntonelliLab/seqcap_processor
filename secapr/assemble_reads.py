#author: Tobias Andermann, tobias.andermann@bioenv.gu.se

'''
Assemble trimmed Illumina read files (fastq)
'''

import os
import sys
import re
import glob
import shutil
import argparse
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
import time

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
        help='Call the folder that contains the trimmed reads, organized in a separate subfolder for each sample. The name of the subfolder has to start with the sample name, delimited with an underscore [_] (default output of secapr clean_reads function)'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CompletePath,
        default=None,
        help='The output directory where results will be saved'
    )
    parser.add_argument(
        '--assembler',
        choices=["spades","abyss","trinity"],
        default="spades",
        help="""The assembler to use (default = spades)."""
    )
    parser.add_argument(
        '--kmer',
        type=str,
        help='Set the kmer value (only available for Abyss and Spades). Provide single value for Abyss, or list of kmers for Spades, e.g. "--kmer 21,33,55". Default for Abyss is 35, and for spades it is 21,33,55,77,99,127. Note that Spades only accepts uneven kmer values.'
    )
    parser.add_argument(
        '--contig_length',
        type=int,
        default=200,
        help='Set the minimum contig length for the assembly. Contigs that are shorter than this threshold will be discarded.'
    )
    parser.add_argument(
        '--max_memory',
        type=str,
        help='Set the maximum memory to be used during assembly in GB (only available for Spades and Trinity). This can be necessary when working with computing nodes with limited memory or to avoid over-allocation of computing resources on clusters which can in some cases cause your assembly to be stopped or interrupted.'
    )
    parser.add_argument(
        '--single_reads',
        action='store_true',
        default=False,
        help='Use this flag if you additionally want to use single reads for the assembly'
    )
    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help='For parallel processing you can set the number of cores you want to run the assembly on.'
    )


def assembly_trinity(forw,backw,output_folder,id_sample,cores,min_length,max_memory):
    print(("De-novo assembly with Trinity of sample %s:" %id_sample))
    #print(output_folder)
    if not max_memory:
        max_memory = '8G'
    else:
        max_memory = '%sG'%max_memory
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
        #"8G",
        "--max_memory",
        max_memory, 
        #"--bypass_java_version_check",
        #"--normalize_reads",
        "--output",
        output_folder
    ]
    print ("Building contigs........")
    with open(os.path.join(output_folder, "%s_trinity_screen_out.txt" %id_sample), 'w') as log_err_file:
        p = subprocess.Popen(command, stdout=log_err_file, stderr=log_err_file)
        p.communicate()
    filename = os.path.join(output_folder, "%s_trinity_screen_out.txt" %id_sample)
    file_object  = open(filename, 'r')
    for line in file_object:
        if line.startswith('Error'):
            print(line)
            print ('SECAPR NOTE:\nTrinity is currently only functional in the Linux distribution of SECAPR due to Java incompatibilities.\n')
                #'However, the environment on MacOS machines can be easily altered by hand in order to properly run Trinity.\n',
                #'This might however compromise the functionality of other parts of the SECAPR pipeline, therefore we recommend to undo the changes made in the envrionment after using Trinity by following the instructions below.\n\n',
                #'In order to run the Trinity assembly on MacOS do the following:\n',
                #'1. within the SECAPR conda envrionment type: "conda install openjdk=7"\n',
                #'2. run the secapr assemble_reads function with Trinity (using the "--assembler trinity" flag)\n',
                #'3. after assembly rebuild the SECAPR default environment by typing "conda install trimmomatic=0.33"\n'
            sys.exit()
        elif line.startswith('Trinity run failed.'):
            print(filename)
            print ('SECAPR NOTE:\nTrinity is currently only functional in the Linux distribution of SECAPR.\n')
            sys.exit()

    print(("%s assembled. Trinity-stats are printed into %s" %(id_sample, os.path.join(output_folder, "%s_trinity_screen_out.txt" %id_sample))))
    #except:
    #    print ("Trinity failed, maybe due to limited stack-size. Try increase stacksize with command 'zsh | ulimit -s unlimited | sh' and run again.")

def assembly_abyss(forw,backw,singlef,singleb,output_folder,id_sample,kmer,cores,args):
    print("WARNING: Abyss is very memory heavy and depending on the size of your read files may throw an error because it's running out of memory. If running on a cluster, ask your system administrator how to allocate more memory to your abyss job.")
    if cores > 1:
        print('WARNING: You chose to run Abyss on more than 1 core. This can cause problems on some systems and will make the script crash. In that case try running Abyss on a sinlge core instead.')
    print(("De-novo assembly with abyss of sample %s:" %id_sample))
    try:
        kmer = int(kmer)
    except:
        quit('\n\nError: Provided kmer value could not be formatted as integer. Please provide single numeric kmer value when choosing the Abyss assembler.')
    command = [
        "abyss-pe",
        "--directory={}".format(output_folder),
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
            p.wait()
        print(("%s assembled. Statistics are printed into %s" %(id_sample, os.path.join(output_folder, "%s_abyss_screen_out.txt" %id_sample))))
    except:
        print(("Could not assemble %s" %id_sample))

def assembly_spades(forw,backw,singlef,singleb,output_folder,id_sample,kmer,cores,max_memory,args):
    print(("De-novo assembly with spades of sample %s:" %id_sample))
    kmer = str(kmer)        
    command = [
        "spades.py",
        "-k",
        kmer,
        "--only-assembler",
        "--pe1-1",
        forw,
        "--pe1-2",
        backw,
        "-o",
        output_folder
    ]
    if args.single_reads:
        command+=["--pe1-s", singlef, "--pe1-s",singleb]
    if args.cores > 1:
        command+=["--threads", str(args.cores)]
    if args.max_memory:
        command+=["--memory", str(args.max_memory)]
    # try:
    print ("Building contigs........")
    with open(os.path.join(output_folder, "%s_spades_screen_out.txt" %id_sample), 'w') as log_err_file:
        p = subprocess.Popen(command, stdout=log_err_file)
        p.communicate()
        p.wait()
    print(("%s assembled. Statistics are printed into %s" %(id_sample, os.path.join(output_folder, "%s_spades_screen_out.txt" %id_sample))))
    # except:
    #     print(("Could not assemble %s" %id_sample))

def get_trinity_stats(sample_output_folder,sample_id,sample_contig_count_dict):
    print(("Extracting statistics for %s" %str(sample_id)))
    contig_file = "%s/Trinity.fasta" %sample_output_folder
    new_contig_file = "%s/Trinity_formatted.fasta" %sample_output_folder
    edit_trinity_headers(contig_file,new_contig_file)
    contig_count = count_contigs(new_contig_file)
    sample_contig_count_dict.setdefault(sample_id,contig_count)
    stats_df=pd.DataFrame.from_dict(sample_contig_count_dict, orient='index').reset_index()
    stats_df.columns = ['sample', 'total_contig_count']
    print(('#'*50))
    print(stats_df)
    return(stats_df)

def edit_trinity_headers(contig_file,new_contig_file):
    fasta =  open(contig_file,'r')
    new_fasta = open(new_contig_file,'w')
    counter = 0
    for line in fasta:
        if line.startswith('>'):
            readcount = int(re.sub(r'.*len=([0-9]*).*','\\1',line).strip())
            new_header = '>%i %i XXX\n' %(counter,readcount)
            counter += 1
            new_fasta.write(new_header)
        else:
            new_fasta.write(line)
    new_fasta.close()

def count_contigs(contig_file):
    """Return a count of contigs from a fasta file"""
    return sum([1 for line in open(contig_file, 'r').readlines() if line.startswith('>')])

def remove_short_contigs(contig_file,min_length):
    fasta_content = list(SeqIO.parse(open(contig_file),'fasta'))
    new_fasta_content = []
    for record in fasta_content:
        contig_length = len(str(record.seq))
        if contig_length < min_length:
            pass
        else:
            new_fasta_content.append(record)
    SeqIO.write(new_fasta_content,contig_file, 'fasta-2line')
            

def get_stats_abyss(sample_output_folder,sample_id,sample_contig_count_dict):
    #contig_count_cmd = subprocess.Popen(["tail", "-n", "2", "%s/%s.fa" %('/'.join(sample_output_folder.split('/')[:-2]),sample_id)], stdout=subprocess.PIPE)
    #contig_count_pre = contig_count_cmd.communicate()[0]
    contig_file = "%s/%s.fa" %('/'.join(sample_output_folder.split('/')[:-2]),sample_id)
    contig_count = count_contigs(contig_file)
    #contig_count = contig_count_pre.split(' ')[0].replace('>','')
    sample_contig_count_dict.setdefault(sample_id,contig_count)
    stats_df=pd.DataFrame.from_dict(sample_contig_count_dict, orient='index').reset_index()
    stats_df.columns = ['sample', 'total_contig_count']
    print(('#'*50))
    print(stats_df)
    return(stats_df,contig_file)
    #contig_count, header, percent, sequence = contig_count_pre.split("\t")    

def get_stats_spades(contig_file,sample_id,sample_contig_count_dict):
    #contig_count_cmd = subprocess.Popen(["tail", "-n", "2", "%s/%s.fa" %('/'.join(sample_output_folder.split('/')[:-2]),sample_id)], stdout=subprocess.PIPE)
    #contig_count_pre = contig_count_cmd.communicate()[0]
    contig_count = count_contigs(contig_file)
    #contig_count = contig_count_pre.split(' ')[0].replace('>','')
    sample_contig_count_dict.setdefault(sample_id,contig_count)
    stats_df=pd.DataFrame.from_dict(sample_contig_count_dict, orient='index').reset_index()
    stats_df.columns = ['sample', 'total_contig_count']
    print(('#'*50))
    print(stats_df)
    return(stats_df,contig_file)
    #contig_count, header, percent, sequence = contig_count_pre.split("\t")    

def cleanup_trinity_assembly_folder(sample_output_folder, sample_id):
# This function is copied (and slightly modified) from phyluce, written by Brant Faircloth
    print(("Removing unnecessary files from the Trinity folder for %s" %sample_id))
    files = glob.glob(os.path.join(sample_output_folder, '*'))
    # check the names to make sure we're not deleting something improperly
    names = [os.path.basename(f) for f in files]
    try:
        assert "Trinity.fasta" in names
        assert "%s_trinity_screen_out.txt" %sample_id in names
    except:
        raise IOError("Neither Trinity.fasta nor %s_trinity_screen_out.txt were found in output." %sample_id)
    for file in files:
        if not os.path.basename(file) in ("Trinity.fasta","Trinity_formatted.fasta", "%s_trinity_screen_out.txt" %sample_id, "%s_stats.txt" %sample_id):
            if os.path.isfile(file) or os.path.islink(file):
                os.remove(file)
            elif os.path.isdir(file):
                shutil.rmtree(file)


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
    assembler = args.assembler

    if args.kmer:
        kmer = str(args.kmer)
    else:
        if assembler == 'spades':
            kmer = '21,33,55,77,99,127'
        else:
            kmer = 35
    if args.max_memory:
        max_memory = args.max_memory
    else:
        max_memory = None
    #home_dir = os.getcwd()
    sample_contig_count_dict = {}
    if cores > 1:
        print(("Running %s parallel on %d cores" %(assembler, cores)))
    for subfolder, dirs, files in os.walk(input_folder):
        subfolder_path_elements = re.split("%s/" %input_folder, subfolder)
        if subfolder_path_elements[-1] != input_folder:
            sample_folder = subfolder_path_elements[-1]
            sample_id = re.split("_clean", sample_folder)[0]
            # Loop through each sample-folder and find read-files
            sample_output_folder = "%s/%s" %(out_dir, sample_id)
            if assembler == "trinity":
                sample_output_folder = '%s_trinity'%sample_output_folder            
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
                    print(('#' * 50))
                    print(("Processing sample %s" %sample_id))
                    start = time.time()
                    if assembler == "trinity":
                        assembly_trinity(forward,backward,sample_output_folder,sample_id,cores,min_length,max_memory)
                        contig_count_df = get_trinity_stats(sample_output_folder,sample_id,sample_contig_count_dict)
                        cleanup_trinity_assembly_folder(sample_output_folder,sample_id)
                        print(("#" * 50))
                        mv_cmd = "mv %s/Trinity_formatted.fasta %s/%s.fasta" %(sample_output_folder,out_folder,sample_id)
                        os.system(mv_cmd)
                    elif assembler == "abyss":
                        assembly_abyss(forward,backward,single_f,single_b,sample_output_folder,sample_id,kmer,cores,args)
                        files = glob.glob(os.path.join(sample_output_folder,'*'))
                        links = [f for f in files if os.path.islink(f)]
                        for l in links:
                            if l.endswith("-contigs.fa"):
                                contig_file = os.path.realpath(l)
                                mv_contig = "mv %s %s/../../%s.fa" %(contig_file,sample_output_folder,sample_id)
                                os.system(mv_contig)
                        #mv_cmd1 = "mv %s/%s* %s" %(home_dir,sample_id,sample_output_folder)
                        #os.system(mv_cmd1)
                        #mv_cmd2 = "mv %s/coverage.hist %s" %(home_dir,sample_output_folder)
                        #os.system(mv_cmd2)
                        contig_count_df,contig_file = get_stats_abyss(sample_output_folder,sample_id,sample_contig_count_dict)
                        remove_short_contigs(contig_file,min_length)
                    elif assembler == 'spades':
                        assembly_spades(forward,backward,single_f,single_b,sample_output_folder,sample_id,kmer,cores,max_memory,args)
                        contig_file = os.path.join(sample_output_folder,'contigs.fasta')
                        new_contig_file = '%s/../../%s.fa'%(sample_output_folder,sample_id)
                        mv_contig = "cp %s %s" %(contig_file,new_contig_file)
                        os.system(mv_contig)
                        contig_count_df,contig_file = get_stats_spades(new_contig_file,sample_id,sample_contig_count_dict)
                        remove_short_contigs(new_contig_file,min_length)
                    end = time.time()
                    print('Assembled contigs for sample %s in %i minutes' %(sample_id,int(np.round((end-start)/60))))
                else:
                    print(("Error: Read-files for sample %s could not be found.Please check if fastq file names end with 'READ1.fastq' and 'READ2.fastq' respectively and if all files are unzipped." %sample_id))
                    raise SystemExit

    try:
        previous_stats_df = pd.read_csv(os.path.join(input_folder,'sample_stats.txt'),sep='\t')
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
                    new_values_previous = np.vstack([new_values_previous, new_values])
                counter += 1
        new_stats_df = pd.DataFrame(data=new_values_previous,columns=new_index)

    except:
        print('No previous stats file found, creating new stats file.')
        new_stats_df = contig_count_df

    new_stats_df.to_csv(os.path.join(out_folder,'sample_stats.txt'),sep="\t",index=False)



