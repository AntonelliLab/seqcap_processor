# encoding: utf-8
'''
Extract the n loci with the best read-coverage from you reference-based assembly (bam-files)
'''
import os
import glob
import argparse
import subprocess
import csv
from .utils import CompletePath


# Get arguments
def add_arguments(parser):
	parser.add_argument(
		'--input',
		required=True,
		action=CompletePath,
		default=None,
		help='The folder with the results of the reference based assembly or the phasing results.'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be safed.'
	)
	parser.add_argument(
		'--n',
		type=int,
		default=30,
		help='The n loci that are best represented accross all samples will be extracted.'
	)


# find all subfolders in input folder
# if subfolder ends with "_remapped" or "_phased" find the bam-file
# run the function on the bam file
# make output dict with subfolder for every sample (_locus_selection)
# print a text file with an overview over the read-depth accross all loci
# choose the top n loci that have best coverage accross all samples (sum of scores accross all loci)
# store selected loci in separate bam file (if possible)
# store the consensus sequences for selected loci in sample_subfolder
# create a final_fasta_file with the selected sequences form all samples


def get_bam_path_dict(input_dir):
    type_input = ''
    subdirs = os.listdir(input_dir)
    sample_bam_dict = {}
    for subd in subdirs:
        if subd.endswith('_remapped'):
            type_input = 'unphased'
            sample_id = subd.split('_')[0]
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
            sample_id = subd.split('_')[0]
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
    sample_base = sample_base.split("_")[0]
    sample_base = sample_base.split(".")[0]
    print ('Reading read-depth info for %s.........' %sample_base)
    sample_dir = os.path.join(output_folder,'%s_locus_selection' %sample_base)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)


    get_read_depth = ["samtools", "depth", bam]
    read_depth_file = os.path.join(sample_dir,"%s_read_depth_per_position.txt" %sample_base)

    with open(read_depth_file, 'w') as logfile:
        sp1 = subprocess.Popen(get_read_depth, shell=False, stderr = subprocess.STDOUT, stdout=logfile)
        sp1.wait()
    return sample_dir,read_depth_file


def get_complete_loci_list(subfolder_file_dict):
    print('Generating locus database.........')
    locus_list = []
    for subfolder in subfolder_file_dict:
        read_depth_file = subfolder_file_dict[subfolder]
        #sample_id = read_depth_file.split('/')[-1].split('_')[0]
        with open(read_depth_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            reader = list(reader)
            for row in reader:
                locus_name = row[0]
                if locus_name not in locus_list:
                    locus_list.append(locus_name)
    return sorted(locus_list)



def summarize_read_depth_files(subfolder,read_depth_file,complete_locus_list,locus_dict_all_samples,sample_list):
    # 1. get a list of all read_depth files
    # iterate through list and calculate the average per locus
    # make sure loci are in same order for all samples
    # store results in list for each sample and then use zip() to combine as rows
    # write results into joint csv file for all samples
    sample_id = read_depth_file.split('/')[-1].split('_')[0]
    sample_list.append(sample_id)
    print ('Calculating coverage for all loci from bam files for %s.........' %sample_id)

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
            avg_read_depth = sum(sample_loci_dict[locus_name])/len(sample_loci_dict[locus_name])
            locus_dict_all_samples.setdefault(locus_name,[])
            locus_dict_all_samples[locus_name].append(avg_read_depth)
        else:
            avg_read_depth = 0.0
            locus_dict_all_samples.setdefault(locus_name,[])
            locus_dict_all_samples[locus_name].append(avg_read_depth)

    return locus_dict_all_samples 


def main(args):
    input_dir = args.input
    output_folder = args.output
    if not os.path.exists(output_folder):
	    os.makedirs(output_folder)
    else:
	    raise IOError("The directory {} already exists.  Please check and remove by hand.".format(output_folder))


    # Create a dictionary containing the bam-file paths for each sample and tell if data is phased or unphased
    sample_bam_dict, input_type = get_bam_path_dict(input_dir)
    if input_type == 'unphased':
        subfolder_list = []
        subfolder_file_dict = {}
        for key in sample_bam_dict:
            if key.endswith('_remapped'):
                bam = sample_bam_dict[key][0]
                sample_dir, read_depth_file = get_bam_read_cov(bam,output_folder)
                subfolder_file_dict.setdefault(sample_dir,read_depth_file)

        locus_list = get_complete_loci_list(subfolder_file_dict)
        locus_dict_all_samples = {}
        sample_id_list = []
        for subfolder in subfolder_file_dict:
            read_depth_file = subfolder_file_dict[subfolder]
            locus_dict_all_samples = summarize_read_depth_files(subfolder,read_depth_file,locus_list,locus_dict_all_samples,sample_id_list)
        
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
        for column in output_dict:
            data = output_dict[column]
            final_data_list.append(data)
        
        output = open("%s/average_cov_per_locus.csv" %output_folder, "w")
        outlog=csv.writer(output, delimiter='\t')
        transformed_data = zip(*final_data_list)
        for row in transformed_data:
            outlog.writerow(row)










	

#
#bam = args.bam
#bam_name = bam.split("/")[-1]
#sample_base = bam_name.split(".bam")[0]
#sample_base = sample_base.split("_")[0]
#sample_base = sample_base.split(".")[0]
#output_folder = args.output
#if not os.path.exists(output_folder):
#	os.makedirs(output_folder)
#sample_dir = os.path.join(output_folder,sample_base)
#if not os.path.exists(sample_dir):
#	os.makedirs(sample_dir)
#
## find the samtools path
#samtools = os.popen("which %s" % "samtools").read().strip()
#
#get_read_depth = [samtools, "depth", bam]
#read_depth_file = os.path.join(sample_dir,"%s_read_depth_per_position.txt" %sample_base)
#
#with open(read_depth_file, 'w') as logfile:
#	sp1 = subprocess.Popen(get_read_depth, shell=False, stderr = subprocess.STDOUT, stdout=logfile)
#	sp1.wait()
#
#loci_dict = {}
#with open(read_depth_file, 'r') as f:
#	reader = csv.reader(f, delimiter='\t')
#	reader = list(reader)
#	for row in reader:
#		locus_name = row[0]
#		position = row[1]
#		coverage = int(row[2])
#		loci_dict.setdefault(locus_name,[])
#		loci_dict[locus_name].append(coverage)
#
#output = open("%s/%s_average_cov_per_locus.csv" %(sample_dir,sample_base), "wb")
#outlog=csv.writer(output, delimiter='\t')
#outlog.writerow(["locus_name",sample_base])
#for locus in loci_dict:
#	locus_list = locus.split("_")[:3]
#	locus_name = "_".join(locus_list)
#	avg_read_depth = sum(loci_dict[locus])/len(loci_dict[locus])
#	outlog.writerow([locus_name, avg_read_depth])

