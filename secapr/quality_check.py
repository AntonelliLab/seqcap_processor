# encoding: utf-8
#author: Tobias Andermann, tobias.andermann@bioenv.gu.se
"""
This script runs a fastqc test on all fastq samples in a user-provided folder and creates an overview plot,
"""

import os
import sys
import glob
import fnmatch
import shutil
import ConfigParser
from .utils import CompletePath
import subprocess
import pandas as pd
import numpy as np

def add_arguments(parser):
    parser.add_argument(
        '--input',
        required=True,
        action=CompletePath,
        default=None,
        help='The directory containing fastq files'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CompletePath,
        default=None,
        help='The output directory where quality-test results will be saved'
    )

def main(args):
    # Set working directory
    out_folder = args.output
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # Get list of all fastq-files
    input_folder = args.input
    matches = []
    for root, dirnames, filenames in os.walk(input_folder):
        for filename in fnmatch.filter(filenames, '*.fastq'):
            matches.append(os.path.join(root, filename))
    fastq_df = pd.DataFrame(index=np.arange(0,len(matches)), columns=['filepaths'])
    fastq_df['filepaths'] = matches
    fastq_list_path = os.path.join(out_folder,'fastq_file_list.txt')
    fastq_df.to_csv(fastq_list_path,index=False,header=False,sep='\t')

    # run FASTQC
    fastqc_cmd = [
        'fastqc -o %s -f fastq $(cat %s)' %(out_folder,fastq_list_path)
    ]
    with open(os.path.join(out_folder, "fastqc_screen_out.txt"), 'w') as log_err_file:
        p = subprocess.Popen(fastqc_cmd, stdout=log_err_file, stderr=log_err_file, shell=True)
        p.communicate()

    # write the r-plotting script to file
    r_plotting_script = 'opt <- c()\nopt$input_folder = workdir\nopt$output_file =paste0(workdir, "/QC_plots.pdf")\n\n#load fastQC summaries and create per test table\ninp <- list.files(opt$input_folder, pattern = ".zip")\n\n\nfastqc_results <- lapply(inp, function(k){\n  unzip(paste(opt$input_folder, k, sep = "/"),exdir = opt$input_folder)\n  inpu <- read_delim(paste(paste(gsub(".zip", "", paste(opt$input_folder,k, sep = "/"))), \n                           "summary.txt", sep = "/"), delim = "\t")\n  out <- as_data_frame(t(inpu[, 1])) %>%\n    mutate(sample.id = names(inpu)[3])\n  names(out) <- c(gsub(" ", "_", unlist(inpu[,2])), "sample_id")\n  unlink(x = paste(opt$input_folder, gsub(".zip", "", k), sep = "/"), \n         recursive = T, force = T)\n  \n  return(out)\n})\n\nn_cols = length(fastqc_results[[1]])\nif (n_cols==11){\n  outp <- do.call("rbind.data.frame", fastqc_results)%>%\n    select(ID = sample_id,\n           PBQ = Per_base_sequence_quality,\n           PTQ = Per_tile_sequence_quality,\n           PSQ = Per_sequence_quality_scores,\n           PBC = Per_base_sequence_content,\n           SGC = Per_sequence_GC_content,\n           PBN = Per_base_N_content,\n           SLD = Sequence_Length_Distribution,\n           SDL = Sequence_Duplication_Levels,\n           ORS = Overrepresented_sequences,\n           AdC = Adapter_Content)\n  #change table format\n  ret <- outp %>% \n    group_by(ID) %>%\n    gather(test, status, PBQ:AdC)\n\n}\nif (n_cols==12){\n  outp <- do.call("rbind.data.frame", fastqc_results)%>%\n    select(ID = sample_id,\n           PBQ = Per_base_sequence_quality,\n           PTQ = Per_tile_sequence_quality,\n           PSQ = Per_sequence_quality_scores,\n           PBC = Per_base_sequence_content,\n           SGC = Per_sequence_GC_content,\n           PBN = Per_base_N_content,\n           SLD = Sequence_Length_Distribution,\n           SDL = Sequence_Duplication_Levels,\n           ORS = Overrepresented_sequences,\n           AdC = Adapter_Content,\n           KmC = Kmer_Content)\n  #change table format\n  ret <- outp %>% \n    group_by(ID) %>%\n    gather(test, status, PBQ:KmC)\n}\n\n#plot how many samples failed the test\nqc.fail <- ggplot()+\n  geom_bar(data = ret, aes(x = test, fill = status), stat = "count", position = "dodge")+\n  theme_bw()\n\n#plot which sample failed which test\nqc.samples <- ggplot()+\n  geom_tile(data = ret, aes(y = ID, x = test, fill = as.factor(status)))+\n  scale_fill_discrete(name = "status")+\n  xlab("FastQC test")+\n  ylab("Samples")+\n  theme_bw()+\n  theme(\n    axis.text.y = element_blank()\n  )\n\n#plot pdf\npdf(opt$output_file)\nprint(qc.fail)\nprint(qc.samples)\ndev.off()\n\npng(gsub(".pdf", "1.png", opt$output_file))\nprint(qc.fail)\ndev.off()\n\npng(gsub(".pdf", "2.png", opt$output_file))\nprint(qc.samples)\ndev.off()\n\n#table with samples that faild a test\nfail <- ret %>%\n  filter(status == "FAIL")\n\n#get the ID number of the failed samples\nfail.samp <- fail %>%\n  filter(!duplicated(ID)) %>%\n  select(ID)%>%\n  unlist() %>%\n  parse_number()%>%\n  unique() %>%\n  sort()'
    
    add_to_script = 'library(tidyverse)\nworkdir = "%s"\n' %out_folder 
    new_r_plotting_script = add_to_script + r_plotting_script
    r_script_path = os.path.join(out_folder,'fastqc_visualization.r')
    text_file = open(r_script_path, "w")
    text_file.write(new_r_plotting_script)
    text_file.close()
    
    # execute r-plotting script
    print('Running R-code for plotting...')
    final_plot = os.path.join(out_folder,'quality_summary_all_samples.pdf')
    plotting_cmd = [
        'Rscript %s -i %s -o %s' %(r_script_path,out_folder,final_plot)
    ]
    with open(os.path.join(out_folder, "r_plotting_screen_out.txt"), 'w') as log_err_file:
        p = subprocess.Popen(plotting_cmd, stdout=log_err_file, stderr=log_err_file, shell=True)
        p.communicate()