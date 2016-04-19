#!/usr/bin/python2.7
import os
import argparse


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Input %%%


# Complete path function
class CompletePath(argparse.Action):
    """give the full path of an input file/folder"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


# Get arguments
def get_args():
    parser = argparse.ArgumentParser(description="Convert a mpileup file generated with samtools into fasta format",		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input',required=True,action=CompletePath,default=None,help='The mpileup file that you want to transform into fasta format')

    parser.add_argument('--fasta_out',required=True,action=CompletePath,default=None,help='The name of the output fasta file')

    parser.add_argument('--coverage',type=int,default=4,help='Minimum coverage that is required for building consensus. Anything below this threshold will be skipped.')

    return parser.parse_args()

# Preparation for calling input variables and files
args = get_args()


pileup = args.input
out_fasta = args.fasta_out
cov = args.coverage
with open(pileup) as f:
    content = [x.strip('\n') for x in f.readlines()]

loci = []
seq_dict = {}
for line in content:
    if '#' not in line:
        element = line.split('\t')
        # Create a list of all loci names
        seq_name = element[0]
        if seq_name not in loci:
            loci.append(seq_name)
        # By default call every position a uncertainty
        basecall = "N"
        # Only if sufficient read depth is given, make an informaed basecall
        if int(element[3]) >= cov:
            reference = element[2]
            #print "Reference is", reference
            basecall = reference

            match_ref = "." ","
            bases = "A" "G" "C" "T"
            sample = element[4].upper()
            #print sample
            # make a directory with all different basecalls and count their occurences
            calls = dict((letter,sample.count(letter)) for letter in set(sample))
            #print calls
            # find out how many agreements with reference are among the basecalls. These are all . and , basecalls.
            # reset the counter before every round (every position)
            list_matches = 0
            for key,value in calls.items():
                if key in match_ref:
                    list_matches += value

            for key,value in calls.items():
                if key in bases:
                    if value >= list_matches:
                        basecall = key
                        #print sample

        #print "Basecall is", basecall, "\n"


        seq_dict.setdefault(element[0],[])
        seq_dict[element[0]].append(basecall)

# Join all basecalls for each key (=locus) into one sequence and deposit in new dictionary
concat_basecalls = {}
for key, value in seq_dict.items():
	concat_basecalls[key] = "".join(value)
#print concat_basecalls

with open(out_fasta, "wb") as f:
	for k, v in concat_basecalls.items():
		f.write(">" + k+ "\n")
		f.write(v+ "\n")
