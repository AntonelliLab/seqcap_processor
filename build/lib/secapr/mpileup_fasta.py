import os
import argparse

from .utils import CompletePath


# Get arguments
def get_args():
    parser = argparse.ArgumentParser(description="Convert a mpileup file generated with samtools into fasta format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input',required=True,action=CompletePath,default=None,help='The mpileup file that you want to transform into fasta format')

    parser.add_argument('--fasta_out',required=True,action=CompletePath,default=None,help='The name of the output fasta file')

    parser.add_argument('--coverage',type=int,default=3,help='Minimum coverage-support that is required for making a base call. Anything below this threshold will be masked as ambiguity.')

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

with open(out_fasta, "wb") as f:
	for k, v in list(concat_basecalls.items()):
		f.write(">" + k+ "\n")
		f.write(v+ "\n")
