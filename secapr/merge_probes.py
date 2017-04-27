import os
import sys
import glob
import shutil
import argparse
import commands
import subprocess

from .utils import CompletePath


# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="merge bait sequences from a sequence capture file (fasta) into full sequences",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--probe_file',
		required=True,
		action=CompletePath,
		help='the probe/bait file in fasta format'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		help='The output file with the new fasta sequences'
	)
	return parser.parse_args()


def read_fasta(fasta):
	name, seq = None, []
	for line in fasta:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))


def longest_common_substring(s1, s2):
	m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
	longest, x_longest = 0, 0
	for x in xrange(1, 1 + len(s1)):
		for y in xrange(1, 1 + len(s2)):
			if s1[x - 1] == s2[y - 1]:
				m[x][y] = m[x - 1][y - 1] + 1
				if m[x][y] > longest:
					longest = m[x][y]
					x_longest = x
			else:
				m[x][y] = 0
	return s1[x_longest - longest: x_longest]


args = get_args()
fasta_file = args.probe_file
out_file = args.output

header_list = []
sequence_list = []

with open(fasta_file) as f:
	for name, seq in read_fasta(f):
		header_list.append(name)
		sequence_list.append(seq)


sequence_dictionary = {}
for seq in sequence_list:
	place_in_list = sequence_list.index(seq)
	header = header_list[place_in_list]
	if place_in_list > 0:
		match = longest_common_substring(sequence_list[place_in_list],sequence_list[place_in_list-1])
		if not sequence_list[place_in_list-1].endswith(match):
			# sometimes there are cases where the matching substring is not at the end of the previous sequence,
			# since there are some slightly diffeent bases in the sequence of the following probe, even though the
			# sequences are clearly homologous. In that case we want to remove everything following the match from 
			# the previous probe and add this complete probe sequence into the final dictionary
			if len(match) > 20:
				splitstring = sequence_list[place_in_list-1].rsplit(match, 1)
				substring = match+splitstring[1]
				sequence_dictionary[new_header][-1] = sequence_dictionary[new_header][-1].replace(substring,'')
				sequence_dictionary[new_header].append(seq)
			else:
				new_header = header
				sequence_dictionary.setdefault(new_header,[])
				sequence_dictionary[new_header].append(seq)
				
		elif sequence_list[place_in_list-1].endswith(match) and sequence_list[place_in_list].startswith(match):
			if len(match)> 9:
				seq = seq.replace(match,'')
				sequence_dictionary[new_header].append(seq)
			else:
				print header, seq,  match, '\n'
				sequence_dictionary.setdefault(new_header,[])
				sequence_dictionary[new_header].append(seq)
		else:
			print 'unidentifiable probe:', header, seq
	elif place_in_list == 0:
		new_header = header
		sequence_dictionary.setdefault(new_header,[])
		sequence_dictionary[new_header].append(seq)



out_fasta = open(out_file, 'w')
for header in sorted(sequence_dictionary):
	sequence = ''.join(sequence_dictionary[header])
	out_fasta.write(header+"\n")
	out_fasta.write(sequence+"\n")
out_fasta.close()
