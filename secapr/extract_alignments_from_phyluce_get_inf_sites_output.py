#!/usr/bin/python2.7
import os
import argparse
import csv
import random

# Complete path function
class CompletePath(argparse.Action):
	"""give the full path of an input file/folder"""
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

# Get arguments
def get_args():
	parser = argparse.ArgumentParser(description="Use the phyluce_align_get_informative_sites output to extract UCE alignments with a certain number of informative sites (or random ones)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('--input',required=True,action=CompletePath,default=None,help='The phyluce_align_get_informative_sites screen output as text file')

	parser.add_argument('--output',required=True,action=CompletePath,default=None,help='The name of the file, which will eb a list of taxa fulfilling the requirement')
	
	parser.add_argument('--mode',choices=["top","bottom","cutoff","random"],default="top",help='Choose which alignments you want to extract: top = the x most informative alignments   bottom = the x least informative alignments   cutoff = all alignments with more than x informative sites   random = randomly chooses x alignments   x is specified with the --threshold flag')
	
	parser.add_argument('--threshold',type=int,default=15,help='Minimum coverage-support that is required for making a base call. Anything below this threshold will be masked as ambiguity.')

	return parser.parse_args()

# Preparation for calling input variables and files
args = get_args()


input = args.input
output = args.output
out_file = output.split("/")[-1]
out_dir = '/'.join(output.split("/")[:-1])
mode = args.mode
threshold = args.threshold

def getkey(string):
	locus, misc, values = string.partition('\t')
	x, y, values = values.partition('\t')
	insites, miscx, miscy = values.partition('\t')
	return int(insites)


output_file = open("%s/%s_%s_%s" %(out_dir,mode,threshold,out_file), "wb")
uce_list=csv.writer(output_file)

with open(input) as f:
	content = [x.strip('\n') for x in f.readlines()]
	header = content[0]
	tail = content[-0]
	body = content[1:-1]
	body = sorted(body, key=getkey)
	if mode == 'cutoff':
		for line in body:
			element = line.split('\t')
			if int(element[2]) >= int(threshold):
				print line
				uce_list.writerow([element[0]])
	elif mode == 'top':
		for line in body[-threshold:]:
			element = line.split('\t')
			print line
			uce_list.writerow([element[0]])
	elif mode == 'bottom':
		for line in body[:threshold]:
			element = line.split('\t')
			print line
			uce_list.writerow([element[0]])
	elif mode =='random':
		random_body = random.sample(body, threshold)
		for line in random_body:
			element = line.split('\t')
			print line
			uce_list.writerow([element[0]])