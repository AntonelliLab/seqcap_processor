import os
import argparse
import subprocess
import csv


# Complete path function
class CompletePath(argparse.Action):
	"""give the full path of an input file/folder"""
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="Get average read depth for each locus from a bam file",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--bam',
		required=True,
		action=CompletePath,
		default=None,
		help='Your !sorted! bam input file'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be safed.'
	)
	return parser.parse_args()

args = get_args()

bam = args.bam
bam_name = bam.split("/")[-1]
sample_base = bam_name.split(".bam")[0]
sample_base = sample_base.split("_")[0]
sample_base = sample_base.split(".")[0]
output_folder = args.output
if not os.path.exists(output_folder):
	os.makedirs(output_folder)
sample_dir = os.path.join(output_folder,sample_base)
if not os.path.exists(sample_dir):
	os.makedirs(sample_dir)

# find the samtools path
samtools = os.popen("which %s" % "samtools").read().strip()

get_read_depth = [samtools, "depth", bam]
read_depth_file = os.path.join(sample_dir,"%s_read_depth_per_position.txt" %sample_base)

with open(read_depth_file, 'w') as logfile:
	sp1 = subprocess.Popen(get_read_depth, shell=False, stderr = subprocess.STDOUT, stdout=logfile)
	sp1.wait()

loci_dict = {}
with open(read_depth_file, 'r') as f:
	reader = csv.reader(f, delimiter='\t')
	reader = list(reader)
	for row in reader:
		locus_name = row[0]
		position = row[1]
		coverage = int(row[2])
		loci_dict.setdefault(locus_name,[])
		loci_dict[locus_name].append(coverage)

output = open("%s/%s_average_cov_per_locus.csv" %(sample_dir,sample_base), "wb")
outlog=csv.writer(output, delimiter='\t')
outlog.writerow(["locus_name",sample_base])
for locus in loci_dict:
	locus_list = locus.split("_")[:3]
	locus_name = "_".join(locus_list)
	avg_read_depth = sum(loci_dict[locus])/len(loci_dict[locus])
	outlog.writerow([locus_name, avg_read_depth])

