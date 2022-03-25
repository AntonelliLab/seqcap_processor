'''Join exon-alignment files belonging to the same gene
'''

#author: Tobias Andermann, tobias.andermann@bioenv.gu.se

import os
import re
#from cogent import LoadSeqs, DNA

from .utils import CompletePath


def add_arguments(parser):
	parser.add_argument(
		'--input',
		required=True,
		action=CompletePath,
		default=None,
		help='The directory containing the fasta-alignment-files'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be saved'
	)

	
def read_fasta(fasta):
	name, seq = None, []
	for line in fasta:
		line = line.rstrip()
		if line.startswith(">"):
			if name:
				yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name:
		yield (name, ''.join(seq))


def main(args):
	work_dir = args.input
	out_dir = args.output
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	# Create a dictionary with the name-pattern as key and all file-names sharing that name-pattern
	fasta_dict = {}
	for fasta in os.listdir(work_dir):
		if fasta.endswith(".fasta") or fasta.endswith(".fa"):
			fasta_split = re.split("_", fasta)
			name_pattern = "%s_%s" %(fasta_split[0], fasta_split[1])
			fasta_dict.setdefault(name_pattern,[]).append(fasta)

		else:
			print("didn't work for", fasta)

	# Get the list of taxa names (headers) for each locus, key is out-file, values are in-files
	for key, value in fasta_dict.items():
		print(key)
		list_headers = []
		# Each k is a separate fasta input file belonging to the same locus (to be joined)
		for k in sorted(value):
			with open("%s/%s" %(work_dir, k)) as f:
				for name, seq in read_fasta(f):
					if name not in list_headers:
						list_headers.append(name)

		# Find the missing taxa in each fasta input file and simulate a sequence of correct length (only "?")
		in_fasta = os.path.join(work_dir, fasta)
		# Each k is a separate fasta input file belonging to the same locus (to be joined)
		all_seq_dict = {}
		for k in sorted(value):
			taxa_names_single = []
			present_seq = []
			length_alignment = ""
			with open("%s/%s" %(work_dir,k)) as f:
				for name, seq in read_fasta(f):
					taxa_names_single.append(name)
					present_seq.append((name,seq))
					length_alignment = len(seq)
			# Make a list of all missing taxa in each fasta input file
			missing_taxa = []
			for header in list_headers:
				if header not in taxa_names_single:
					missing_taxa.append(header)
			simulated_seq = []
			for mistax in missing_taxa:
				fake_string = "?" * length_alignment
				simulated_seq.append((mistax,fake_string))
			all_seq = sorted(simulated_seq+present_seq)

			for seq_header, sequence in all_seq:
				all_seq_dict.setdefault(seq_header,[]).append(sequence)

		with open(os.path.join(out_dir, "%s.fasta" %key), 'w') as out_fasta:
			for seqname, sequences in all_seq_dict.items():
				final_sequence = "".join(sequences)
				out_fasta.write(seqname+"\n")
				out_fasta.write(final_sequence+"\n")
