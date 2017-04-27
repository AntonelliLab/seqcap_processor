# encoding: utf-8

"""
Copyright (c) 2010-2012, Brant C. Faircloth All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

* Neither the name of the University of California, Los Angeles nor the names
of its contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Modified by Tobias Hofmann (tobiashofmann@gmx.net):
# Modifications include: 	- removal of unnecessary flag-options
# 				- renaming of uce-related parameters into exon-terminology
#				- modification of regex-patterns to match Trinity-contigs and palm-exon set probe-order file
# Additions include:		- automatic generation of config file including all recovered exon names for further processing in this pipeline
# 				- automatic generation of the match-table in tab-delimited text-format to be opened in e.g. Excel for match overview
# 				- user-input choice for path to sqlite3, which is necesarry to access the database and generate the match-text-file


from __future__ import print_function
import re
import os
import sys
import glob
import copy
import operator
import itertools
import logging
import sqlite3
import argparse
from phyluce import lastz
from phyluce.helpers import is_dir, is_file, FullPaths
from collections import defaultdict
from Bio import SeqIO

#import pdb
log = logging.getLogger(__name__)

def add_arguments(parser):
        parser.add_argument(
		'--contigs',
		required=True,
		type=is_dir,
		action=FullPaths,
		help="The directory containing the assembled contigs in fasta format."
	)
	parser.add_argument(
		'--reference',
		required=True,
		type=is_file,
		action=FullPaths,
		help="The fasta-file containign the reference sequences (probe-order-file)."
	)
	parser.add_argument(
		'--output',
		required=True,
		action=FullPaths,
		help="The directory in which to store the resulting SQL database and LASTZ files."
	)
	parser.add_argument(
		"--verbosity",
		type=str,
		choices=["INFO", "WARN", "CRITICAL"],
		default="INFO",
		help="The logging level to use."
	)
	parser.add_argument(
		"--log-path",
		action=FullPaths,
		type=is_dir,
		default=None,
		help="The path to a directory to hold logs."
	)
	parser.add_argument(
		'--min-coverage',
		default=80,
		type=int,
		help="The minimum percent coverage required for a match [default=80]."
	)
	parser.add_argument(
		'--min-identity',
		default=80,
		type=int,
		help="The minimum percent identity required for a match [default=80]."
	)
	parser.add_argument(
		'--dupefile',
		help="Path to self-to-self lastz results for baits to remove potential duplicate probes."
	)
	parser.add_argument(
		"--regex",
		type=str,
		default="\w+_\d+_\d+",
		help="A regular expression to apply to the reference sequence names as tags in the output table.",
	)
	parser.add_argument(
		"--keep-duplicates",
		type=str,
		default=None,
		help="A optional output file in which to store those loci that appear to be duplicates.",
	)
	parser.add_argument(
		'--sqlite3',
		type=str,
		default="/usr/bin/sqlite3",
		action=FullPaths,
		help="The complete path to sqlite3"
	)
	parser.add_argument(
		'--assembler',
		choices=["trinity", "abyss"],
		default="abyss",
		help="""Please specify which assembler was used to generate the input contigs"""
	)
	


def create_probe_database(log, db, organisms, exons):
	"""Create the exon-match database"""
	log.info("Creating the exon-match database")
	conn = sqlite3.connect(db)
	c = conn.cursor()
	c.execute("PRAGMA foreign_keys = ON")
	try:
		create_string = [org + ' text' for org in organisms]
		query = "CREATE TABLE matches (exon text primary key, {0})".format(','.join(create_string))
		c.execute(query)
		query = "CREATE TABLE match_map (exon text primary key, {0})".format(','.join(create_string))
		c.execute(query)
		# convert exons to list of tuples for executemany
		all_exons = [(exon,) for exon in exons]
		c.executemany("INSERT INTO matches(exon) values (?)", all_exons)
		c.executemany("INSERT INTO match_map(exon) values (?)", all_exons)
	except sqlite3.OperationalError, e:
		log.critical("Database already exists")
		if e[0] == 'table matches already exists':
			answer = raw_input("Database already exists.  Overwrite [Y/n]? ")
			if answer == "Y" or "YES":
				os.remove(db)
				conn, c = create_probe_database(db, organisms, exons)
			else:
				sys.exit(2)
		else:
			log.critical("Cannot create database")
			raise sqlite3.OperationalError("Cannot create database")
	return conn, c


def store_lastz_results_in_db(c, matches, orientation, critter):
	"""enter matched loci in database"""
	for key, match in matches.iteritems():
		# We should have dropped all duplicates at this point
		assert len(match) == 1, "More than one match"
		item = list(match)[0]
		insert_string = "UPDATE matches SET {0} = 1 WHERE exon = '{1}'".format(critter, item)
		c.execute(insert_string)
		#pdb.set_trace()
		orient_key = "{0}({1})".format(key, list(orientation[item])[0])
		insert_string = "UPDATE match_map SET {0} = '{1}' WHERE exon = '{2}'".format(critter, orient_key, item)
		c.execute(insert_string)


def get_dupes(log, lastz_file, regex):
	"""Given a lastz_file of probes aligned to themselves, get duplicates"""
	log.info("Checking probe/bait sequences for duplicates")
	matches = defaultdict(list)
	dupes = set()
	# get names and strip probe designation since loci are the same
	for lz in lastz.Reader(lastz_file):
		target_name = new_get_probe_name(lz.name1, regex)
		query_name = new_get_probe_name(lz.name2, regex)
		matches[target_name].append(query_name)
	# see if one probe matches any other probes
	# other than the children of the locus
	for k, v in matches.iteritems():
		# if the probe doesn't match itself, we have
		# problems
		if len(v) > 1:
			for i in v:
				if i != k:
					dupes.add(k)
					dupes.add(i)
		elif k != v[0]:
			dupes.add(k)
	# make sure all names are lowercase
	return set([d.lower() for d in dupes])


def contig_count(contig):
	"""Return a count of contigs from a fasta file"""
	return sum([1 for line in open(contig, 'rU').readlines() if line.startswith('>')])


def get_organism_names_from_fasta_files(ff):
	"""Given a fasta file name, parse taxon name from file name"""
	return [os.path.basename(f).split('.')[0] for f in ff]


def check_contigs_for_dupes(matches):
	"""check for contigs that match more than 1 exon locus"""
	node_dupes = defaultdict(list)
	for node in matches:
		node_dupes[node] = len(set(matches[node]))
	dupe_set = set([node for node in node_dupes if node_dupes[node] > 1])
	return dupe_set


def check_loci_for_dupes(revmatches):
	"""Check for exon probes that match more than one contig"""
	dupe_contigs = []
	dupe_exons = []
	for exon, node in revmatches.iteritems():
		if len(node) > 1:
			dupe_contigs.extend(node)
			dupe_exons.append(exon)
	#dupe_contigs = set([i for exon, node in revmatches.iteritems() if len(node) > 1 for i in list(node)])
	#pdb.set_trace()
	return set(dupe_contigs), set(dupe_exons)


def pretty_log_output(log, critter, matches, contigs, pd, mc, exon_dupe_exons):
	"""Write some nice output to the logfile/stdout"""
	unique_matches = sum([1 for node, exon in matches.iteritems()])
	out = "{0}: {1} ({2:.2f}%) uniques of {3} contigs, {4} dupe probe matches, " + \
		"{5} exon loci removed for matching multiple contigs, {6} contigs " + \
		"removed for matching multiple exon loci"
	log.info(
		out.format(
			critter,
			unique_matches,
			float(unique_matches) / contigs * 100,
			contigs,
			len(pd),
			len(exon_dupe_exons),
			len(mc)
		)
	)


def get_contig_name(header,args):
	#parse the contig name from the header of Trinity assembled contigs"
	#args = get_args()
	match = ""
	if args.assembler == "trinity":
		match = re.search("^(c\d+_g\d+_i\d+).*", header)
	elif args.assembler == "abyss":
		match = re.search("^(\d+).*", header)
	#print "match:", match
	return match.groups()[0]


def get_kmer_value(header):
	match = re.search("^\d*\s\d*\s(\d*).*", header)
	return match.groups()[0]


def new_get_probe_name(header, regex):
	match = re.search(regex, header)
	#print match
	return match.groups()[0]


def main(args):
	#args = get_args()
	pre_regex = args.regex
	regex = re.compile("^(%s)(?:.*)" %pre_regex)
	if not os.path.isdir(args.output):
		os.makedirs(args.output)
	else:
		raise IOError("The directory {} already exists.  Please check and remove by hand.".format(args.output))
	exons = set(new_get_probe_name(seq.id, regex) for seq in SeqIO.parse(open(args.reference, 'rU'), 'fasta'))
	#print exons
	if args.dupefile:
		dupes = get_dupes(log, args.dupefile, regex)
	else:
		dupes = set()
	fasta_files = glob.glob(os.path.join(args.contigs, '*.fa*'))
	for f in fasta_files:
		replace_bad_fasta_chars = "sed -i -e '/>/! s=[K,Y,R,S,M,W,B,D,H,V,k,y,r,s,m,w,b,d,h,v]=N=g' %s" %f
		os.system(replace_bad_fasta_chars)
	#print fasta_files
	organisms = get_organism_names_from_fasta_files(fasta_files)
	#print organisms
	conn, c = create_probe_database(
		log,
		os.path.join(args.output, 'probe.matches.sqlite'),
		organisms,
		exons
	)
	log.info("Processing contig data")
	# open a file for duplicate writing, if we're interested
	if args.keep_duplicates is not None:
		dupefile = open(args.keep_duplicates, 'w')
	else:
		dupefile = None
	log.info("{}".format("-" * 65))
	kmers = {}
	for contig in sorted(fasta_files):
		critter = os.path.basename(contig).split('.')[0].replace('-', "_")
		output = os.path.join(
			args.output,
			os.path.splitext(os.path.basename(contig))[0] + '.lastz'
		)
		contigs = contig_count(contig)
		# align the probes to the contigs
		alignment = lastz.Align(
			contig,
			args.reference,
			args.min_coverage,
			args.min_identity,
			output
		)
		lzstdout, lztstderr = alignment.run()
		if lztstderr:
			raise EnvironmentError("lastz: {}".format(lztstderr))
		# parse the lastz results of the alignment
		matches = defaultdict(set)
		orientation = defaultdict(set)
		revmatches = defaultdict(set)
		probe_dupes = set()
		if not lztstderr:
			for lz in lastz.Reader(output):
				contig_name = get_contig_name(lz.name1,args)
				exon_name = new_get_probe_name(lz.name2, regex)
				if args.dupefile and exon_name in dupes:
					probe_dupes.add(exon_name)
				else:
					matches[contig_name].add(exon_name)
					orientation[exon_name].add(lz.strand2)
					revmatches[exon_name].add(contig_name)

		# we need to check nodes for dupe matches to the same probes
		contigs_matching_mult_exons = check_contigs_for_dupes(matches)
		exon_dupe_contigs, exon_dupe_exons = check_loci_for_dupes(revmatches)
		nodes_to_drop = contigs_matching_mult_exons.union(exon_dupe_contigs)
		# write out duplicates if requested
		if dupefile is not None:
			log.info("Writing duplicates file for {}".format(critter))
			if len(exon_dupe_exons) != 0:
				dupefile.write("[{} - probes hitting multiple contigs]\n".format(critter))
				for exon in exon_dupe_exons:
					dupefile.write("{}:{}\n".format(exon, ', '.join(revmatches[exon])))
				dupefile.write("\n")
			if len(contigs_matching_mult_exons) != 0:
				dupefile.write("[{} - contigs hitting multiple probes]\n".format(critter))
				for dupe in contigs_matching_mult_exons:
					dupefile.write("{}:{}\n".format(dupe, ', '.join(matches[dupe])))
				dupefile.write("\n")

		# remove dupe and/or dubious nodes/contigs
		match_copy = copy.deepcopy(matches)
		for k in match_copy.keys():
			if k in nodes_to_drop:
				del matches[k]
		#print matches
		#print lz.name1
		#get contig id
		#contig_id = re.search("^(\d*)\s\d*\s\d*.*", lz.name1).groups()[0]
		#print matches

		#added function to return the kmer count (sum of all kmers of target contigs)
		for lz in lastz.Reader(output):
			for element in matches:
				#print element, "has to match", lz[1]
				if re.search("^(\d*)\s\d*\s\d*.*", lz[1]).groups()[0] == element:
					kmer_value = get_kmer_value(lz.name1)
					kmers.setdefault(contig,[])
					kmers[contig].append(kmer_value)
		store_lastz_results_in_db(c, matches, orientation, critter)
		conn.commit()
		pretty_log_output(
			log,
			critter,
			matches,
			contigs,
			probe_dupes,
			contigs_matching_mult_exons,
			exon_dupe_exons
		)

	kmerfile = open(os.path.join(args.output,'kmer_count.txt'), 'w')

	for key in kmers:
		count = 0
		for element in kmers[key]:
			count += int(element)
		kmerfile.write("%s : %d\n" %(os.path.basename(key).split('.')[0],count))


	if dupefile is not None:
		dupefile.close()
	log.info("{}".format("-" * 65))
	log.info("The LASTZ alignments are in {}".format(args.output))
	log.info("The exon match database is in {}".format(os.path.join(args.output, "probes.matches.sqlite")))
	text = "Completed"
        
	log.info(text.center(65, "="))

	# Access the SQL file and export tab-separated text-file
	sql_file = os.path.join(args.output, 'probe.matches.sqlite')
	tsf_out = os.path.join(args.output, 'match_table.txt')
	sql_cmd = "%s -header -nullvalue '.' -separator '\t' %s \"select * from matches;\" > %s" %(args.sqlite3,sql_file,tsf_out)
	os.system(sql_cmd)

	# Create the config file for the extraction of the desired loci
	output_folder = args.output
	
        with open(os.path.join(output_folder, 'config'), 'w') as f:
		print('[Organisms]', file=f)
		for aln in glob.glob(os.path.join(output_folder, '*.lastz')):
			aln = os.path.basename(aln)
			aln = aln.split('_')[0]
			aln = aln.replace('.lastz', '')
			print(aln, file=f)

		print('\n[Loci]', file=f)
		with open(os.path.join(output_folder, 'match_table.txt')) as match_table:
			lines = match_table.readlines()
		for line in lines[1:]:
			print(line.split('\t')[0], file=f)




#if __name__ == '__main__':
#	main()
