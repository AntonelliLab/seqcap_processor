# encoding: utf-8

"""
Extract all contigs that were flagged by the find_target_contigs function. Enables to also include loci flagged as duplicates.

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


________________________________________
Modified by Tobias Hofmann (tobias.hofmann@bioenv.gu.se):
Additions include:
- function that allows to include duplicates (--include_duplicates). This enables to read the dupe-file created by the script find_target_contigs and include all those contigs that match several reference sequences.
- Regex patterns in script modified to match contigs created by Trinity and Abyss
- added --assembler flag for user to specify how contigs where created
- compatibility with improved find_target_contigs function (output of that function will be directly compatible with this script, user is not required to alter any of the intermedate files)
_________________________________________
"""

import os
import re
import sqlite3
import argparse
import itertools
import ConfigParser
import logging #add
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from phyluce.helpers import FullPaths, is_dir, is_file, get_names_from_config
#from phyluce.log import setup_logging

#import pdb
log = logging.getLogger(__name__)

def add_arguments(parser):
    parser.add_argument(
        '--contigs',
        required=True,
        action=FullPaths,
        type=is_dir,
        help='The directory containing the assembled contigs in which you searched for exon loci.',
    )
    parser.add_argument(
        '--locus-db',
        required=True,
        action=FullPaths,
        type=is_file,
        help='The SQL database file holding probe matches to targeted loci (usually "probe.matches.sqlite").'
    )
    parser.add_argument(
        '--config',
        required=True,
        action=FullPaths,
        type=is_file,
        help='The config file (named "contig") containing taxa and loci, which can be found in the folder together with the match-database sql-file.'
    )
    parser.add_argument(
        '--output',
        required=True,
        action=FullPaths,
        help='The path to the output FASTA file you want to create.'
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
        '--assembler',
        choices=["trinity", "abyss"],
        default="abyss",
        help="""Please specify which assembler was used to generate the input contigs"""
    )
    parser.add_argument(
        '--include_duplicates',
        type = str,
        default = None,
        help = "Use this flag and add the path to the duplicate-log file resulting from the 'find_target_contigs' function if you want to extract contigs which matched multiple probes. These kinds of duplicate may be of interest as they may represent long sequences spannign accross several target regions."
    )    
    
    
def flatten(lst):
    new_list = []
    for element in lst:
        element = element.split(', ')
        for subelement in element:
            new_list.append(subelement)
    return new_list

def get_nodes_for_exons(c, organism, exons, args, organism_dupe_dict, organism_orientation_dict, extend=False, notstrict=False):
    #args = get_args()
    # get only those exons we know are in the set
    exons = [("\'{0}\'").format(u) for u in exons]
    if not extend:
        query = "SELECT lower({0}), exon FROM match_map where exon in ({1})".format(organism, ','.join(exons))
    else:
        query = "SELECT lower({0}), exon FROM extended.match_map where exon in ({1})".format(organism, ','.join(exons))
    c.execute(query)
    rows = c.fetchall()
    #print(rows)
    node_dict = defaultdict()
    missing = []
    for node in rows:
        dupe_loci = flatten(organism_dupe_dict.values())
        if node[0] is not None:
            match = ""
            if args.assembler == "trinity":
                match = re.search('^(c\d+_g\d+_i\d+)\(([+-])\)', node[0])
            elif args.assembler == "abyss":
                match = re.search('^(\d+)\(([+-])\)', node[0])
            #print "match:", match.groups()
            #print node[1]
            node_dict[match.groups()[0]] = (node[1], match.groups()[1])
        elif node[0] is None and str(node[1]) in dupe_loci:
            dupe_dict = {}
            for contig, loci in organism_dupe_dict.iteritems():
                if node[1] in loci:
                    dupe_dict.setdefault(node[1],[])
                    dupe_dict[node[1]].append(contig)
            for key in dupe_dict:
                if len(dupe_dict[key]) == 1:
                    #print(str(key),dupe_dict[key][0])
                    contig = unicode(dupe_dict[key][0], "utf-8")
                    locus = unicode(str(key), "utf-8")
                    orientation = unicode(organism_orientation_dict[str(key)], "utf-8")
                    node_dict[contig] = (locus, orientation)
        elif notstrict:
            missing.append(node[1])
        else:
            raise IOError("Complete matrices should have no missing data")
    return node_dict, missing


def find_file(contigs, name):
    extensions = ['.fa', '.fasta', '.contigs.fasta', '.contigs.fa', '.gz', '.fasta.gz', '.fa.gz']
    for ext in extensions:
        reads1 = os.path.join(contigs, name) + ext
        reads2 = os.path.join(contigs, name.replace('-', '_')) + ext
        for reads in [reads1, reads2]:
            if os.path.isfile(reads):
                break
            elif os.path.isfile(reads.lower()):
                reads = reads.lower()
                break
            else:
                reads = None
        if reads is not None:
            break
    if reads is None:
        raise ValueError("Cannot find a fasta file for {} with any of the extensions ({}) ".format(
            name,
            ', '.join(extensions)
        ))
    return reads


def get_contig_name(header, args):
    #args = get_args()
    """parse the contig name from the header of either abyss/trinity assembled contigs"""
    match = ""
    if args.assembler == "trinity":
        match = re.search("^(c\d+_g\d+_i\d+).*", header)
    elif args.assembler == "abyss":
        match = re.search("^(\d+).*", header)
    return match.groups()[0]


def replace_and_remove_bases(regex, seq, count):
    new_seq_string = str(seq.seq)
    if regex.search(new_seq_string):
        new_seq_string = re.sub(regex, "", new_seq_string)
        #print "\tReplaced < 20 ambiguous bases in {0}".format(seq.id)
        count += 1
    new_seq_string = re.sub("^[acgtn]+", "", new_seq_string)
    new_seq_string = re.sub("[acgtn]+$", "", new_seq_string)
    new_seq = Seq(new_seq_string)
    new_seq_record = SeqRecord(new_seq, id=seq.id, name='', description='')
    return new_seq_record, count


def main(args):
    #args = get_args()
    # setup logging
    #log, my_name = setup_logging(args)
    # parse the config file - allowing no values (e.g. no ":" in config file)
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.optionxform = str
    config.read(args.config)
    # connect to the database
    conn = sqlite3.connect(args.locus_db)
    c = conn.cursor()
    # attach to external database, if passed as option
    organisms = get_names_from_config(config, 'Organisms')
    log.info("There are {} taxa in the match-count-config file named {}".format(
        len(organisms),
        os.path.basename(args.config)
    ))
    exons = get_names_from_config(config, 'Loci')

    dupefile = None
    if args.include_duplicates is not None:
    	dupefile = args.include_duplicates

    dupe_config = ConfigParser.RawConfigParser(allow_no_value=True)
    dupe_config.optionxform = str
    if args.include_duplicates is not None:
        dupe_config.read(dupefile)

    log.info("There are {} exon loci in the matrix".format(len(exons)))
    regex = re.compile("[N,n]{1,21}")
    out_dir = '/'.join(args.output.split('/')[:-1])
    temp_conf = os.path.join(out_dir, 'config_extended')
    incomplete_outf = open(temp_conf, 'w')
    with open(args.output, 'w') as exon_fasta_out:
        for organism in organisms:
            organism_dupe_dict = {}
            if args.include_duplicates is not None:
                dupes = dupe_config.items('%s - contigs hitting multiple probes' %organism)
                for element in dupes:
                    organism_dupe_dict.setdefault(element[0],element[1])
            organism_orientation_dict = {}
            if args.include_duplicates is not None:
                locus_orientation = dupe_config.items('%s - contig orientation' %organism)
                for element in locus_orientation:
                    organism_orientation_dict.setdefault(element[0],element[1])            
            #print (organism_dupe_dict)
            text = "Getting exon loci for {0}".format(organism)
            log.info(text.center(65, "-"))
            written = []
            # going to need to do something more generic w/ suffixes
            name = organism.replace('_', '-')
            if not organism.endswith('*'):
                reads = find_file(args.contigs, name)
                node_dict, missing = get_nodes_for_exons(c, organism, exons, args, organism_dupe_dict, organism_orientation_dict, extend=False, notstrict=True)
            count = 0
            log.info("There are {} exon loci for {}".format(len(node_dict), organism))
            log.info("Parsing and renaming contigs for {}".format(organism))
            for seq in SeqIO.parse(open(reads, 'rU'), 'fasta'):
                name = get_contig_name(seq.id,args).lower()
                #print "name:", name
                #print node_dict.keys()
                if name in node_dict.keys():
                    seq.id = "{0}_{1} |{0}".format(node_dict[name][0], organism.rstrip('*'))
                    seq.name = ''
                    seq.description = ''
                    # deal with strandedness because aligners sometimes dont, which
                    # is annoying
                    if node_dict[name][1] == '-':
                        seq.seq = seq.seq.reverse_complement()
                    # Replace any occurrences of <21 Ns in a given sequence with
                    # blanks.  These should gap out during alignment. Also, replace
                    # leading/trailing lowercase bases from velvet assemblies.
                    # Lowercase bases indicate low coverage, and these
                    # have been problematic in downstream alignments).
                    seq, count = replace_and_remove_bases(regex, seq, count)
                    exon_fasta_out.write(seq.format('fasta'))
                    #print "node_dict:", node_dict[name][0]
                    written.append(str(node_dict[name][0]))
                else:
                    pass
            if count > 0:
                log.info("Replaced <20 ambiguous bases (N) in {} contigs for {}".format(count, organism))
            if missing:
                log.info("Writing missing locus information to {}".format(temp_conf))
                incomplete_outf.write("[{0}]\n".format(organism))
                for name in missing:
                    incomplete_outf.write("{0}\n".format(name))
                    written.append(name)
            #print written
            #print exons
            # This test will result in an error if duplicates are included
            #assert set(written) == set(exons), "exon names do not match"
    text = " Completed! "
    log.info(text.center(65, "="))

#if __name__ == '__main__':
#    main()
