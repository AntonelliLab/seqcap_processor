#!/usr/local/anaconda/bin/python
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

# Modified by Tobias Hofmann:
# Modifications include: 	- Standardizing script for incomplete data 
#							- Setting suitable alignment settings for standard palm contigs as default, to avoid discarding too many loci
#							- Format the sequence headers of the output alignment files to simply the sample name (no locus information in the header, only in the filename)

import os
import sys
import copy
import argparse
import tempfile
import multiprocessing
from Bio import SeqIO
from collections import defaultdict

from phyluce.helpers import FullPaths, CreateDir, is_dir, is_file, write_alignments_to_outdir
from phyluce.log import setup_logging

#import pdb

def get_args():
    parser = argparse.ArgumentParser(
        description="""Align and trim records in a monolothic FASTA file containing all loci for all taxa""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--extracted-contigs",
        required=True,
        action=FullPaths,
        type=is_file,
        help="""The fasta file containing the extracted contigs that match the target loci"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=CreateDir,
        help="""The directory in which to store the resulting alignments."""
    )
    parser.add_argument(
        "--taxa",
        required=True,
        type=int,
        help="""Number of taxa expected in each alignment."""
    )
    parser.add_argument(
        "--aligner",
        choices=["dialign", "muscle", "mafft"],
        default="mafft",
        help="""The alignment engine to use."""
    )
    parser.add_argument(
        "--output-format",
        choices=["fasta", "nexus", "phylip", "clustal", "emboss", "stockholm"],
        default="fasta",
        help="""The output alignment format.""",
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use."""
    )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    parser.add_argument(
        "--no-trim",
        action="store_true",
        default=False,
        help="""Align, but DO NOT trim alignments."""
    )
    parser.add_argument(
        "--window",
        type=int,
        default=10,
        help="""Sliding window size for trimming."""
    )
    parser.add_argument(
        "--proportion",
        type=float,
        default=0.5,
        help="""The proportion of taxa required to have sequence at alignment ends."""
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="""The proportion of residues required across the window in """ +
        """proportion of taxa."""
    )
    parser.add_argument(
        "--max-divergence",
        type=float,
        default=0.40,
        help="""The max proportion of sequence divergence allowed between any row """ +
        """of the alignment and the alignment consensus."""
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=100,
        help="""The minimum length of alignments to keep."""
    )
    parser.add_argument(
        "--ambiguous",
        action="store_true",
        default=False,
        help="""Allow reads in alignments containing N-bases."""
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="""Process alignments in parallel using --cores for alignment. """ +
        """This is the number of PHYSICAL CPUs."""
    )
    return parser.parse_args()


def build_locus_dict(log, loci, locus, record, ambiguous=False):
    if not ambiguous:
        if not "N" in record.seq:
            loci[locus].append(record)
        else:
            log.warn("Skipping {} because it contains ambiguous bases".format(locus))
    else:
        loci[locus].append(record)
    return loci


def create_locus_specific_fasta(sequences):
    fd, fasta_file = tempfile.mkstemp(suffix=".fasta")
    for seq in sequences:
        os.write(fd, seq.format("fasta"))
    os.close(fd)
    return fasta_file


def align(params):
    locus, opts = params
    name, sequences = locus
    # get additional params from params tuple
    window, threshold, notrim, proportion, divergence, min_len = opts
    fasta = create_locus_specific_fasta(sequences)
    aln = Align(fasta)
    aln.run_alignment()
    if notrim:
        aln.trim_alignment(
                method="notrim"
            )
    else:
        aln.trim_alignment(
                method="running",
                window_size=window,
                proportion=proportion,
                threshold=threshold,
                max_divergence=divergence,
                min_len=min_len
            )
    if aln.trimmed:
        sys.stdout.write(".")
    else:
        sys.stdout.write("X")
    sys.stdout.flush()
    return (name, aln)


def get_fasta_dict(log, args):
    log.info("Building the locus dictionary")
    if args.ambiguous:
        log.info("NOT removing sequences with ambiguous bases...")
    else:
        log.info("Removing ALL sequences with ambiguous bases...")
    loci = defaultdict(list)
    with open(args.extracted_contigs, "rU") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            locus = record.description.split("|")[1]
            loci = build_locus_dict(log, loci, locus, record, args.ambiguous)
    # workon a copy so we can iterate and delete
    snapshot = copy.deepcopy(loci)
    # iterate over loci to check for all species at a locus
    for locus, data in snapshot.iteritems():
        if len(data) < 3:
            del loci[locus]
            log.warn("DROPPED locus {0}. Too few taxa (N < 3).".format(locus))
    return loci



def main(args):
    # setup logging
    log, my_name = setup_logging(args)
    # create the fasta dictionary
    loci = get_fasta_dict(log, args)
    log.info("Aligning with {}".format(str(args.aligner).upper()))
    opts = [[args.window, args.threshold, args.no_trim, args.proportion, args.max_divergence, args.min_length] \
            for i in range(len(loci))]
    # combine loci and options
    params = zip(loci.items(), opts)
    log.info("Alignment begins. 'X' indicates dropped alignments (these are reported after alignment)")
    # During alignment, drop into sys.stdout for progress indicator
    # because logging in multiprocessing is more painful than what
    # we really need.  Return to logging when alignment completes.
    if args.cores > 1:
        assert args.cores <= multiprocessing.cpu_count(), "You've specified more cores than you have"
        pool = multiprocessing.Pool(args.cores)
        alignments = pool.map(align, params)
    else:
        alignments = map(align, params)
    # kick the stdout down one line since we were using sys.stdout
    print("")
    # drop back into logging
    log.info("Alignment ends")
    # write the output files
    write_alignments_to_outdir(log, args.output, alignments, args.output_format)
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))


if __name__ == '__main__':
    args = get_args()
    # globally import Align method
    if args.aligner == "muscle":
        from phyluce.muscle import Align
    elif args.aligner == "mafft":
        from phyluce.mafft import Align
    elif args.aligner == "dialign":
        from phyluce.dialign import Align
    main(args)
    output_folder = args.output
    file_format = args.output_format
    cmd = "for file in $(ls %s/*.%s); do sed -i -e 's/>\w*_[0-9]*_[0-9]*_/>/g' $file; done" %(output_folder,file_format)
    os.system(cmd)
