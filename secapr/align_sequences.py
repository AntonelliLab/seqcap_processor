# encoding: utf-8
"""
Align sequences and produce separate alignment file for each locus, containing the seqeunces of all taxa.

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
- Standardizing script for incomplete data 
- More forgiving default options for non-UCE datasets
- Format the sequence headers of the output alignment files to simply the sample name (no locus information in the header, only in the filename)
________________________________________
"""

from __future__ import print_function
import os
import sys
import copy
import tempfile
import multiprocessing
import logging
import pickle
from Bio import SeqIO
from collections import defaultdict
from phyluce.helpers import FullPaths, CreateDir, is_dir, is_file, write_alignments_to_outdir


log = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "--sequences",
        required=True,
        action=FullPaths,
        type=is_file,
        help="""The fasta file containing individual sequences for several samples and loci"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=CreateDir,
        help="""The directory in which to store the resulting alignments."""
    )
    parser.add_argument(
        "--aligner",
        choices=["muscle", "mafft"],
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
        default=80,
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
    window, threshold, notrim, proportion, divergence, min_len, align_class = opts
    fasta = create_locus_specific_fasta(sequences)
    aln = align_class(fasta)
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
    with open(args.sequences, "rU") as infile:
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
    if args.aligner == "muscle":
        from phyluce.muscle import Align as align_class
    elif args.aligner == "mafft":
        from phyluce.mafft import Align as align_class
    
    # create the fasta dictionary
    loci = get_fasta_dict(log, args)
    log.info("Aligning with {}".format(str(args.aligner).upper()))
    opts = [[args.window, args.threshold, args.no_trim, args.proportion, args.max_divergence, args.min_length, align_class] \
            for _ in loci]
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

    #import pickle
    #with open('/Users/tobias/Desktop/alignments.pickle', 'wb') as handle:
    #    pickle.dump(alignments, handle, protocol=pickle.HIGHEST_PROTOCOL)
    #with open('/Users/tobias/Desktop/alignments.pickle', 'rb') as handle:
    #    alignments = pickle.load(handle)

    # kick the stdout down one line since we were using sys.stdout
    print("")
    # drop back into logging
    log.info("Alignment ends")
    # write the output files
    for name, alignment in alignments:
        if alignment.trimmed:
            for t in alignment.trimmed:
                t.id = t.description.split('|')[0].split('_')[-1]
                t.name = t.id
                t.description = ''
    write_alignments_to_outdir(log, args.output, alignments, args.output_format)
    try:
        #input_folder = '/'.join(args.sequences.split('/')[:-2])
        pickle_path = os.path.join(args.output,'.secapr_files')
        if not os.path.exists(pickle_path):
            os.makedirs(pickle_path)
        with open(os.path.join(pickle_path,'sequence_origin.pickle'), 'wb') as handle:
            pickle.dump(args.sequences, handle, protocol=pickle.HIGHEST_PROTOCOL)
        # end
        text = " Completed! "
        log.info(text.center(65, "="))
    except:
        print('Could not pass origin of sequences to %s'%pickle_path)
  
