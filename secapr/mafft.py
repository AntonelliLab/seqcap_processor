#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 0 March 2012 09:03 PST (-0800)
"""


import os
import tempfile
import subprocess

from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

from .generic_align import GenericAlign


class Align(GenericAlign):
    """ MAFFT alignment class.  Subclass of GenericAlign which
    contains a majority of the alignment-related helper functions
    (trimming, etc.) """

    def __init__(self, input):
        """initialize, calling superclass __init__ also"""
        super(Align, self).__init__(input)

    def run_alignment(self, clean=True):
        # create results file
        fd, aln = tempfile.mkstemp(suffix='.mafft')
        # run MAFFT on the temp file
        cmd = ["mafft", "--adjustdirection", "--op",(self.input)[0], "--ep", (self.input)[1],"--maxiterate", "1000", (self.input)[2]]
        with open(aln, 'w') as aln_stdout:
            # just pass all ENV params
            proc = subprocess.Popen(cmd,
                    stderr=subprocess.PIPE,
                    stdout=aln_stdout
                )
            stderr = proc.communicate()
        #print(stderr)
        self.alignment = AlignIO.read(open(aln, 'r'), "fasta", \
                alphabet=Gapped(IUPAC.unambiguous_dna, "-"))
        if clean:
            self._clean(aln)
            self._clean((self.input)[2])
