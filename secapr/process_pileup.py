#author: Tobias Andermann, tobias.andermann@bioenv.gu.se

import csv
import os
import sys
import re
import glob
import shutil
import argparse
import ConfigParser
import commands
import subprocess
from Bio import SeqIO

from .utils import CompletePath


# Get arguments
def get_args():
    parser = argparse.ArgumentParser(
        description="Mask all positions with low read coverage or strange coverage (many reads beginning or ending at same position) as uncertainties.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--pileup',
        required=True,
        action=CompletePath,
        default=None,
        help='The name of the file containing the samtools mpileup output'
    )
    parser.add_argument(
        '--cutoff',
        type=int,
        default=6,
        help='The minimum read depth that you want to accept'
    	)
    return parser.parse_args()

args = get_args()

file = args.pileup


def count_letters(word):
    GOOD_LETTERS = "actgACTG"
    return len([letter for letter in word if letter in GOOD_LETTERS])


with open(file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        if count_letters(row[4]) < args.cutoff:
            row[4] = "N" * count_letters(row[4])
        print row
