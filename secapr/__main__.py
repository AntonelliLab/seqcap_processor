# encoding: utf-8

#author: Tobias Andermann, tobias.andermann@bioenv.gu.se
#__main__.py created by Estelle, based on IgDiscover (https://bitbucket.org/igdiscover/igdiscover)

import os
import sys
from argparse import ArgumentParser
import logging
import warnings
from . import __version__
import importlib


__author__ = "Tobias Andermann"

# List of all subcommands. A module of the given name must exist and define
# add_arguments() and main() functions.

COMMANDS = [
		'quality_check',
        'clean_reads',
        'assemble_reads',
        'find_target_contigs',
        'align_sequences',
        'join_exons',
        'reference_assembly',
        'phase_alleles',
        'add_missing_sequences',
        'locus_selection',
		'automate_all',
		'concatenate_alignments',
		'paralogs_to_ref'
]


def main(arguments=None):
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	parser = ArgumentParser(description=__doc__, prog='secapr')
	parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

	subparsers = parser.add_subparsers()
	for command_name in COMMANDS:
		module = importlib.import_module('.' + command_name, 'secapr')
		subparser = subparsers.add_parser(command_name,
			help=module.__doc__.split('\n')[1], description=module.__doc__)
		subparser.set_defaults(func=module.main)
		module.add_arguments(subparser)
        
	args = parser.parse_args(arguments)
	if not hasattr(args, 'func'):
		parser.error('Please provide the name of a subcommand to run')
	else:
		args.func(args)


if __name__ == '__main__':
	main()
