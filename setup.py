import sys
from setuptools import setup
import versioneer

if sys.version_info < (2, 7):
	sys.stdout.write("At least Python 2.7 is required.\n")
	sys.exit(1)

with open('README.md') as f:
	long_description = f.read()


setup(
	name = 'secapr',
	version = versioneer.get_version(),
	cmdclass = versioneer.get_cmdclass(),
	author = 'Tobias Hofmann',
	author_email = 'tobias.hofmann@bioenv.gu.se',
	url = 'https://github.com/AntonelliLab/seqcap_processor',
	description = 'Process sequence-capture fastq files into alignments for phylogenetic analyses',
	long_description = long_description,
	#license = 'BSD',
	#entry_points = {'console_scripts': ['igdiscover = igdiscover.__main__:main']},
	packages = ['secapr'],
	#package_data = {'igdiscover': ['igdiscover.yaml', 'Snakefile', 'empty.aux']},

	#install_requires = [ 	#are these python packages?
	#	'sqt>=0.8.0',
	#	'pandas>=0.16.2',
	#	'numpy',
	#	'matplotlib>=1.5.0',
	#	'snakemake>=3.9.0',
	#	'cutadapt',
	#	'seaborn>=0.6.0',
	#	'scipy>=0.16.1',
	#	'xopen>=0.1.1',
	#	'PyYAML',
	#],
	

        ## Estelle
        install_requires = [ 
                'Bio>=1.68 ',
                'phyluce>=1.5.0',
                'cogent>=0.7.7', #?
                'sqlite3>=3.8.6',
        ],
                

	classifiers = [
		"Development Status :: 4 - Beta",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 2.7",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
