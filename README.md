# SEquence CApture PRocessor (SECAPR)

## See the [documentation](./documentation.ipynb)

## Installation instructions are also found in the [documentation](./documentation.ipynb) or under [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/secapr/README.html)

This semi-automated pipeline aims to make the processing and analysis of sequence capture data simple and straight forward for all users. The detailed workflow and simple installation setup makes this workflow accessible also for users with limited biofinformatic knowledge, while enabling user-defined processing options for the more experienced users. This workflow covers the processing from raw Illumina read data into different data-types and -formats, commonly used in phylogenetic analyses (e.g. SNPs, multiple sequence alignments, haplotype sequences, etc.).

This workflow can be applied to any Illumina dataset, independently of the underlying bait set and organism group. All scripts are written in python2.7 and in some cases have to be run with that specific python version for proper functioning!

Some functions in this pipeline are modified versions of the scripts from the [Phyluce pipeline](https://github.com/faircloth-lab/phyluce) by Braint Faircloth (in those cases you find the license agreement at the beginning of the script), which is aimed at the processing of UCE data for Vertebrates. To honour some of the ingenious ideas belonging to Brant Faircloth for efficiently scripted NGS read processing, and the generous sharing of all is code as open-source, we ask you to cite his Phyluce pipeline (Faircloth et al. 2015) alongside with ours, when using these scripts.  

#### References:
Faircloth BC. 2015. PHYLUCE is a software package for the analysis of conserved genomic loci. bioRxiv. doi: 10.1101/027904.
