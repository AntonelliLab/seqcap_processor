# <img src="images/secapr_logo.png" width="100"> SEquence CApture PRocessor (SECAPR)

**Original Publication: https://doi.org/10.7717/peerj.5175**

___

*We are now teaching a **free 1-week intensive course** on target enrichment data, including practical exercises for all functionalities of the SECAPR pipeline. See course overview from last year and contact tobias.andermann@bioenv.gu.se in case you are interested.*

http://www.forbio.uio.no/events/courses/2018/target_capture.html

___

## Installation instructions and documentation [click here](./documentation.html)

[![downloads](https://anaconda.org/bioconda/secapr/badges/downloads.svg)](http://bioconda.github.io/recipes/secapr/README.html)


This semi-automated pipeline aims to make the processing and analysis of sequence capture (= target enrichment) data simple and straight forward for all users. The detailed documentation and simple installation makes this pipeline accessible also for users with limited biofinformatic knowledge, while enabling user-defined processing options for the more experienced users.

We included an empirical data tutorial in the [pipeline documentation](./documentation.html), which covers the processing from raw Illumina read data into multiple seqeunce alignments (MSAs) for phylogenetic analyses, including the compiling of allele sequence MSAs. This workflow can be applied to any Illumina dataset, independently of the underlying bait set and organism group.

Some functions in this pipeline are inspired by the scripts from the [Phyluce pipeline](https://github.com/faircloth-lab/phyluce) by Braint Faircloth, which is aimed at the processing of Ultraconserved Elements (UCEs). To honour some of the ingenious ideas belonging to Brant Faircloth, and the generous sharing of all his code as open-source, we ask you to cite the Phyluce pipeline (Faircloth 2016) alongside with ours (Andermann et al. 2018), when using SECAPR.  

#### Please cite:

**Andermann T.**, Cano Á., Zizka A., Bacon C., Antonelli A. 2018. SECAPR—a bioinformatics pipeline for the rapid and user-friendly processing of targeted enriched Illumina sequences, from raw reads to alignments. PeerJ 6:e5175 https://doi.org/10.7717/peerj.5175

Faircloth B.C., PHYLUCE is a software package for the analysis of conserved genomic loci, Bioinformatics, Volume 32, Issue 5, 1 March 2016, Pages 786–788, https://doi.org/10.1093/bioinformatics/btv646
