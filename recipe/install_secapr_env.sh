conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda  
conda create -y -n secapr_env
sleep 5
source activate secapr_env
sleep 5
conda activate secapr_env
sleep 5
conda install -y python=3.8
conda install -y pandas
conda install -y matplotlib-base
conda install -y biopython
conda install -y fastqc
conda install -y fastp=0.23
conda install -y spades=3.15.2
conda install -y blast
conda install -y mafft
conda install -y muscle
conda install -y emboss
conda install -y bwa
conda install -y samtools==1.3.1
conda install -y trimal
conda install -y secapr
#pip install https://github.com/AntonelliLab/seqcap_processor/archive/refs/tags/v2.2.5.tar.gz
