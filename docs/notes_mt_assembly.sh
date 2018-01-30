# downloaded mitohcondrial genome sequence form NCBI nucleotide database
# https://www.ncbi.nlm.nih.gov/nuccore/KX028885.1?report=fasta
# Cocos nucifera mitochondrion, complete genome
# GenBank: KX028885.1

mkdir data/processed/assembling_mt_genome
cp data/raw/mitochondrial_genome_cocos_nucifera.fasta data/processed/assembling_mt_genome
bwa index data/processed/assembling_mt_genome/mitochondrial_genome_cocos_nucifera.fasta

# bwa mem default settings
mkdir data/processed/assembling_mt_genome/conservative
bwa mem -M data/processed/assembling_mt_genome/mitochondrial_genome_cocos_nucifera.fasta data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ1.fastq data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ2.fastq > data/processed/assembling_mt_genome/mapping_results_1074.sam
samtools view -b -o data/processed/assembling_mt_genome/mapping_results_1074.bam -S data/processed/assembling_mt_genome/mapping_results_1074.sam
samtools sort -o data/processed/assembling_mt_genome/1074_mitochondrion_remapped_sorted_bam.bam data/processed/assembling_mt_genome/mapping_results_1074.bam
samtools index data/processed/assembling_mt_genome/1074_mitochondrion_remapped_sorted_bam.bam
mv data/processed/assembling_mt_genome/*1074* data/processed/assembling_mt_genome/conservative
#--> too conservative

# liberal
mkdir data/processed/assembling_mt_genome/liberal
# bwa mem -B 1
mkdir data/processed/assembling_mt_genome/liberal/B_1
bwa mem -B 1 -M data/processed/assembling_mt_genome/mitochondrial_genome_cocos_nucifera.fasta data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ1.fastq data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ2.fastq > data/processed/assembling_mt_genome/liberal/B_1/mapping_results_1074.sam
samtools view -b -o data/processed/assembling_mt_genome/liberal/B_1/mapping_results_1074.bam -S data/processed/assembling_mt_genome/liberal/B_1/mapping_results_1074.sam
samtools sort -o data/processed/assembling_mt_genome/liberal/B_1/1074_mitochondrion_remapped_sorted_bam.bam data/processed/assembling_mt_genome/liberal/B_1/mapping_results_1074.bam
samtools index data/processed/assembling_mt_genome/liberal/B_1/1074_mitochondrion_remapped_sorted_bam.bam

# bwa mem -B 2 -k 12 -L 25
mkdir data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25
bwa mem -B 2 -k 12 -L 25 -M data/processed/assembling_mt_genome/mitochondrial_genome_cocos_nucifera.fasta data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ1.fastq data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ2.fastq > data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25/mapping_results_1074.sam
samtools view -b -o data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25/mapping_results_1074.bam -S data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25/mapping_results_1074.sam
samtools sort -o data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25/1074_mitochondrion_remapped_sorted_bam.bam data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25/mapping_results_1074.bam
samtools index data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25/1074_mitochondrion_remapped_sorted_bam.bam


# bwa mem -B 1 -k 12 -L 25 -r 0.5
mkdir data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25_r_05
bwa mem -B 1 -k 12 -L 25 -r 0.5 -M data/processed/assembling_mt_genome/mitochondrial_genome_cocos_nucifera.fasta data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ1.fastq data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ2.fastq > data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25_r_05/mapping_results_1074.sam
samtools view -b -o data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25_r_05/mapping_results_1074.bam -S data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25_r_05/mapping_results_1074.sam
samtools sort -o data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25_r_05/1074_mitochondrion_remapped_sorted_bam.bam data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25_r_05/mapping_results_1074.bam
samtools index data/processed/assembling_mt_genome/liberal/B_2_k_12_L_25_r_05/1074_mitochondrion_remapped_sorted_bam.bam


# bwa mem -B 1 -k 20 -L 50
mkdir data/processed/assembling_mt_genome/liberal/B_2_k_20_L_50
bwa mem -B 1 -k 12 -L 50 -M data/processed/assembling_mt_genome/mitochondrial_genome_cocos_nucifera.fasta data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ1.fastq data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ2.fastq > data/processed/assembling_mt_genome/liberal/B_2_k_20_L_50/mapping_results_1074.sam
samtools view -b -o data/processed/assembling_mt_genome/liberal/B_2_k_20_L_50/mapping_results_1074.bam -S data/processed/assembling_mt_genome/liberal/B_2_k_20_L_50/mapping_results_1074.sam
samtools sort -o data/processed/assembling_mt_genome/liberal/B_2_k_20_L_50/1074_mitochondrion_remapped_sorted_bam.bam data/processed/assembling_mt_genome/liberal/B_2_k_20_L_50/mapping_results_1074.bam
samtools index data/processed/assembling_mt_genome/liberal/B_2_k_20_L_50/1074_mitochondrion_remapped_sorted_bam.bam

# bwa mem -B 4 -k 20 -L 50
mkdir data/processed/assembling_mt_genome/liberal/B_4_k_20_L_50
bwa mem -B 4 -k 12 -L 50 -M data/processed/assembling_mt_genome/mitochondrial_genome_cocos_nucifera.fasta data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ1.fastq data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ2.fastq > data/processed/assembling_mt_genome/liberal/B_4_k_20_L_50/mapping_results_1074.sam
samtools view -b -o data/processed/assembling_mt_genome/liberal/B_4_k_20_L_50/mapping_results_1074.bam -S data/processed/assembling_mt_genome/liberal/B_4_k_20_L_50/mapping_results_1074.sam
samtools sort -o data/processed/assembling_mt_genome/liberal/B_4_k_20_L_50/1074_mitochondrion_remapped_sorted_bam.bam data/processed/assembling_mt_genome/liberal/B_4_k_20_L_50/mapping_results_1074.bam
samtools index data/processed/assembling_mt_genome/liberal/B_4_k_20_L_50/1074_mitochondrion_remapped_sorted_bam.bam


# bwa mem -B 4 -k 20 -L 100
mkdir data/processed/assembling_mt_genome/liberal/B_4_k_20_L_100
bwa mem -B 4 -k 20 -L 100 -M data/processed/assembling_mt_genome/mitochondrial_genome_cocos_nucifera.fasta data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ1.fastq data/processed/cleaned_trimmed_reads/1074_clean/1074_clean-READ2.fastq > data/processed/assembling_mt_genome/liberal/B_4_k_20_L_100/mapping_results_1074.sam
samtools view -b -o data/processed/assembling_mt_genome/liberal/B_4_k_20_L_100/mapping_results_1074.bam -S data/processed/assembling_mt_genome/liberal/B_4_k_20_L_100/mapping_results_1074.sam
samtools sort -o data/processed/assembling_mt_genome/liberal/B_4_k_20_L_100/1074_mitochondrion_remapped_sorted_bam.bam data/processed/assembling_mt_genome/liberal/B_4_k_20_L_100/mapping_results_1074.bam
samtools index data/processed/assembling_mt_genome/liberal/B_4_k_20_L_100/1074_mitochondrion_remapped_sorted_bam.bam
