#!/bin/bash

cd /escratch4/s_11/s_11_Aug_17/refassem

time /usr/local/bwa/latest/bwa
time /usr/local/samtools/latest/bin/samtools


# download reference genome and unzip reference genome
wget -q -o ref.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz | gunzip -c > ecoli_MG1655.fa


# download reads
wget -q http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/s_6_1.fastq.gz
wget -q http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/s_6_2.fastq.gz


# create index of reference genome for BWA
bwa index ecoli_MG1655.fa


# map reads to reference with BWA & output as BAM file
bwa mem -t 4 ecoli_MG1655.fa s_6_1.fastq.gz s_6_2.fastq.gz | samtools view -b - > aln.bam


# sort BAM file
samtools sort -@ 4 -O bam -T tmp aln.bam > aln.sort.bam


# index BAM file
samtools index aln.sort.bam


# generate genotype likelihoods with `samtools mpileup` and base calls with `bcftools call`, output as VCF
samtools mpileup -u -f ecoli_MG1655.fa aln.sort.bam | bcftools call -c -v -O z -o aln.sort.vcf.gz


# index VCF file
bcftools index aln.sort.vcf.gz


# generate consensus sequence
bcftools consensus -f ecoli_MG1655.fa aln.sort.vcf.gz > consensus.fa
