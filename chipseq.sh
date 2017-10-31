#!/bin/bash

#$ -pe thread 6

cd /escratch4/s_11/s_11_Aug_17/chipseq

#path to bowtie2 2.2.3
export PATH=/usr/local/bowtie2/2.2.3:$PATH

#path to Macs 1.4.2, "latest" points to this version
export PATH=/usr/local/macs/latest/bin:$PATH

#path to samtools
export PATH=/usr/local/samtools/1.2:$PATH

#path to BEDTools 2.25.0
export PATH=/usr/local/bedtools/2.25.0/bin:$PATH

#path to MEME 4.10.0
export LD_LIBRARY_PATH=/usr/local/mpich2/1.4.1p1/gcc_4.5.3/lib:$PATH
export PATH=/usr/local/meme/latest/bin:$PATH

# download Ensembl MG1655 reference genome
wget -q -O ref.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz

# unzip reference genome
gunzip -c ref.fa.gz > ref.fa

# download FNR anaerobic chip DNA from EBI as fastq: SRR576933  from https://www.ebi.ac.uk/ena/data/view/SRR576933
wget -q -O SRR576933.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR576/SRR576933/SRR576933.fastq.gz

# download FNR anaerobic input DNA from EBI as fastq: SRR576938   from https://www.ebi.ac.uk/ena/data/view/SRR576938)
wget -q -O SRR576938.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR576/SRR576938/SRR576938.fastq.gz

# build forward and backward indices of reference genome for mapping reads with bowtie
bowtie2-build -f ref.fa ref_index_prefix

# map the ChIP-seq reads to the reference genome using bowtie2
bowtie2 -x ref_index_prefix -U SRR576933.fastq.gz -k 1 --threads 6 | samtools view -b - > chip.bam

# map the input reads to the reference genome using bowtie2
bowtie2 -x ref_index_prefix -U SRR576938.fastq.gz -k 1 --threads 6 | samtools view -b - > input.bam

# predict FNR binding regions using MACS
macs14 -t chip.bam -c input.bam -f BAM -g 4641652 -n "FNR" --bw=400 --keep-dup=1 --bdg --single-profile &> MACS.out

# get fasta sequences of the extended peak binding regions in the summit file using BEDtools:
#index the reference genome with `samtools faidx`
samtools faidx ref.fa

# add ±100 bp to the FNR_summits.bed using `bedtools slop`
bedtools slop -b 100 -i FNR_summits.bed -g ref.fa > FNR_summits_100.bed

#extract fasta sequences for ±100 bp FNR_summits using `bedtools getfasta`
bedtools getfasta -fo FNR_summits_100.fa -fi ref.fa -bed FNR_summits_100.bed

# predict FNR binding motif using MEME
meme -dna -mod zoops -revcomp FNR_summits_100.fa -minw 5 -maxw 15 -o FNR
