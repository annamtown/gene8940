#!/bin/bash

cd /escratch4/s_11/s_11_Aug_17/project

export PATH=/usr/local/samtools/1.2/:$PATH
export PATH=/usr/local/bwa/0.7.10/:$PATH

# get reference genome
wget -q -O ref.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_4_collection/salmonella_enterica_subsp_enterica_serovar_enteritidis_str_p125109/dna/Salmonella_enterica_subsp_enterica_serovar_enteritidis_str_p125109.ASM950v1.dna.chromosome.Chromosome.fa.gz

# unzip reference genome
gunzip -c ref.fa.gz > ref.fa

# download isolates' DNA from EBI as fastq
wget -q -O ERR369378_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR369/ERR369378/ERR369378_1.fastq.gz
wget -q -O ERR369378_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR369/ERR369378/ERR369378_2.fastq.gz

wget -q -O ERR369348_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR369/ERR369348/ERR369348_1.fastq.gz
wget -q -O ERR369348_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR369/ERR369348/ERR369348_2.fastq.gz

wget -q -O ERR338264_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338264/ERR338264_1.fastq.gz
wget -q -O ERR338264_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338264/ERR338264_2.fastq.gz

wget -q -O ERR338265_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338265/ERR338265_1.fastq.gz
wget -q -O ERR338265_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338265/ERR338265_2.fastq.gz

wget -q -O ERR338272_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338272/ERR338272_1.fastq.gz
wget -q -O ERR338272_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338272/ERR338272_2.fastq.gz

wget -q -O ERR338273_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338273/ERR338273_1.fastq.gz
wget -q -O ERR338273_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338273/ERR338273_2.fastq.gz

wget -q -O ERR338280_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338280/ERR338280_1.fastq.gz
wget -q -O ERR338280_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338280/ERR338280_2.fastq.gz

wget -q -O ERR338281_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338281/ERR338281_1.fastq.gz
wget -q -O ERR338281_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338281/ERR338281_2.fastq.gz

wget -q -O ERR338288_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338288/ERR338288_1.fastq.gz
wget -q -O ERR338288_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338288/ERR338288_2.fastq.gz

wget -q -O ERR338255_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338255/ERR338255_1.fastq.gz
wget -q -O ERR338255_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR338/ERR338255/ERR338255_2.fastq.gz

# create index of reference genome for BWA
bwa index ref.fa

# map reads to reference with BWA & output as BAM file
bwa mem -t 4 ref.fa ERR369378_1.fastq.gz ERR369378_2.fastq.gz | samtools view -b - > ERR369378.aln.bam

bwa mem -t 4 ref.fa ERR369348_1.fastq.gz ERR369348_2.fastq.gz | samtools view -b - > ERR369348.aln.bam

bwa mem -t 4 ref.fa ERR338264_1.fastq.gz ERR338264_2.fastq.gz | samtools view -b - > ERR338264.aln.bam

bwa mem -t 4 ref.fa ERR338265_1.fastq.gz ERR338265_2.fastq.gz | samtools view -b - > ERR338265.aln.bam

bwa mem -t 4 ref.fa ERR338272_1.fastq.gz ERR338272_2.fastq.gz | samtools view -b - > ERR338272.aln.bam

bwa mem -t 4 ref.fa ERR338273_1.fastq.gz ERR338273_2.fastq.gz | samtools view -b - > ERR338273.aln.bam

bwa mem -t 4 ref.fa ERR338280_1.fastq.gz ERR338280_2.fastq.gz | samtools view -b - > ERR338280.aln.bam

bwa mem -t 4 ref.fa ERR338281_1.fastq.gz ERR338281_2.fastq.gz | samtools view -b - > ERR338281.aln.bam

bwa mem -t 4 ref.fa ERR338288_1.fastq.gz ERR338288_2.fastq.gz | samtools view -b - > ERR338288.aln.bam

bwa mem -t 4 ref.fa ERR338255_1.fastq.gz ERR338255_2.fastq.gz | samtools view -b - > ERR338255.aln.bam


# sort BAM file
samtools sort -@ 4 -O bam -T tmp ERR369378.aln.bam > ERR369378.aln.sort.bam

samtools sort -@ 4 -O bam -T tmp ERR369348.aln.bam > ERR369348.aln.sort.bam

samtools sort -@ 4 -O bam -T tmp ERR338264.aln.bam > ERR338264.aln.sort.bam

samtools sort -@ 4 -O bam -T tmp ERR338265.aln.bam > ERR338265.aln.sort.bam

samtools sort -@ 4 -O bam -T tmp ERR338272.aln.bam > ERR338272.aln.sort.bam

samtools sort -@ 4 -O bam -T tmp ERR338273.aln.bam > ERR338273.aln.sort.bam

samtools sort -@ 4 -O bam -T tmp ERR338280.aln.bam > ERR338280.aln.sort.bam

samtools sort -@ 4 -O bam -T tmp ERR338281.aln.bam > ERR338281.aln.sort.bam

samtools sort -@ 4 -O bam -T tmp ERR338288.aln.bam > ERR338288.aln.sort.bam

samtools sort -@ 4 -O bam -T tmp ERR338255.aln.bam > ERR338255.aln.sort.bam


# index BAM file
samtools index ERR369378.aln.sort.bam

samtools index ERR369348.aln.sort.bam

samtools index ERR338264.aln.sort.bam

samtools index ERR338265.aln.sort.bam

samtools index ERR338272.aln.sort.bam

samtools index ERR338273.aln.sort.bam

samtools index ERR338280.aln.sort.bam

samtools index ERR338281.aln.sort.bam

samtools index ERR338288.aln.sort.bam

samtools index ERR338255.aln.sort.bam


# generate genotype likelihoods with `samtools mpileup` and base calls with `bcftools call`, output as VCF;
samtools mpileup -u -f ref.fa ERR369378.aln.sort.bam | bcftools call -c -O z -o ERR369378.aln.sort.vcf.gz

samtools mpileup -u -f ref.fa ERR369348.aln.sort.bam | bcftools call -c -O z -o ERR369348.aln.sort.vcf.gz

samtools mpileup -u -f ref.fa ERR338264.aln.sort.bam | bcftools call -c -O z -o ERR338264.aln.sort.vcf.gz

samtools mpileup -u -f ref.fa ERR338265.aln.sort.bam | bcftools call -c -O z -o ERR338265.aln.sort.vcf.gz

samtools mpileup -u -f ref.fa ERR338272.aln.sort.bam | bcftools call -c -O z -o ERR338272.aln.sort.vcf.gz

samtools mpileup -u -f ref.fa ERR338273.aln.sort.bam | bcftools call -c -O z -o ERR338273.aln.sort.vcf.gz

samtools mpileup -u -f ref.fa ERR338280.aln.sort.bam | bcftools call -c -O z -o ERR338280.aln.sort.vcf.gz

samtools mpileup -u -f ref.fa ERR338281.aln.sort.bam | bcftools call -c -O z -o ERR338281.aln.sort.vcf.gz

samtools mpileup -u -f ref.fa ERR338288.aln.sort.bam | bcftools call -c -O z -o ERR338288.aln.sort.vcf.gz

samtools mpileup -u -f ref.fa ERR338255.aln.sort.bam | bcftools call -c -O z -o ERR338255.aln.sort.vcf.gz

# index VCF files
bcftools index ERR369378.aln.sort.vcf.gz

bcftools index ERR369348.aln.sort.vcf.gz

bcftools index ERR338264.aln.sort.vcf.gz

bcftools index ERR338265.aln.sort.vcf.gz

bcftools index ERR338272.aln.sort.vcf.gz

bcftools index ERR338273.aln.sort.vcf.gz

bcftools index ERR338280.aln.sort.vcf.gz

bcftools index ERR338281.aln.sort.vcf.gz

bcftools index ERR338288.aln.sort.vcf.gz

bcftools index ERR338255.aln.sort.vcf.gz


# generate consensus sequences
bcftools consensus -f ref.fa ERR369378.aln.sort.vcf.gz > ERR369378.consensus.fa

bcftools consensus -f ref.fa ERR369348.aln.sort.vcf.gz > ERR369348.consensus.fa

bcftools consensus -f ref.fa ERR338264.aln.sort.vcf.gz > ERR338264.consensus.fa

bcftools consensus -f ref.fa ERR338265.aln.sort.vcf.gz > ERR338265.consensus.fa

bcftools consensus -f ref.fa ERR338272.aln.sort.vcf.gz > ERR338272.consensus.fa

bcftools consensus -f ref.fa ERR338273.aln.sort.vcf.gz > ERR338273.consensus.fa

bcftools consensus -f ref.fa ERR338280.aln.sort.vcf.gz > ERR338280.consensus.fa

bcftools consensus -f ref.fa ERR338281.aln.sort.vcf.gz > ERR338281.consensus.fa

bcftools consensus -f ref.fa ERR338288.aln.sort.vcf.gz > ERR338288.consensus.fa

bcftools consensus -f ref.fa ERR338255.aln.sort.vcf.gz > ERR338255.consensus.fa
