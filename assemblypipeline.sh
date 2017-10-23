#!/bin/bash
#$ -pe thread 4

cd /escratch4/s_11/s_11_Aug_17/assemblypipeline

#path for pacbio executables
export PATH=/usr/local/canu/1.4/Linux-amd64/bin:/usr/local/java/jdk1.8.0_74/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/gcc/5.3.0/lib64:${LD_LIBRARY_PATH}
export JAVA_HOME=/usr/local/java/jdk1.8.0_74

#path for reference assembly executables
export PATH=/usr/local/samtools/1.2/:$PATH
export PATH=/usr/local/bwa/0.7.10/:$PATH



#PACBIO ASSEMBLY

curl -L -o pacbio.fastq http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq

canu -p ecoli -d ecoli-pacbio genomeSize=4.8m -pacbio-raw pacbio.fastq useGrid=false




#REFERENCE ASSEMBLY

# download reference genome and unzip reference genome
wget -q -O ref.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz

# download reads
wget -q http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/s_6_1.fastq.gz
wget -q http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/s_6_2.fastq.gz

#unzip reference genome
gunzip -c ref.fa.gz > ref.fa

# create index of reference genome for BWA
bwa index ref.fa

# map reads to reference with BWA & output as BAM file
bwa mem -t 4 ref.fa s_6_1.fastq.gz s_6_2.fastq.gz | samtools view -b - > aln.bam

# sort BAM file
samtools sort -@ 4 -O bam -T tmp aln.bam > aln.sort.bam

# index BAM file
samtools index aln.sort.bam

# generate genotype likelihoods with `samtools mpileup` and base calls with `bcftools call`, output as VCF
samtools mpileup -u -f ref.fa aln.sort.bam | bcftools call -c -O z -o aln.sort.vcf.gz

# index VCF file
bcftools index aln.sort.vcf.gz

# generate consensus sequence
bcftools consensus -f ref.fa aln.sort.vcf.gz > consensus.fa



cd /escratch4/s_11/s_11_Aug_17/assemblypipeline


# put the mummer executables in your $PATH
export PATH=/usr/local/mummer/3.22/:$PATH

# put the Prokka executable and its dependencies in your $PATH
export PATH=/usr/local/prokka/1.11/bin/:/usr/local/hmmer/2.3.2/bin/:/usr/local/rnammer/latest/:/usr/local/tbl2asn/01052015:/usr/local/signalp/4.1c/:/usr/local/parallel/20150822/bin:$PATH

#GET SPADES ASSEMBLY

# download Spades Illumina E. coli assembly
wget -q http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/scaffolds.fasta


# get reference genome and decompress
# download Ensembl MG1655 reference genome
wget -q -O ref_ecoli.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz

#unzip reference genome
gunzip -c ref_ecoli.fa.gz > ref_ecoli.fa



#ANALYSIS OF ASSEMBLIES
# run QUAST 3.1 on reference-based, Pacbio and Spades assemblies using Ensembl MG1655 as reference
python2.7 /usr/local/quast/3.1/quast.py -o /escratch4/s_11/s_11_Aug_17/assemblypipeline/quast_output -R /escratch4/s_11/s_11_Aug_17/assemblypipeline/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/assemblypipeline/consensus.fa /escratch4/s_11/s_11_Aug_17/ecoli/assemblypipeline/ecolicontigs.fasta /escratch4/s_11/s_11_Aug_17/assemblypipeline/scaffolds.fasta

# make mummerplots for reference-based, Pacbio and Spades assemblies using Ensembl MG1655 as reference
# need to run the 3 following lines for each assembly
#reference based (refassem)
nucmer -o /escratch4/s_11/s_11_Aug_17/assemblypipeline/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/assemblypipeline/consensus.fa -p outputref_prefix
delta-filter -1 outputref_prefix.delta > outputref_prefix.1delta
mummerplot --size large -fat --color -f --png outputref_prefix.1delta -p outputref_prefix

#pacbio (canu)
nucmer -o /escratch4/s_11/s_11_Aug_17/assemblypipeline/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/ecoli/assemblypipeline/ecolicontigs.fasta -p outputpacbio_prefix
delta-filter -1 outputpacbio_prefix.delta > outputpacbio_prefix.1delta
mummerplot --size large -fat --color -f --png outputpacbio_prefix.1delta -p outputpacbio_prefix

#spades
nucmer -o /escratch4/s_11/s_11_Aug_17/assemblypipeline/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/assemblypipeline/scaffolds.fasta -p outputspades_prefix
delta-filter -1 outputspades_prefix.delta > outputspades_prefix.1delta
mummerplot --size large -fat --color -f --png outputspades_prefix.1delta -p outputspades_prefix

# generate Prokka genome annotations for reference-based, Pacbio and Spades assemblies using Ensembl MG1655 as reference
# need to run the following lines for each assembly
#reference-based
prokka /escratch4/s_11/s_11_Aug_17/assemblypipeline/consensus.fa --outdir prokka_ref

#pacbio
prokka /escratch4/s_11/s_11_Aug_17/ecoli/assemblypipeline/ecolicontigs.fasta --outdir prokka_pacbio

#spades
prokka /escratch4/s_11/s_11_Aug_17/assemblypipeline/scaffolds.fasta --outdir prokka_spades
