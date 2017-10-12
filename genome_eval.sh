#!/bin/bash

cd /escratch4/s_11/s_11_Aug_17/genome_eval

# put the mummer executables in your $PATH
export PATH=/usr/local/mummer/latest/mummer/:$PATH

# put the Prokka executable and its dependencies in your $PATH
export PATH=/usr/local/prokka/1.11/:/usr/local/hmmer/2.3.2/:/usr/local/rnammer/latest/:/usr/local/tbl2asn/01052015:/usr/local/signalp/4.1c/:/usr/local/parallel/20150822/:$PATH



# download Ensembl MG1655 reference genome
wget -q -O ref_ecoli.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/
Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz

#unzip reference genome
gunzip -c ref_ecoli.fa.gz > ref_ecoli.fa

# download Spades Illumina E. coli assembly
wget -q http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/scaffolds.fasta

# run QUAST 3.1 on reference-based, Pacbio and Spades assemblies using Ensembl MG1655 as reference
python2.7 /usr/local/quast/3.1/quast.py -o /escratch4/s_11/s_11_Aug_17/genome_eval/quast_output -R /escratch4/s_11/s_11_Aug_17/genome_eval/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/genome_eval/samtools_consensus.fa /escratch4/s_11/s_11_Aug_17/genome_eval/pacbio_contigs.fasta
/escratch4/s_11/s_11_Aug_17/genome_eval/spades_scaffolds.fasta

# make mummerplots for reference-based, Pacbio and Spades assemblies using Ensembl MG1655 as reference
# need to run the 3 following lines for each assembly
nucmer -o /escratch4/s_11/s_11_Aug_17/genome_eval/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/genome_eval -p outputfile_prefix
delta-filter -1 outputfile_prefix.delta > outputfile_prefix.1delta
mummerplot --size large -fat --color -f --png outputfile_prefix.1delta -p outputfile_prefix

# generate Prokka genome annotations for reference-based, Pacbio and Spades assemblies using Ensembl MG1655 as reference
# need to run the 3 following lines for each assembly
prokka /escratch4/s_11/s_11_Aug_17/genome_eval --outdir prokka_directory
