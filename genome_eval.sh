#!/bin/bash

cd /escratch4/s_11/s_11_Aug_17/genome_eval

# put the mummer executables in your $PATH
export PATH=/usr/local/mummer/3.22/:$PATH

# put the Prokka executable and its dependencies in your $PATH
PATH=/usr/local/prokka/1.11/bin/:/usr/local/hmmer/2.3.2/bin/:/usr/local/rnammer/latest/:/usr/local/tbl2asn/01052015:/usr/local/signalp/4.1c/:/usr/local/parallel/20150822/bin:$PATH

# download Ensembl MG1655 reference genome
wget -q -O ref_ecoli.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz

#unzip reference genome
gunzip -c ref_ecoli.fa.gz > ref_ecoli.fa

# download Spades Illumina E. coli assembly
wget -q http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/scaffolds.fasta

# run QUAST 3.1 on reference-based, Pacbio and Spades assemblies using Ensembl MG1655 as reference
python2.7 /usr/local/quast/3.1/quast.py -o /escratch4/s_11/s_11_Aug_17/genome_eval/quast_output -R /escratch4/s_11/s_11_Aug_17/genome_eval/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/refassem/consensus.fa /escratch4/s_11/s_11_Aug_17/ecoli/ecoli-pacbio/ecoli.contigs.fasta /escratch4/s_11/s_11_Aug_17/genome_eval/scaffolds.fasta

# make mummerplots for reference-based, Pacbio and Spades assemblies using Ensembl MG1655 as reference
# need to run the 3 following lines for each assembly
#reference based (refassem)
nucmer -o /escratch4/s_11/s_11_Aug_17/genome_eval/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/refassem/consensus.fa -p outputref_prefix
delta-filter -1 outputref_prefix.delta > outputref_prefix.1delta
mummerplot --size large -fat --color -f --png outputref_prefix.1delta -p outputref_prefix

#pacbio (canu)
nucmer -o /escratch4/s_11/s_11_Aug_17/genome_eval/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/ecoli/ecoli-pacbio/ecoli.contigs.fasta -p outputpacbio_prefix
delta-filter -1 outputpacbio_prefix.delta > outputpacbio_prefix.1delta
mummerplot --size large -fat --color -f --png outputpacbio_prefix.1delta -p outputpacbio_prefix

#spades
nucmer -o /escratch4/s_11/s_11_Aug_17/genome_eval/ref_ecoli.fa /escratch4/s_11/s_11_Aug_17/genome_eval/scaffolds.fasta -p outputspades_prefix
delta-filter -1 outputspades_prefix.delta > outputspades_prefix.1delta
mummerplot --size large -fat --color -f --png outputspades_prefix.1delta -p outputspades_prefix

# generate Prokka genome annotations for reference-based, Pacbio and Spades assemblies using Ensembl MG1655 as reference
# need to run the following lines for each assembly
#reference-based
prokka /escratch4/s_11/s_11_Aug_17/refassem/consensus.fa --outdir prokka_directory_ref

#pacbio
prokka /escratch4/s_11/s_11_Aug_17/ecoli/ecoli-pacbio/ecoli.contigs.fasta --outdir prokka_directory_pacbio

#spades
prokka /escratch4/s_11/s_11_Aug_17/genome_eval/scaffolds.fasta --outdir prokka_directory_spades

#clean up script by removing files
