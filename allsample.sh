#!/bin/bash
#$ -pe thread 6

cd /escratch4/s_11/s_11_Aug_17/allsample

#path for reference assembly executables
export PATH=/usr/local/samtools/1.2/:$PATH
export PATH=/usr/local/bwa/0.7.10/:$PATH

# path for emboss points to version 6.5.7
export PATH=/usr/local/emboss/latest/bin/:$PATH

# path for RAxML
export PATH=/usr/local/raxml/8.2.4/:$PATH

# get reference genome
wget -q -O ref.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_4_collection/salmonella_enterica_subsp_enterica_serovar_enteritidis_str_p125109/dna/Salmonella_enterica_subsp_enterica_serovar_enteritidis_str_p125109.ASM950v1.dna.chromosome.Chromosome.fa.gz

#curl -o meta.tsv "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB634&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,sample_alias,experiment_accession,run_accession,scientific_name,instrument_model,library_layout,library_source,library_selection,read_count,base_count,experiment_title,fastq_ftp" | grep -v PacBio > meta.tsv

# unzip reference genome
gunzip -c ref.fa.gz > ref.fa

# create index of reference genome for BWA
bwa index ref.fa

#for i in ERR369378 ERR369348 ERR338264 ERR338265 ERR338272 ERR338273 ERR338280 ERR338281 ERR338288 ERR338255 ERR338296 ERR338256 ERR338304 ERR338257 ERR338914 ERR338922 ERR338929 ERR338930 ERR338956 ERR338957 ERR338959 ERR338960 ERR338967 ERR338968 ERR338969 ERR338970 ERR338913 ERR353345 ERR374229 ERR311263 ERR311256 ERR311257 ERR311258 ERR311259 ERR311261 ERR311262 ERR338266 ERR338274 ERR338282 ERR338289 ERR338303 ERR338306 ERR338312 ERR338314 ERR338258 ERR338915 ERR338916 ERR338923 ERR338906 ERR338928
for i in `cut -f6 meta.tsv | grep -v run_accession`
do
  # get reads from EBI
  wget -q -O ${i}_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${i:0:6}/${i}/${i}_1.fastq.gz
  wget -q -O ${i}_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${i:0:6}/${i}/${i}_2.fastq.gz

  # map reads to reference with BWA & output as BAM file
  bwa mem -t 4 ref.fa ${i}_1.fastq.gz ${i}_2.fastq.gz | samtools view -b - > ${i}.aln.bam

  # sort BAM file
  samtools sort -@ 4 -O bam -T tmp ${i}.aln.bam > ${i}.aln.sort.bam

  # index BAM file
  samtools index ${i}.aln.sort.bam

  # generate genotype likelihoods with `samtools mpileup` and base calls with `bcftools call`, output as VCF
  samtools mpileup -I -u -f ref.fa ${i}.aln.sort.bam | bcftools call -c -O z -o ${i}.aln.sort.vcf.gz

  # index VCF files
  bcftools index ${i}.aln.sort.vcf.gz

  # generate consensus sequences
  bcftools consensus -f ref.fa ${i}.aln.sort.vcf.gz | sed "s/>.*/>$i/g" > ${i}.consensus.fa

done


# concatenate consensus files into one multi-fasta files
cat /escratch4/s_11/s_11_Aug_17/allsample/*.consensus.fa > /escratch4/s_11/s_11_Aug_17/allsample/allconsensus.fasta


# convert allconsensus.fasta to .phy format
seqret -sequence fasta::allconsensus.fasta -outseq phylip::allconsensus.phy


# run RAxML GTR with + I + G model and 100 bootstrap pseudoreplicate analyses of the alignment data
#-T is number of threads
#-f a is rapid Bootstrap analysis and search for best­scoring ML tree in one program run 
#-G enables the ML­based evolutionary placement algorithm heuristics by specifying a threshold value of 0.1 (10% of branches considered for thorough insertions)
#-I is a posteriori bootstopping analysis
#-x specifies an integer number (random seed) and turn on rapid bootstrapping
#-p specifies a random number seed for the parsimony inferences
#-m is the model
#-n is the name of the new files
#-s is the input file
/usr/local/raxml/latest/raxmlHPC-PTHREADS -T 6 -f a -x 12345 -p 4523 -m GTRGAMMAI -n enteritidis -s allconsensus.phy -# 100
