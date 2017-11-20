#!/bin/bash
#$ -pe thread 4

cd /escratch4/s_11/s_11_Aug_17/project

#path for reference assembly executables
export PATH=/usr/local/samtools/1.2/:$PATH
export PATH=/usr/local/bwa/0.7.10/:$PATH

# path for emboss points to version 6.5.7
export PATH=/usr/local/emboss/latest/bin/:$PATH

# path for RAxML
export PATH=/usr/local/raxml/8.2.4/:$PATH

# get reference genome
wget -q -O ref.fa.gz ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_4_collection/salmonella_enterica_subsp_enterica_serovar_enteritidis_str_p125109/dna/Salmonella_enterica_subsp_enterica_serovar_enteritidis_str_p125109.ASM950v1.dna.chromosome.Chromosome.fa.gz

curl -s "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB634&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,sample_alias,experiment_accession,run_accession,scientific_name,instrument_model,library_layout,library_source,library_selection,read_count,base_count,experiment_title,fastq_ftp" | grep -v PacBio > meta.tsv

# unzip reference genome
gunzip -c ref.fa.gz > ref.fa

# create index of reference genome for BWA
bwa index ref.fa

for i in ERR369378 ERS255514 ERS255513 ERS235836 ERS235837 ERS235844 ERS235845 ERS235852 ERS235853 ERS235860 ERS235827 ERS235868 ERS235828 ERS235876 ERS235829 ERS235900 ERS235908 ERS235915 ERS235916 ERS235942 ERS235943 ERS235945 ERS235946 ERS235953 ERS235954 ERS235955 ERS235956 ERS235899 ERS325071 ERS248547 ERS217265 ERS217258 ERS217259 ERS217260 ERS217261 ERS217263 ERS217264 ERS235838 ERS235846 ERS235854 ERS235861 ERS235875 ERS235878 ERS235884 ERS235886 ERS235830 ERS235901 ERS235902 ERS235909 ERS235892 ERS235914
do
  # get reads from EBI
  wget -q -O ${i}_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR369/${i:0:6}/${i}_1.fastq.gz
  wget -q -O ${i}_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR369/${i:0:6}/${i}_2.fastq.gz

  # map reads to reference with BWA & output as BAM file
  bwa mem -t 4 ref.fa ${i}_1.fastq.gz ${i}_2.fastq.gz | samtools view -b - > ${i}.aln.bam

  # sort BAM file
  samtools sort -@ 4 -O bam -T tmp ${i}.aln.bam > ${i}.aln.sort.bam

  # index BAM file
  samtools index ${i}.aln.sort.bam

  # generate genotype likelihoods with `samtools mpileup` and base calls with `bcftools call`, output as VCF
  samtools mpileup -u -f ref.fa ${i}.aln.sort.bam | bcftools call -c -O z -o ${i}.aln.sort.vcf.gz

  # index VCF files
  bcftools index ${i}.aln.sort.vcf.gz

  # generate consensus sequences
  bcftools consensus -f ref.fa ${i}.aln.sort.vcf.gz | sed "s/>.*/>$i/g" > ${i}.consensus.fa

done


# concatenate consensus files into one multi-fasta files
cat /escratch4/s_11/s_11_Aug_17/project/*.consensus.fa > /escratch4/s_11/s_11_Aug_17/project/allconsensus.fasta


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
