#!/bin/bash
#$ -pe thread 4

cd /escratch4/s_11/s_11_Aug_17/project

#path for reference assembly executables
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

for i in ERR369378 ERR369348 ERR338264 ERR338265 ERR338272 ERR338273 ERR338280 ERR338281 ERR338288 ERR338255
do
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
  bcftools consensus -f ref.fa ${i}.aln.sort.vcf.gz > ${i}.consensus.fa

  #ANALYSIS OF ASSEMBLIES

  # run QUAST 3.1 on reference-based assemblies using reference
  python2.7 /usr/local/quast/3.1/quast.py -o /escratch4/s_11/s_11_Aug_17/project/quast_output -R /escratch4/s_11/s_11_Aug_17/project/ref.fa /escratch4/s_11/s_11_Aug_17/project/${i}.consensus.fa

  # make mummerplots for reference-based assemblies using Ensembl MG1655 as reference
  nucmer -o /escratch4/s_11/s_11_Aug_17/project/ref.fa /escratch4/s_11/s_11_Aug_17/project/${i}.consensus.fa -p outputref_${i}
  delta-filter -1 outputref_${i}.delta > outputref_${i}.1delta
  mummerplot --size large -fat --color -f --png outputref_${i}.1delta -p outputref_${i}

  # generate Prokka genome annotations for reference-based assemblies using Ensembl MG1655 as reference
  prokka /escratch4/s_11/s_11_Aug_17/project/${i}.consensus.fa --outdir prokka_${i}

done
