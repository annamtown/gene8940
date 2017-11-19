#!/bin/bash

cd /escratch4/s_11/s_11_Aug_17/samplestest

for i in ERR369378
do
  wget -q -O ${i}_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR369/${i:0:6}/${i}_1.fastq.gz
  wget -q -O ${i}_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR369/${i:0:6}/${i}_2.fastq.gz
done
