#!/bin/bash

cd /escratch4/s_11/s_11_Aug_17/ecoli

export PATH=/usr/local/canu/1.4/Linux-amd64/bin:/usr/local/java/jdk1.8.0_74/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/gcc/5.3.0/lib64:${LD_LIBRARY_PATH}
export JAVA_HOME=/usr/local/java/jdk1.8.0_74

curl -L -o pacbio.fastq http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq

canu \
 -p ecoli -d ecoli-pacbio genomeSize=4.8m -pacbio-raw pacbio.fastq useGrid=false
