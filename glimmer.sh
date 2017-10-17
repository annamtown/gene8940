#!/bin/bash
#$ -o /home/student/pbio4550/s_74/vcholerae_genome
#$ -e /home/student/pbio4550/s_74/vcholerae_genome


cd /home/student/pbio4550/s_74/vcholerae_genome/Vcholera_spades

export /usr/local/glimmer/3.0.2/:$PATH

glimmer3 scaffolds.fasta scaffolds.output scaffolds.predict

glimmer3 contigs.fasta contigs.output contigs.predict
