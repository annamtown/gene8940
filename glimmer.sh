#!/bin/bash

cd /home/student/pbio4550/s_74/vcholerae_genome/Vcholera_spades

export /usr/local/glimmer/latest/bin/glimmer3:$PATH

glimmer3 scaffolds.fasta scaffolds.icm scaffolds

glimmer3 contigs.fasta contigs.com contigs
