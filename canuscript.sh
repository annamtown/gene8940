canu \
 -p ecoli -d ecoli-pacbio \
 genomeSize=4.8m \
 -pacbio-raw pacbio.fastq

 curl -L -o mix.tar.gz http://gembox.cbcb.umd.edu/mhap/raw/ecoliP6Oxford.tar.gz
tar xvzf mix.tar.gz

canu \
 -p ecoli -d ecoli-mix \
 genomeSize=4.8m \
 -pacbio-raw pacbio.part?.fastq.gz \
 -nanopore-raw oxford.fasta.gz

 canu -correct \
  -p ecoli -d ecoli \
  genomeSize=4.8m \
  -pacbio-raw  pacbio.fastq

  canu -trim \
  -p ecoli -d ecoli \
  genomeSize=4.8m \
  -pacbio-corrected ecoli/ecoli.correctedReads.fasta.gz

  canu -assemble \
  -p ecoli -d ecoli-erate-0.039 \
  genomeSize=4.8m \
  correctedErrorRate=0.039 \
  -pacbio-corrected ecoli/ecoli.trimmedReads.fasta.gz

canu -assemble \
  -p ecoli -d ecoli-erate-0.075 \
  genomeSize=4.8m \
  correctedErrorRate=0.075 \
  -pacbio-corrected ecoli/ecoli.trimmedReads.fasta.gz
