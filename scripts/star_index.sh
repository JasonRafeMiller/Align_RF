#!/bin/sh

echo module
module --force purge
module load StdEnv 
module load GCC/11.3.0
module load STAR/2.7.10b-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module list

date
echo star
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --limitGenomeGenerateRAM 70079829600 \
     --genomeSAindexNbases 11 \
     --genomeDir . \
     --genomeFastaFiles diploid.fasta
echo -n $?
echo " exit status"
date


