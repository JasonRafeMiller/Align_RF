#!/bin/sh

echo make hisat index
echo assume module load
date
hisat2-build --threads 4 diploid.fasta diploid
echo -n $?
echo " exit status"
date

