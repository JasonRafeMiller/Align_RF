#!/bin/sh

SRC=..

module --force purge
module load StdEnv 
module load GCC/11.3.0
module load Python/3.10.4-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

# Put oleracea genome on the right
# so is_primary means is_oleracea

chmod 444 map_rapa_to_rapa/Sorted.bam
chmod 444 map_rapa_to_oleracea/Sorted.bam
chmod 444 map_oleracea_to_rapa/Sorted.bam
chmod 444 map_oleracea_to_oleracea/Sorted.bam

## samtools flags: -f 2 "read mapped in proper pair", -F 256 not "not primary alignment"

date
cd map_rapa_to_rapa
pwd
samtools view -b -f 2 -F 256 -q 1 -o Filtered.bam Sorted.bam
cd ..
cd map_rapa_to_oleracea
pwd
samtools view -b -f 2 -F 256 -q 1 -o Filtered.bam Sorted.bam
cd ..
cd map_oleracea_to_rapa
pwd
samtools view -b -f 2 -F 256 -q 1 -o Filtered.bam Sorted.bam
cd ..
cd map_oleracea_to_oleracea
pwd
samtools view -b -f 2 -F 256 -q 1 -o Filtered.bam Sorted.bam
cd ..

date
echo 'Extract rapa read stats from its two BAM files'
python ${SRC}/bam_two_targets.py \
       --diff \
       --maxq F \
       map_rapa_to_rapa/Filtered.bam \
       map_rapa_to_oleracea/Filtered.bam \
       > filtered.rapa_read_stats.csv
echo -n $?
echo " exit status"

date
echo 'Extract oleracea read stats from its two BAM files'
python ${SRC}/bam_two_targets.py \
       --diff \
       --maxq F \
       map_oleracea_to_rapa/Filtered.bam \
       map_oleracea_to_oleracea/Filtered.bam \
       > filtered.oleracea_read_stats.csv
echo -n $?
echo " exit status"

date
gzip -v *stats.csv

date
ls -l
