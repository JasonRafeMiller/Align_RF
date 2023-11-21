#!/bin/sh

# The reduced.bam files contain the first several million lines
# of the Sorted.bam files ... more than enough for training models.

SRC=..

module --force purge
module load StdEnv 
module load GCC/11.3.0
module load Python/3.10.4-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

# Put oleracea genome on the right
# so is_primary means is_oleracea

date
echo 'Extract rapa read stats from its two BAM files'
python ${SRC}/bam_two_targets.py \
       --diff \
       --maxq F \
       map_rapa_to_rapa/Sorted.bam \
       map_rapa_to_oleracea/Sorted.bam \
       > unfiltered.rapa_read_stats.csv
echo -n $?
echo " exit status"

date
echo 'Extract oleracea read stats from its two BAM files'
python ${SRC}/bam_two_targets.py \
       --diff \
       --maxq F \
       map_oleracea_to_rapa/Sorted.bam \
       map_oleracea_to_oleracea/Sorted.bam \
       > unfiltered.oleracea_read_stats.csv
echo -n $?
echo " exit status"

date
gzip -v *stats.csv

date
ls -l
