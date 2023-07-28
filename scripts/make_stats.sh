#!/bin/sh

SRC=/Align_RF/scripts  # adjust this for your environment

module --force purge
module load StdEnv 
module load GCC/11.3.0
module load Python/3.10.4-GCCcore-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

# Put caballus genome on the right
# so is_primary means is_caballus

date
echo 'Extract asinus read stats from its two BAM files'
python ${SRC}/bam_two_targets.py \
       --diff \
       map_asinus_to_asinus/Sorted.bam \
       map_asinus_to_caballus/Sorted.bam \
       > asinus_read_stats.csv
echo -n $?
echo " exit status"

date
echo 'Extract caballus read stats from its two BAM files'
python ${SRC}/bam_two_targets.py \
       --diff \
       map_caballus_to_asinus/Sorted.bam \
       map_caballus_to_caballus/Sorted.bam \
       > caballus_read_stats.csv
echo -n $?
echo " exit status"

date
gzip -v *stats.csv

date
ls -l
