#!/bin/sh

# Out of first 2 million read pairs mapped, how may mapped to their parent?
echo 'count reads (not pairs)'
echo 'remember to divide by two'

echo 'asinus'
samtools view map_asinus_to_diploid/Primary.bam | head -n 2000000 | cut -f 3 | grep -c '^AS'
echo 'caballus'
samtools view map_caballus_to_diploid/Primary.bam | head -n 2000000 | cut -f 3 | grep -c '^CA'
