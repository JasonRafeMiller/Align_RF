#!/bin/sh


function count_all() {
    pwd
    echo no filter

    echo rapa numerator
    samtools view map_rapa_to_diploid/Sorted.bam | cut -f 3 | grep '^RA_' | wc -l
    echo rapa denominator
    samtools view map_rapa_to_diploid/Sorted.bam | wc -l

    echo oleracea numerator
    samtools view map_oleracea_to_diploid/Sorted.bam | cut -f 3 | grep '^OL_' | wc -l
    echo oleracea denominator
    samtools view map_oleracea_to_diploid/Sorted.bam | wc -l

    echo pair filter

    echo rapa numerator
    samtools view -f 66 -F 256 map_rapa_to_diploid/Sorted.bam | cut -f 3 | grep '^RA_' | wc -l
    echo rapa denominator
    samtools view -f 66 -F 256 map_rapa_to_diploid/Sorted.bam | wc -l

    echo oleracea numerator
    samtools view -f 66 -F 256 map_oleracea_to_diploid/Sorted.bam | cut -f 3 | grep '^OL_' | wc -l
    echo oleracea denominator
    samtools view -f 66 -F 256 map_oleracea_to_diploid/Sorted.bam | wc -l

    echo plus mapq filter

    echo rapa numerator
    samtools view -f 66 -F 256 -q 1 map_rapa_to_diploid/Sorted.bam | cut -f 3 | grep '^RA_' | wc -l
    echo rapa denominator
    samtools view -f 66 -F 256 -q 1 map_rapa_to_diploid/Sorted.bam | wc -l

    echo oleracea numerator
    samtools view -f 66 -F 256 -q 1 map_oleracea_to_diploid/Sorted.bam | cut -f 3  | grep '^OL_' | wc -l
    echo oleracea denominator
    samtools view -f 66 -F 256 -q 1 map_oleracea_to_diploid/Sorted.bam | wc -l
}

count_all
