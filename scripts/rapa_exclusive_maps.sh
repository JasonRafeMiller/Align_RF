#!/bin/sh

cd map_rapa_to_rapa
echo "Warning: awk substr(13) assumes IDs like SRR21735970.51"
samtools view Filtered.bam | cut -f 1 | head -n 5
cd ..

echo rapa mapped reads
cd map_rapa_to_rapa
samtools view Filtered.bam | awk '{print sprintf("%09d",substr($1,13))}' | uniq > mapped.ids
cd ..
echo oleracea mapped reads
cd map_rapa_to_oleracea
samtools view Filtered.bam | awk '{print sprintf("%09d",substr($1,13))}' | uniq > mapped.ids
cd ..

echo mapped exclusively to rapa
comm -23 map_rapa_to_rapa/mapped.ids map_rapa_to_oleracea/mapped.ids > rapa_to_rapa.ids

echo mapped exclusively to oleracea
comm -13 map_rapa_to_rapa/mapped.ids map_rapa_to_oleracea/mapped.ids > rapa_to_oleracea.ids

echo mapped to both
comm -12 map_rapa_to_rapa/mapped.ids map_rapa_to_oleracea/mapped.ids > rapa_to_both.ids

wc -l map*/*.ids
wc -l *.ids
echo done
