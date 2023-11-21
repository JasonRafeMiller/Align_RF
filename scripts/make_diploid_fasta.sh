#!/bin/sh

# Samtools cannot handle all references having the same ID.

gunzip -c ../Transcripts/rapa/Brassica_rapa.Brapa_1.0.cdna.all.fa.gz | sed 's/^>/>RA_/' > diploid.fasta

gunzip -c ../Transcripts/oleracea/Brassica_oleracea.BOL.cdna.all.fa.gz  | sed 's/^>/>OL_/' >> diploid.fasta

