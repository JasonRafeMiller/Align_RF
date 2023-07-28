#!/bin/sh

# Samtools cannot handle all references having the same ID.

gunzip -c ../Transcriptomes/Equus_asinus.ASM1607732v2.cdna.all.fa.gz | sed 's/^>/>AS_/' > diploid.fasta

gunzip -c ../Transcriptomes/Equus_caballus.EquCab3.0.cdna.all.fa.gz | sed 's/^>/>CA_/' >> diploid.fasta

