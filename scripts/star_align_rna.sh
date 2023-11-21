#!/bin/sh
#SBATCH --account=${ACCOUNT}
#SBATCH --job-name=star
#SBATCH --time=08:00:00   
#SBATCH --mem-per-cpu=8G  # STAR out of memory using 4G/cpu=16GB total
#SBATCH --cpus-per-task=4  # 4 cpu is optimal for 4 threads
set -o errexit # exit on errors
#savefile *.bam
#savefile *.db
#savefile *.log
#savefile *.SN

# How to submit this script
# sbatch --account=${ACCOUNT} ../star_align.sh

echo MODULES
module --force purge
module load StdEnv 

module load GCC/11.3.0
module load STAR/2.7.10b-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module list

echo
echo LD_LIBRARY_PATH $LD_LIBRARY_PATH
echo

date
INITIALDIR=`pwd`
echo INITIALDIR ${INITIALDIR}
echo THREADS ${THREADS}

# Location of STAR index file called Genome
GENOMEDIR='..'
echo GENOMEDIR $GENOMEDIR

# One thread is slower but read order is preserved in every BAM.
THREADS=4
echo 

R1=*_R1_*.fq.gz
R2=*_R2_*.fq.gz
echo R1 ${R1}
echo R2 ${R2}
echo

ls -l
echo

date
echo run STAR
STAR --runThreadN ${THREADS} \
     --genomeDir ${GENOMEDIR} \
     --readFilesIn ${R1} ${R2} \
     --outSAMattributes NH AS nM NM MD \
     --outSAMtype BAM Unsorted \
     --alignIntronMin 100000 \
     --alignIntronMax 0 \
     --alignEndsType EndToEnd \
     --readFilesCommand gunzip -c 
echo -n $?
echo " exit status"
date

#echo compress sam to bam
#samtools view -b -o Primary.bam Aligned.out.sam
mv Aligned.out.bam Primary.bam

#echo sort
#samtools sort -n -T . -@ 4 -o Sorted.bam Primary.bam

echo stats
samtools flagstat Primary.bam > samtools.flagstat
samtools stats Primary.bam | grep '^SN' | cut -f 2- > samtools.stats.SN

date
ls -l
echo DONE
