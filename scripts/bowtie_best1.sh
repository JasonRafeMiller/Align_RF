#!/bin/sh

module --force purge
module load StdEnv 
module load GCC/11.3.0
module load Bowtie2/2.4.5-GCC-11.3.0

if [ "$#" -ne  1 ]
then
  echo "Argument required: a name for the index"
  exit 1
fi

NAME=$1
echo NAME ${NAME}

date
ls -l unzip_this
gunzip -v -c unzip_this > ${NAME}.fa

date
bowtie2-build --quiet --threads 4 ${NAME}.fa ${NAME}
echo -n $?
echo " exit status"
date

rm -v ${NAME}.fa
echo "done"
[login-4:bowtie]$ ca bowtie_best1.sh 
-bash: ca: command not found
[login-4:bowtie]$ cat bowtie_best1.sh 
#!/bin/sh
#SBATCH --account=${ACCOUNT}
#SBATCH --job-name=bowtie
#SBATCH --time=08:00:00   # Bowtie takes over 4 hr to find best 2
#SBATCH --mem-per-cpu=4G  # 16 GB total
#SBATCH --cpus-per-task=4  # 4 cpu is optimal for 4 threads
## source /cluster/bin/jobsetup   ## abel only
set -o errexit # exit on errors
# Our python will generate smaller Aligned.bam which we keep.
#savefile *.bam
#savefile *.db
#savefile *.log
#savefile *.SN

echo 'Launch this in the sub-directory with the index and fastq files.'
# sbatch --account=${ACCOUNT} ../bowtie_best1.sh

echo MODULES
module --force purge
module load StdEnv 
module load GCC/11.3.0
module load Bowtie2/2.4.5-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module list

echo
echo LD_LIBRARY_PATH $LD_LIBRARY_PATH
echo

date
INITIALDIR=`pwd`
echo INITIALDIR ${INITIALDIR}
echo

THREADS=4
echo THREADS ${THREADS}

# expect bowtie index file like lyrata.1.bt2 
IDX=*1.bt2
TARGET=` echo $IDX | sed 's/\..*//g' `
echo "TARGET ${TARGET}"

R1=*_R1_*.fq.gz
R2=*_R2_*.fq.gz
echo R1 ${R1}
echo R2 ${R2}
echo

ls -l
echo

date
echo run Bowtie
# assumes *.fasta and *.bt2 have same first name
# we want one alignment per read, so do not use k2
bowtie2  \
    --no-unal \
    --no-mixed \
    --no-discordant \
    --sensitive \
    --end-to-end \
    --threads ${THREADS} \
    -x ${TARGET} \
    -1 ${R1} -2 ${R2} \
    -S Aligned.sam 2>> run_bowtie.log
echo -n $?
echo " exit status"
date

echo compress sam to bam
samtools view -b -o Primary.bam Aligned.sam

echo sort
samtools sort -n -T . -@ 4 -o Sorted.bam Primary.bam

echo stats
samtools flagstat Primary.bam > samtools.flagstat
samtools stats Primary.bam | grep '^SN' | cut -f 2- > samtools.stats.SN

date
ls -l
echo DONE
