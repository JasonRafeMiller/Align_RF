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
