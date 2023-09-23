#! /bin/bash

#$ -N bcfmerge
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

file_list=${1}
output=${2}

bcftools merge           \
    -l ${file_list}      \
    -o ${output}.vcf.gz  \
    -O z -m all

tabix -p vcf ${output}.vcf.gz


