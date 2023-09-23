#! /bin/bash

#$ -N pipeline
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -pe smp 8

#set -e
#set -u
#set -o pipefail

if [ $HOSTNAME == 'sofia.udd.cl' ] || [[ $HOSTNAME == compute*-1-*.local ]]
then
    genomes="/hpcudd/ICIM/shared/genomes"
    bundle="/hpcudd/ICIM/shared/gatk-bundle"
    annDir="/hpcudd/ICIM/boris/annotation"
    export PATH="/hpcudd/home/boris/miniconda3/bin:$PATH"
elif [ $HOSTNAME == 'mendel' ]
then
    genomes="/storage/shared/references"
    bundle="/storage/shared/gatk-bundle"
else
    echo "Unrecognized host $HOSTNAME"
    echo "can't locate genome references"
    exit 1
fi


index="${genomes}/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"
refdata="${genomes}/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"
dbsnp="${bundle}/b37/dbsnp_138.b37.vcf.gz"
indels="${bundle}/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"

bash decipherd_01_alignment.sh
bash decipherd_02_GATK.sh
bash decipherd_03_annovar.sh