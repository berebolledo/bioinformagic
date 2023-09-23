#! /bin/bash

#$ -N fq2annvcf
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

while getopts '1:2:i:l:h' ARGS; do
        case "$ARGS" in
        1)
          read1="$OPTARG"
          ;;
        2)
          read2="$OPTARG"
          ;;
        i)
          RGID="$OPTARG"
          ;;
        l)
          RGLB="$OPTARG"
          ;;
        h)
          echo "script usage: $(basename $0) [-1 read1.fq] [-2 read2.fq] [-i readgroup ID] [-l readgroup LIB]" >&2
          exit 0
          ;;
        ?)
          echo "script usage: $(basename $0) [-1 read1.fq] [-2 read2.fq] [-i readgroup ID] [-l readgroup LIB]" >&2
          exit 1
          ;;
    esac
done
    
shift "$(($OPTIND - 1))"
mkdir -p ${RGID}_tmpdir

index="${genomes}/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"
refdata="${genomes}/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"
dbsnp="${bundle}/b37/dbsnp_138.b37.vcf.gz"
indels="${bundle}/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"


if [ $exit_left -eq 0 ] && [ -s left-aligned.bqsr.markDups.sorted.${RGID}.bam ]
then
    gatk -Xms4g -Xmx8g     \
        -T HaplotypeCaller \
        -R $genome         \
        -I left-aligned.bqsr.markDups.sorted.${RGID}.bam \
        --genotyping_mode DISCOVERY      \
        --emitRefConfidence GVCF         \
        --variant_index_type LINEAR      \
        --variant_index_parameter 128000 \
        -o ${RGID}.gvcf
    exit_gvcf=$?
fi

if [ $exit_gvcf -eq 0 ] && [ -s ${RGID}.gvcf ]
then
    gatk -Xms4g -Xmx8g   \
        -T GenotypeGVCFs \
        -R $genome       \
        -V ${RGID}.gvcf  \
        -o ${RGID}.raw.vcf
    exit_vcf=$?
fi

