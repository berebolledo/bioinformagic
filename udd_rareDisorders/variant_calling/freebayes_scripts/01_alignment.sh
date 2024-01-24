#! /bin/bash

refdir="/home/administrador/nfs_icim/1kGP-trio-data/1kGP-trio-ref"
reference_fasta_file=${refdir}/GRCh38_full_analysis_set_plus_decoy_hla.fa
rg_string="@RG\tID:${1}\tLB:22qpairs\tSM:${1}\tPL:ILLUMINA"
fastq1=${1}_1.fq.gz
fastq2=${1}_2.fq.gz
bam_file=${1}.bam

bwa mem -Y \
	-K 100000000 \
	-t 8 \
	-R $rg_string \
	$reference_fasta_file \
	$fastq1 \
	$fastq2 | samtools view -Shb -o $bam_file -
