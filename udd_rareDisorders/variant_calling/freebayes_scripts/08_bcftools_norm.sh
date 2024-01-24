#! /bin/bash

refdir="/home/administrador/nfs_icim/1kGP-trio-data/1kGP-trio-ref"
reference_fasta_file=${refdir}/GRCh38_full_analysis_set_plus_decoy_hla.fa

vcf=${1}
out=norm.${vcf}

bcftools norm                  \
	-f ${reference_fasta_file} \
	-c s                       \
	-m -both                   \
	-O z                       \
	-o ${out}                  \
	${vcf}
