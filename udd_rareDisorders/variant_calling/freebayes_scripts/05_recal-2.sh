#! /bin/bash

refdir="/home/administrador/nfs_icim/1kGP-trio-data/1kGP-trio-ref"
reference_fasta_file=${refdir}/GRCh38_full_analysis_set_plus_decoy_hla.fa

bam_dedup_sorted=dedup.sorted.fixed.${1}.bam
recalibrated_bam=recal.dedup.sorted.fixed.${1}.bam
recal_data_table=recal_data.${1}

gatk3  \
	-T PrintReads -Xmx8g \
	-nct 8 \
	--disable_indel_quals \
	--preserve_qscores_less_than 6 \
	-SQQ 10 \
	-SQQ 20 \
	-SQQ 30 \
	-rf BadCigar \
	-R $reference_fasta_file \
	-o $recalibrated_bam \
	-I $bam_dedup_sorted \
	-BQSR $recal_data_table
