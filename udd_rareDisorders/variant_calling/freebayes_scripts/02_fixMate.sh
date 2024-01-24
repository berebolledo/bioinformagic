#! /bin/bash

bam_file=${1}.bam
bam_file_fixedmate=fixed.${1}.bam
mkdir -p temp

picard \
	FixMateInformation -Xmx8g \
	MAX_RECORDS_IN_RAM=2000000 \
	VALIDATION_STRINGENCY=SILENT \
	ADD_MATE_CIGAR=true \
	ASSUME_SORTED=false \
	SORT_ORDER=coordinate \
	TMP_DIR=temp \
	I=$bam_file \
	O=$bam_file_fixedmate
