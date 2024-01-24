#! /bin/bash

bam_file_fixedmate=fixed.${1}.bam
dedup_metrics=dedup.metrics.${1}
bam_dedup_sorted=dedup.sorted.fixed.${1}.bam 
mkdir -p temp

picard \
	MarkDuplicates -Xmx8g \
	MAX_RECORDS_IN_RAM=2000000 \
	VALIDATION_STRINGENCY=SILENT \
	TMP_DIR=temp \
	ASSUME_SORTED=true \
	CREATE_INDEX=true \
	M=$dedup_metrics \
	I=$bam_file_fixedmate \
	O=$bam_dedup_sorted 
