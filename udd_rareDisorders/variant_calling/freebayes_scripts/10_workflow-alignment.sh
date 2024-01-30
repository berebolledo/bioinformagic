#! /bin/bash

# 01_alignment.sh

refdir="/home/administrador/nfs_icim/1kGP-trio-data/1kGP-trio-ref"
reference_fasta_file=${refdir}/GRCh38_full_analysis_set_plus_decoy_hla.fa
rg_string="@RG\tID:${1}\tLB:mylibrary\tSM:${1}\tPL:ILLUMINA"
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

# 02_fixMate.sh

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

# 03_dedup.sh

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

# 04_recal-1.sh

known_indels_from_mills_1000genomes=${refdir}/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz
known_snps_from_dbSNP142=${refdir}/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz

for i in $(seq 1 22);do echo chr$i;done > autosomes.intervals

bam_dedup_sorted=dedup.sorted.fixed.${1}.bam
autosomes=autosomes.intervals
recal_data_table=recal_data.${1}
mkdir -p temp

gatk3  \
	-T BaseRecalibrator -Xmx8g \
	--downsample_to_fraction 0.1 \
	-nct 8 \
	--preserve_qscores_less_than 6 \
	-L $autosomes \
	-R $reference_fasta_file \
	-o $recal_data_table \
	-I $bam_dedup_sorted \
	--knownSites $known_snps_from_dbSNP142 \
	--knownSites $known_indels_from_mills_1000genomes


# 05_recal-2.sh

recalibrated_bam=recal.dedup.sorted.fixed.${1}.bam

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


