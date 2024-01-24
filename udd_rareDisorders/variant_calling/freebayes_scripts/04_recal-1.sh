#! /bin/bash 

refdir="/home/administrador/nfs_icim/1kGP-trio-data/1kGP-trio-ref"
reference_fasta_file=${refdir}/GRCh38_full_analysis_set_plus_decoy_hla.fa

known_indels_from_mills_1000genomes=${refdir}/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz
known_snps_from_dbSNP142=${refdir}/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz
#known_indels=${refdir}/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz
#-knownSites $known_indels

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
