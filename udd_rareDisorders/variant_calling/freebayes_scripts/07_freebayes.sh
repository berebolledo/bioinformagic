#! /bin/bash

refdir="/home/administrador/nfs_icim/1kGP-trio-data/1kGP-trio-ref"
reference_fasta_file=${refdir}/GRCh38_full_analysis_set_plus_decoy_hla.fa



bam=${1}
chrom=chr${2}
output=${chrom}_${3}

targets=targets.txt

#   -m --min-mapping-quality Q
#                   Exclude alignments from analysis if they have a mapping
#                   quality less than Q.  default: 1
#
#   -q --min-base-quality Q
#                   Exclude alleles from analysis if their supporting base
#                   quality is less than Q.  default: 0
#
#   -z --read-max-mismatch-fraction N
#                   Exclude reads with more than N [0,1] fraction of mismatches where
#                   each mismatch has base quality >= mismatch-base-quality-threshold
#                   default: 1.0
#
#   -C --min-alternate-count N
#                   Require at least this count of observations supporting
#                   an alternate allele within a single individual in order
#                   to evaluate the position.  default: 2
#
#   --min-coverage N
#                   Require at least this coverage to process a site. default: 0
#
#   -N --exclude-unobserved-genotypes
#                   Skip sample genotypings for which the sample has no supporting reads.
#
#   -G --min-alternate-total N
#                   Require at least this count of observations supporting
#                   an alternate allele within the total population in order
#                   to use the allele in analysis.  default: 1


freebayes \
	-m 30             \
	-q 30             \
	-z 0.5            \
	-C 3              \
	--min-coverage 20 \
	-G 3              \
	-f ${reference_fasta_file} \
	-r ${chrom} ${bam} |bgzip -c > ${output}.vcf.gz
tabix -p vcf ${output}.vcf.gz
