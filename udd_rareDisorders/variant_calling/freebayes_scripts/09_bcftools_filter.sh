#!/bin/bash

vcf=${1}

bcftools filter                                                             \
	-i "QUAL > 30 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
	-O z                                                                    \
	-o filt.${vcf}                                                          \
	${vcf}

tabix -p vcf filt.${vcf}

