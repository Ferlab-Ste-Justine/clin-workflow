#!/usr/bin/env bash
cat ${1} | java -jar ~/bin/snpEff/SnpSift.jar extractFields -s "," -e " " - CHROM POS ID REF ALT QUAL FILTER DP MQ MQRankSum ReadPosRankSum LOD FractionInformativeReads\
 SNP MNP INS DEL MIXED HOM HET "GEN[*].GT" "GEN[*].AD" "GEN[*].AF" "GEN[*].F1R2" "GEN[*].F2R1" "GEN[*].DP" "GEN[*].SB" "GEN[*].MB" CSQ
