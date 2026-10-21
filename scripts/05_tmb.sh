#!/usr/bin/env bash
# 05_tmb.sh - tumor mutational burden: nonsynonymous somatic mutations per Mb of the exome target.
# Usage: bash 05_tmb.sh <PATIENT>
# Why: high TMB means more mutations, so more chances for neoantigens; it is also a biomarker for immunotherapy response.
set -euo pipefail                                                   # stop at the first error
source "$(dirname "$0")/../config.sh"                               # load paths and settings
P=$1; D=$OUT/$P/tmb; mkdir -p "$D"                                  # patient and output folder
# exome size in Mb: sort the targets, merge overlapping intervals (so no base is counted twice), add up the lengths
MB=$(sort -k1,1 -k2,2n "$TARGETS" | awk '$1 != c || $2 > e { s += e - b; c = $1; b = $2; e = $3; next } $3 > e { e = $3 } END { printf "%.4f", (s + e - b) / 1e6 }')
# one row per variant: worst VEP consequence + read counts; 'view -s T,N' fixes the column order to tumor then normal
bcftools view -s "${P}_T,${P}_N" "$OUT/$P/pvacseq_input/$P.somatic.vep.vcf.gz" \
  | bcftools +split-vep -s worst -f '%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%SYMBOL[\t%AD]\n' > "$D/variants.tsv"
awk -F'\t' -v mb="$MB" -v depth="$TMB_DEPTH" -v alt="$TMB_ALT" -v vaf="$TMB_VAF" -v nvaf="$TMB_NVAF" -v high="$TMB_HIGH" \
  -v list="$D/$P.tmb_variants.tsv" -f "$(dirname "$0")/../bin/tmb.awk" "$D/variants.tsv" > "$D/$P.tmb.tsv"   # count and classify
cat "$D/$P.tmb.tsv"                                                 # show the result in the log
