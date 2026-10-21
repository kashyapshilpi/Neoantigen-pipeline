#!/usr/bin/env bash
# 03_hla.sh - HLA class I (OptiType) and class II (T1K) typing from the NORMAL (germline) sample.
# Usage: bash 03_hla.sh <PATIENT> ["clinical alleles; separated by ;"]
# Why the normal: HLA genes are inherited; the tumor can lose HLA alleles, so the germline sample gives the true type.
set -euo pipefail                                                   # stop at the first error
source "$(dirname "$0")/../config.sh"                               # load paths and settings
P=$1; N=$OUT/${P}_N/bam/${P}_N.final.bam; D=$OUT/$P/hla; mkdir -p "$D"   # patient, normal BAM, output folder

# keep only reads that can come from HLA genes: mapped reads in the MHC region + all unmapped reads
# (HLA genes are so variable that many of their reads fail to align to the reference and end up unmapped)
samtools merge -f -o "$D/hla.bam" <(samtools view -b -F 4 "$N" "$HLA_REGION") <(samtools view -b -f 4 "$N")
samtools collate -Ou "$D/hla.bam" "$D/tmp" | samtools fastq -1 "$D/R1.fq" -2 "$D/R2.fq" -0 /dev/null -s /dev/null -n   # back to paired FASTQ
OptiTypePipeline.py -i "$D/R1.fq" "$D/R2.fq" --dna -o "$D/optitype" -p "$P"   # class I (HLA-A, -B, -C): accurate from DNA reads
run-t1k -1 "$D/R1.fq" -2 "$D/R2.fq" --preset hla-wgs -f "$T1K_REF" -o "$P" --od "$D/t1k" -t "$THREADS"   # class II (DR, DQ, DP) + class I backup
# merge both callers into the allele format pVACseq expects; clinical alleles (if given) replace the computed ones of the same class
awk -F'\t' -v override="${2:-}" -f "$(dirname "$0")/../bin/hla_merge.awk" "$D/optitype/${P}_result.tsv" "$D/t1k/${P}_genotype.tsv" > "$D/$P.hla.txt"
[ -s "$D/$P.hla.txt" ] || { echo "No HLA alleles found for $P"; exit 1; }   # without alleles pVACseq cannot run
rm -f "$D/R1.fq" "$D/R2.fq" "$D/hla.bam"; cat "$D/$P.hla.txt"      # remove temporary reads; print the result into the log
