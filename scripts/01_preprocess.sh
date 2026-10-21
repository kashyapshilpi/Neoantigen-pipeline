#!/usr/bin/env bash
# 01_preprocess.sh - FASTQ -> trimmed reads -> aligned, duplicate-marked, recalibrated BAM.
# Usage: bash 01_preprocess.sh <SAMPLE_ID> <R1.fastq.gz> <R2.fastq.gz>
set -euo pipefail                                                   # stop at the first error, unset variable or failed pipe
source "$(dirname "$0")/../config.sh"                               # load paths and settings
S=$1; R1=$(realpath "$2"); R2=$(realpath "$3"); D=$OUT/$S/bam; mkdir -p "$D"   # sample id, absolute FASTQ paths, output folder

# remove adapters and low-quality bases: they would otherwise look like mismatches (false variants)
fastp -i "$R1" -I "$R2" -o "$D/R1.fq.gz" -O "$D/R2.fq.gz" -w "$THREADS" -j "$D/$S.fastp.json" -h "$D/$S.fastp.html"
# align to GRCh38; the read group (SM = sample name) is how Mutect2 tells tumor reads from normal reads; then sort by position
bwa-mem2 mem -t "$THREADS" -R "@RG\tID:$S\tSM:$S\tPL:ILLUMINA\tLB:$S" "$REF" "$D/R1.fq.gz" "$D/R2.fq.gz" | samtools sort -@ 4 -o "$D/$S.sorted.bam" -
# mark PCR duplicates so one DNA fragment is not counted as several independent reads supporting a variant
gatk MarkDuplicates -I "$D/$S.sorted.bam" -O "$D/$S.md.bam" -M "$D/$S.dup_metrics.txt" --CREATE_INDEX true
if [ -n "$KNOWN_SITES" ]; then                                      # BQSR only when known-variant VCFs are configured
  KS=(); while read -r v; do KS+=(--known-sites "$v"); done < <(tr ' ' '\n' <<< "$KNOWN_SITES")   # one --known-sites flag per VCF
  gatk BaseRecalibrator -I "$D/$S.md.bam" -R "$REF" "${KS[@]}" -O "$D/$S.recal.table"   # learn the sequencer's base-quality errors
  gatk ApplyBQSR -I "$D/$S.md.bam" -R "$REF" --bqsr-recal-file "$D/$S.recal.table" -O "$D/$S.final.bam"   # correct the qualities
else ln -sf "$S.md.bam" "$D/$S.final.bam"; ln -sf "$S.md.bai" "$D/$S.final.bai"; fi   # no known sites: use the duplicate-marked BAM
rm -f "$D/R1.fq.gz" "$D/R2.fq.gz" "$D/$S.sorted.bam"                # delete big intermediate files to save disk space
samtools flagstat "$D/$S.final.bam" > "$D/$S.flagstat.txt"          # quick QC: mapped %, properly paired %, duplicates %
