#!/usr/bin/env bash
# dry_run.sh - runs the WHOLE pipeline (run_all.sh, all 7 steps) with pretend tools (test/fake_tool.sh) on empty files.
# Why: proves that every step finds the files the step before it wrote, resume works, and the final report is built,
#      without installing any bioinformatics tool or downloading data. Real results need real tools (see README).
set -euo pipefail                                                   # stop at the first error
R=$(cd "$(dirname "$0")/.." && pwd); W=$(mktemp -d); [ -n "${KEEP:-}" ] || trap 'rm -rf "$W"' EXIT   # repo folder; temporary work folder, deleted at the end
mkdir -p "$W/bin" "$W/fq"                                           # folder for the pretend tools and for fake FASTQ files
echo "bwa-mem2 samtools bcftools tabix gatk fastp OptiTypePipeline.py run-t1k vep whatshap docker" | tr ' ' '\n' | while read -r t; do
  ln -s "$R/test/fake_tool.sh" "$W/bin/$t"                          # every tool name points to the same pretend script
done
printf 'chr1\t0\t1000000\nchr1\t500000\t1500000\nchr2\t0\t500000\n' > "$W/targets.bed"   # tiny exome BED (overlap tests the merge): 2.0 Mb
touch "$W"/fq/{T,N}_R{1,2}.fq.gz                                    # empty tumor and normal FASTQ files
printf 'patient,sample_type,fastq_1,fastq_2,hla_alleles\nPT1,tumor,%s,%s,\nPT1,normal,%s,%s,\n' \
  "$W/fq/T_R1.fq.gz" "$W/fq/T_R2.fq.gz" "$W/fq/N_R1.fq.gz" "$W/fq/N_R2.fq.gz" > "$W/samplesheet.csv"   # one tumor/normal pair
export PATH="$W/bin:$PATH" OUT="$W/results" TARGETS="$W/targets.bed" REF="$W/ref.fa" THREADS=2   # pretend tools first in PATH; settings
bash "$R/run_all.sh" "$W/samplesheet.csv"                           # first run: every step executes
bash "$R/run_all.sh" "$W/samplesheet.csv" | grep -c '\[skip\]' | grep -qx 8 || { echo "FAIL: resume did not skip all 8 steps"; exit 1; }   # second run skips all
grep -q 'TMB:</b> 1.50' "$OUT/PT1/report/PT1.neoantigen_report.html" || { echo "FAIL: TMB not in report"; exit 1; }   # 3 mutations / 2.0 Mb
[ "$(wc -l < "$OUT/PT1/report/PT1.neoantigens_ranked.tsv")" -eq 6 ] || { echo "FAIL: expected 5 ranked candidates"; exit 1; }   # header + 5 rows
echo "PASS dry run: all steps ran, resume works, report built"      # all checks passed
