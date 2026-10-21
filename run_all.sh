#!/usr/bin/env bash
# run_all.sh - ONE command: tumor/normal WES FASTQ -> ranked MHC class I + II neoantigens for every patient.
# Usage: bash run_all.sh samplesheet.csv     (run the same command again to resume: finished steps are skipped)
set -uo pipefail                                                  # stop on unset variables / failed pipes (step errors are handled in step())
H=$(cd "$(dirname "$0")" && pwd); source "$H/config.sh"           # folder of this script, then load every path and setting
[ $# -eq 1 ] && [ -s "$1" ] || { sed -n 2,3p "$0"; exit 1; }      # no samplesheet given -> print the two usage lines above and stop
SHEET=$(realpath "$1"); mkdir -p "$OUT/logs" "$OUT/.done"         # absolute path (still valid after cd) + folders for logs and done-markers

step() {                                                          # step <name> <command...>: run once, log it, skip it next time
  local n=$1; shift                                               # first argument = step name, the rest = the command to run
  [ -f "$OUT/.done/$n" ] && { echo "[skip] $n"; return 0; }       # a done-marker means it already finished -> this makes the run resumable
  echo "[$(date '+%F %T')] START $n"                              # progress line on screen
  if "$@" < /dev/null > "$OUT/logs/$n.log" 2>&1; then             # run it; </dev/null stops tools from eating the samplesheet lines
    touch "$OUT/.done/$n"; echo "[$(date '+%F %T')] DONE  $n"     # success -> write the done-marker
  else echo "FAILED $n, last lines of $OUT/logs/$n.log:"; tail -20 "$OUT/logs/$n.log"; exit 1; fi   # failure -> show the error, stop
}

echo "bwa-mem2 samtools bcftools tabix gatk fastp OptiTypePipeline.py run-t1k vep whatshap docker" | tr ' ' '\n' | while read -r t; do
  command -v "$t" > /dev/null || { echo "MISSING tool: $t (see README: Install)"; exit 1; }   # check tools first, not after hours of work
done || exit 1                                                    # the loop runs in a subshell, so stop the script here if a tool was missing

tail -n +2 "$SHEET" | tr -d '\r' | while IFS=, read -r p type r1 r2 _; do    # step 1 for every FASTQ pair (tail skips the header row)
  case $type in tumor) s=${p}_T ;; normal) s=${p}_N ;; *) echo "sample_type must be tumor or normal, got '$type'"; exit 1 ;; esac   # sample id
  step "$s.preprocess" bash "$H/scripts/01_preprocess.sh" "$s" "$r1" "$r2"   # FASTQ -> analysis-ready BAM
done || exit 1                                                    # stop if any sample failed

tail -n +2 "$SHEET" | tr -d '\r' | cut -d, -f1 | sort -u | while read -r p; do   # steps 2-7 once per patient (needs both BAMs)
  hla=$(tail -n +2 "$SHEET" | tr -d '\r' | awk -F, -v p="$p" '$1 == p && $5 != "" {print $5; exit}')   # clinical HLA typing, if given
  step "$p.variants" bash "$H/scripts/02_variants.sh" "$p"         # somatic (Mutect2) + germline (HaplotypeCaller) variants
  step "$p.hla"      bash "$H/scripts/03_hla.sh" "$p" "$hla"       # HLA class I + II from the normal sample
  step "$p.annotate" bash "$H/scripts/04_annotate_phase.sh" "$p"   # VEP annotation + phased proximal variants = pVACseq input
  step "$p.tmb"      bash "$H/scripts/05_tmb.sh" "$p"              # tumor mutational burden
  step "$p.pvacseq"  bash "$H/scripts/06_pvacseq.sh" "$p"          # MHC class I + II binding prediction
  step "$p.rank"     bash "$H/scripts/07_rank_report.sh" "$p"      # one ranked table + HTML report
done || exit 1                                                    # stop if any patient failed
echo "Finished. Reports: $OUT/<patient>/report/<patient>.neoantigen_report.html"   # where to look
