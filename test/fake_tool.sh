#!/usr/bin/env bash
# fake_tool.sh - pretend version of every bioinformatics tool, used only by test/dry_run.sh.
# It is linked under each tool's name; it creates the files the real tool would create (empty, or copied from test/data),
# so the whole pipeline can be checked end to end - every step, file name and hand-over - in seconds and without real data.
set -euo pipefail                                                   # stop at the first error
T=$(basename "$0"); DATA=$(dirname "$(readlink -f "$0")")/data   # which tool we are pretending to be; folder with test inputs
out() { mkdir -p "$(dirname "$1")"; touch "$1"; }                   # create an (empty) output file and its folder
a=("$@")                                                            # all arguments (positions below follow the calls in scripts/)
case $T in                                                          # tools whose outputs are not named by an -o style option
  OptiTypePipeline.py) mkdir -p "${a[5]}"; cp "$DATA/optitype_result.tsv" "${a[5]}/${a[7]}_result.tsv"; exit 0 ;;   # -o DIR -p PREFIX
  run-t1k) mkdir -p "${a[11]}"; cp "$DATA/t1k_genotype.tsv" "${a[11]}/${a[9]}_genotype.tsv"; exit 0 ;;   # -o PREFIX --od DIR
  docker) d=${a[15]}; s=${a[10]}; mkdir -p "$d/MHC_Class_I" "$d/MHC_Class_II"   # docker run ... pvacseq run VCF SAMPLE HLA A1 A2 A3 OUTDIR
    cp "$DATA/MHC_Class_I/S_T.all_epitopes.aggregated.tsv" "$d/MHC_Class_I/$s.all_epitopes.aggregated.tsv"     # class I result
    cp "$DATA/MHC_Class_II/S_T.all_epitopes.aggregated.tsv" "$d/MHC_Class_II/$s.all_epitopes.aggregated.tsv"; exit 0 ;;   # class II result
  gatk) [ "${1:-}" = MarkDuplicates ] && out "$(dirname "${a[4]}")/$(basename "${a[4]}" .bam).bai" ;;   # MarkDuplicates also writes the .bai index
  bcftools) [ "${1:-}" = +split-vep ] && cat "$DATA/tmb_variants.tsv" ;;   # split-vep prints the variant table that TMB reads
esac
i=0; while [ $i -lt ${#a[@]} ]; do                                  # every other tool: create each file named after an output option
  case ${a[$i]} in -o|-O|--output_file|-1|-2|-M|-j|-h) [ "${a[$i+1]:-/dev/null}" != /dev/null ] && [ "${a[$i+1]}" != - ] && out "${a[$i+1]}" ;; esac
  i=$((i + 1))                                                      # next argument
done
exit 0                                                              # every pretend tool "succeeds"
