#!/usr/bin/env bash
# run_tests.sh - checks the awk helpers on small made-up inputs in test/data (no real data or bioinformatics tools needed).
# Usage: bash test/run_tests.sh        (exit code 0 = all tests passed; CI runs this on every push)
set -uo pipefail                                                    # unset variables / failed pipes are errors; keep going after a failed test
cd "$(dirname "$0")" || exit 1; B=../bin; T=$'\t'; fail=0                     # run from the test folder; B = helper folder; T = tab for -F
TMB="-v mb=2 -v depth=20 -v alt=4 -v vaf=0.05 -v nvaf=0.02 -v high=10"   # TMB settings used by the TMB tests (target = 2 Mb)
while IFS='|' read -r name cmd; do                                  # one test per line below: expected-output file | command
  if diff <(eval "$cmd") "expected/$name" > /dev/null; then echo "PASS $name"   # same output as expected -> pass
  else echo "FAIL $name"; diff <(eval "$cmd") "expected/$name"; fail=1; fi   # otherwise show the difference and remember the failure
done <<'TESTS'
hla_optitype_t1k.txt|awk -F"$T" -f $B/hla_merge.awk data/optitype_result.tsv data/t1k_genotype.tsv
hla_t1k_only.txt|awk -F"$T" -f $B/hla_merge.awk /dev/null data/t1k_genotype.tsv
hla_override.txt|awk -F"$T" -v override="A*01:01; HLA-B*08:01" -f $B/hla_merge.awk data/optitype_result.tsv data/t1k_genotype.tsv
hla_min_quality.txt|awk -F"$T" -v q=39 -f $B/hla_merge.awk /dev/null data/t1k_genotype.tsv
tmb_summary.tsv|awk -F"$T" $TMB -v list=/dev/null -f $B/tmb.awk data/tmb_variants.tsv
tmb_counted.tsv|awk -F"$T" $TMB -v list=/dev/stderr -f $B/tmb.awk data/tmb_variants.tsv 2>&1 > /dev/null
rank_both_classes.tsv|bash $B/rank.sh data/MHC_Class_I/S_T.all_epitopes.aggregated.tsv data/MHC_Class_II/S_T.all_epitopes.aggregated.tsv
rank_class_I_only.tsv|bash $B/rank.sh data/MHC_Class_I/S_T.all_epitopes.aggregated.tsv /dev/null
TESTS
exit $fail                                                          # non-zero exit makes CI fail
