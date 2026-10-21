#!/usr/bin/env bash
# rank.sh - ranked neoantigen table from pVACseq class I + II aggregated reports.
# Usage: bash rank.sh <class_I.aggregated.tsv> <class_II.aggregated.tsv>   (a missing class can be given as /dev/null)
# Order: best tier first -> good in BOTH classes -> strongest binding (lowest IC50) -> highest tumor VAF (more tumor cells carry it).
set -euo pipefail                                                   # stop at the first error
printf 'rank\tvariant\tgene\taa_change\tdna_vaf\tbest_tier\tboth_classes\tI_tier\tI_ic50\tI_allele\tI_peptide\tII_tier\tII_ic50\tII_allele\tII_peptide\n'   # header
awk -F'\t' -f "$(dirname "$0")/rank.awk" "$1" "$2" \
  | sort -t$'\t' -k1,1n -k2,2nr -k3,3g -k4,4gr -k5,5 \
  | cut -f5- | awk '{ print NR "\t" $0 }'                           # join -> sort by the 4 keys (+ variant id as a tie-breaker) -> drop keys -> number
