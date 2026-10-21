#!/usr/bin/env bash
# 07_rank_report.sh - one ranked candidate table (class I + II together) and a simple HTML report with TMB, HLA and top candidates.
# Usage: bash 07_rank_report.sh <PATIENT>
set -euo pipefail                                                   # stop at the first error
source "$(dirname "$0")/../config.sh"                               # load paths and settings
P=$1; V=$OUT/$P/pvacseq; D=$OUT/$P/report; mkdir -p "$D"           # patient, pVACseq output, report folder
f() { [ -s "$V/MHC_Class_$1/${P}_T.all_epitopes.aggregated.tsv" ] && echo "$V/MHC_Class_$1/${P}_T.all_epitopes.aggregated.tsv" || echo /dev/null; }   # a class may be empty
bash "$(dirname "$0")/../bin/rank.sh" "$(f I)" "$(f II)" > "$D/$P.neoantigens_ranked.tsv"   # the main result: every candidate, best first
# HTML report: TMB line + HLA line + the top candidates as a table (awk turns tab-separated rows into <tr><td> cells)
{ echo "<html><head><meta charset='utf-8'><title>$P neoantigens</title><style>body{font-family:sans-serif}td,th{border:1px solid #ccc;padding:3px 6px}table{border-collapse:collapse;font-size:13px}</style></head><body>"
  echo "<h1>Neoantigen report: $P</h1>"                            # title
  awk -F'\t' 'NR == 2 { printf "<p><b>TMB:</b> %s mutations/Mb (%s; %s nonsynonymous mutations in %s Mb)</p>\n", $4, $5, $1, $3 }' "$OUT/$P/tmb/$P.tmb.tsv"   # TMB summary
  echo "<p><b>HLA alleles:</b> $(sed 's/,/, /g' "$OUT/$P/hla/$P.hla.txt")</p><h2>Top $TOP_N candidates</h2><table>"   # HLA used for prediction
  head -n $((TOP_N + 1)) "$D/$P.neoantigens_ranked.tsv" | awk -F'\t' '{ t = NR == 1 ? "th" : "td"; printf "<tr>"; i = 0; while (++i <= NF) printf "<%s>%s</%s>", t, $i, t; print "</tr>" }'   # header + top N rows
  echo "</table><p>Ranking: best pVACseq tier, then good in both MHC classes, then lowest IC50, then highest tumor VAF. Research use only.</p></body></html>"
} > "$D/$P.neoantigen_report.html"                                  # write it all to one HTML file
echo "$(($(wc -l < "$D/$P.neoantigens_ranked.tsv") - 1)) candidates ranked"   # how many candidates (minus the header line)
