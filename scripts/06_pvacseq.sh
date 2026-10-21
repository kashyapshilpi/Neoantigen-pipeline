#!/usr/bin/env bash
# 06_pvacseq.sh - predict which mutant peptides bind the patient's MHC class I and class II molecules (pVACseq).
# Usage: bash 06_pvacseq.sh <PATIENT>
set -euo pipefail                                                   # stop at the first error
source "$(dirname "$0")/../config.sh"                               # load paths and settings
P=$1; I=$OUT/$P/pvacseq_input; D=$OUT/$P/pvacseq; mkdir -p "$D"    # patient, inputs from step 4, output folder
HLA=$(cat "$OUT/$P/hla/$P.hla.txt")                                 # alleles from step 3, e.g. HLA-A*02:01,DRB1*15:01
read -ra AL <<< "$ALGOS"                                            # algorithms as separate words (pVACseq takes them as separate arguments)
# run inside the official image (IEDB tools are pre-installed there); mount the results folder at the same path; run as you, not root
docker run --rm -u "$(id -u):$(id -g)" -v "$OUT:$OUT" "$PVAC_IMAGE" pvacseq run \
  "$I/$P.somatic.vep.vcf.gz" "${P}_T" "$HLA" "${AL[@]}" "$D" \
  -e1 "$LEN_I" -e2 "$LEN_II" -b "$IC50" -p "$I/$P.phased.vep.vcf.gz" --normal-sample-name "${P}_N" \
  --iedb-install-directory /opt/iedb --tdna-vaf 0.25 --normal-vaf 0.02 --pass-only -t "$THREADS"
# what the options mean: -e1/-e2 peptide lengths for class I/II, -b IC50 binder cutoff (nM), -p phased nearby variants,
# --normal-sample-name lets pVACseq check normal VAF, --tdna-vaf 0.25 = clonal cutoff used for the 'Subclonal' tier, --pass-only = PASS calls only
ls "$D"                                                             # show what was produced (MHC_Class_I, MHC_Class_II, combined)
