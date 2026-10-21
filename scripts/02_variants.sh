#!/usr/bin/env bash
# 02_variants.sh - somatic variants (tumor vs normal, Mutect2) and germline variants (normal, HaplotypeCaller).
# Usage: bash 02_variants.sh <PATIENT>
set -euo pipefail                                                   # stop at the first error
source "$(dirname "$0")/../config.sh"                               # load paths and settings
P=$1; T=$OUT/${P}_T/bam/${P}_T.final.bam; N=$OUT/${P}_N/bam/${P}_N.final.bam; D=$OUT/$P/variants; mkdir -p "$D"   # inputs + output folder

M2=(); [ -n "$GNOMAD" ] && M2+=(--germline-resource "$GNOMAD"); [ -n "$PON" ] && M2+=(--panel-of-normals "$PON")   # optional Mutect2 resources
# somatic calling: a variant seen in the tumor but not in the matched normal is acquired by the tumor (only these make neoantigens)
gatk Mutect2 -R "$REF" -I "$T" -I "$N" -normal "${P}_N" -L "$TARGETS" "${M2[@]}" --f1r2-tar-gz "$D/f1r2.tar.gz" -O "$D/raw.vcf.gz"
gatk LearnReadOrientationModel -I "$D/f1r2.tar.gz" -O "$D/orientation.tar.gz"   # models FFPE / oxidation damage that fakes C>T and G>T calls
F=(--ob-priors "$D/orientation.tar.gz")                             # filters to apply: start with the orientation model
if [ -n "$COMMON_SNPS" ]; then                                      # contamination check needs a common-SNP VCF
  gatk GetPileupSummaries -I "$T" -V "$COMMON_SNPS" -L "$COMMON_SNPS" -O "$D/tumor.pileups.table"    # allele counts at common SNPs (tumor)
  gatk GetPileupSummaries -I "$N" -V "$COMMON_SNPS" -L "$COMMON_SNPS" -O "$D/normal.pileups.table"   # same in the normal
  gatk CalculateContamination -I "$D/tumor.pileups.table" -matched "$D/normal.pileups.table" -O "$D/contamination.table" --tumor-segmentation "$D/segments.table"   # % of other-person DNA
  F+=(--contamination-table "$D/contamination.table" --tumor-segmentation "$D/segments.table")      # and use it when filtering
fi
gatk FilterMutectCalls -R "$REF" -V "$D/raw.vcf.gz" "${F[@]}" -O "$D/filtered.vcf.gz"   # mark artefacts, germline and contamination calls
# split multi-allelic sites, left-align indels and keep PASS only: pVACseq needs one clean ALT allele per line
bcftools norm -m -any -f "$REF" "$D/filtered.vcf.gz" | bcftools view -f PASS -Oz -o "$D/$P.somatic.vcf.gz" && tabix -f "$D/$P.somatic.vcf.gz"
# germline calling on the normal: inherited variants near a tumor mutation change the peptide too, so pVACseq needs them for phasing
gatk HaplotypeCaller -R "$REF" -I "$N" -L "$TARGETS" -O "$D/$P.germline.vcf.gz"
bcftools stats "$D/$P.somatic.vcf.gz" > "$D/$P.somatic.stats.txt"   # QC: number of SNVs / indels, Ts/Tv ratio
