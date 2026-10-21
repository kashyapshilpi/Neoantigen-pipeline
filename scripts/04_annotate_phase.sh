#!/usr/bin/env bash
# 04_annotate_phase.sh - build the two pVACseq inputs:
#   1. somatic VCF annotated by VEP with the Frameshift + Wildtype plugins (so pVACseq can build mutant AND normal peptides)
#   2. phased VCF of somatic + germline variants (a nearby variant on the same DNA copy changes the peptide sequence)
# Usage: bash 04_annotate_phase.sh <PATIENT>
set -euo pipefail                                                   # stop at the first error
source "$(dirname "$0")/../config.sh"                               # load paths and settings
P=$1; V=$OUT/$P/variants; D=$OUT/$P/pvacseq_input; mkdir -p "$D"   # patient, variant folder, output folder
T=$OUT/${P}_T/bam/${P}_T.final.bam                                  # tumor BAM: its reads show which variants sit on the same copy

vep_run() {                                                         # vep_run <in.vcf> <out.vcf.gz>: the VEP options pVACseq requires
  vep --input_file "$1" --output_file "$2" --format vcf --vcf --compress_output bgzip --symbol --terms SO --tsl --biotype --hgvs \
    --fasta "$REF" --offline --cache --dir_cache "$VEP_CACHE" --plugin Frameshift --plugin Wildtype --dir_plugins "$VEP_PLUGINS" \
    --fork "$THREADS" --force_overwrite && tabix -f -p vcf "$2"     # --symbol/--tsl/--biotype/--hgvs fill pVACseq's columns; index for fast reading
}
vep_run "$V/$P.somatic.vcf.gz" "$D/$P.somatic.vep.vcf.gz"           # input 1: annotated somatic variants

bcftools view -s "${P}_T" "$V/$P.somatic.vcf.gz" -Oz -o "$D/som.T.vcf.gz" && tabix -f "$D/som.T.vcf.gz"   # somatic: keep the tumor column only
echo "${P}_N ${P}_T" > "$D/rename.txt"                              # germline calls are named after the normal; rename them to the tumor
bcftools reheader -s "$D/rename.txt" "$V/$P.germline.vcf.gz" -o "$D/germ.T.vcf.gz" && tabix -f "$D/germ.T.vcf.gz"   # so both files share one sample
bcftools concat -a -D "$D/som.T.vcf.gz" "$D/germ.T.vcf.gz" | bcftools sort -Oz -o "$D/combined.vcf.gz" && tabix -f "$D/combined.vcf.gz"   # one sorted VCF
# read-based phasing with the tumor reads (WhatsHap replaces GATK3 ReadBackedPhasing from the pVACseq manual, which is discontinued)
whatshap phase --reference "$REF" -o "$D/phased.vcf" "$D/combined.vcf.gz" "$T"
vep_run "$D/phased.vcf" "$D/$P.phased.vep.vcf.gz"                   # input 2: the phased proximal-variants file, also VEP-annotated
rm -f "$D"/som.T.vcf.gz* "$D"/germ.T.vcf.gz* "$D/combined.vcf.gz"* "$D/phased.vcf" "$D/rename.txt"   # clean up temporary files
