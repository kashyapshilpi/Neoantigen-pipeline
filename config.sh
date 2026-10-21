#!/usr/bin/env bash
# config.sh - every path and setting in ONE place. All scripts load it with 'source', so you only edit this file.
# ${VAR:-default} means: use VAR if it is already set in your shell, otherwise use the default.

# ---- reference files (one-time downloads, same for every patient) ----
REF=${REF:-$HOME/ref/GRCh38.fa}                          # genome FASTA; needs .fai, .dict and a bwa-mem2 index next to it
TARGETS=${TARGETS:-$HOME/ref/exome_targets.bed}          # exome capture BED: limits calling to the exome and is the TMB denominator
KNOWN_SITES=${KNOWN_SITES:-}                             # space-separated dbSNP / Mills VCFs for BQSR; empty = skip BQSR
GNOMAD=${GNOMAD:-}                                       # af-only-gnomad VCF: tells Mutect2 which variants are common germline
PON=${PON:-}                                             # panel of normals VCF: removes recurring sequencing artefacts
COMMON_SNPS=${COMMON_SNPS:-}                             # small_exac_common VCF: used to estimate tumor contamination
VEP_CACHE=${VEP_CACHE:-$HOME/ref/vep_cache}              # offline Ensembl VEP cache (no internet needed while annotating)
VEP_PLUGINS=${VEP_PLUGINS:-$HOME/ref/vep_plugins}        # folder with Frameshift.pm and Wildtype.pm (pVACseq needs both)
T1K_REF=${T1K_REF:-$HOME/ref/t1k/hla_dna_seq.fa}         # T1K HLA allele sequences (made once with t1k-build.pl)

# ---- tool settings ----
HLA_REGION=chr6:28510120-33480577                        # GRCh38 MHC region: HLA reads come from here (plus unmapped reads)
PVAC_IMAGE=griffithlab/pvactools:7.1.5                   # official pVACtools Docker image; it already contains the IEDB tools
ALGOS="NetMHCpan MHCflurry NetMHCIIpan"                   # binding predictors: NetMHCpan + MHCflurry for class I, NetMHCIIpan for class II
LEN_I=8,9,10,11                                          # class I peptides are 8-11 amino acids long
LEN_II=12,13,14,15,16,17,18                              # class II peptides are longer (12-18)
IC50=500                                                 # nM; below 500 is the usual "binder" cutoff
TMB_DEPTH=20; TMB_ALT=4; TMB_VAF=0.05; TMB_NVAF=0.02     # a variant counts for TMB only with enough reads, alt reads, VAF and a clean normal
TMB_HIGH=10                                              # mutations per Mb; >= 10 is the common "TMB-high" cutoff
TOP_N=25                                                 # how many top candidates the HTML report shows

# ---- compute and output ----
THREADS=${THREADS:-16}                                   # CPU threads per tool
OUT=${OUT:-$PWD/results}                                 # every result goes under this folder
