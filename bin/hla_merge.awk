# hla_merge.awk - OptiType (class I) + T1K (class II) calls -> one pVACseq allele list, e.g.
#   HLA-A*02:01,HLA-B*07:02,DRB1*15:01,DQA1*01:02-DQB1*06:02
# Usage: awk -F'\t' -v override="HLA-A*01:01;DRB1*03:01" -f hla_merge.awk optitype_result.tsv t1k_genotype.tsv
# Why: pVACseq wants 2-field alleles; class I with an "HLA-" prefix, class II DQ/DP as alpha-beta pairs (the two chains form one molecule).

function two(a,   f) { sub(/^HLA-/, "", a); split(a, f, ":"); return f[1] ":" f[2] }   # HLA-DRB1*15:01:01 -> DRB1*15:01 (pVACseq uses 2 fields)
function add(L, x) { if (x != "" && !(("has", x) in L)) { L[++L["n"]] = x; L["has", x] = 1 } }   # append x to list L once (keeps order)
function join(L, out,   i) { i = 0; while (++i <= L["n"]) out = out (out == "" ? "" : ",") L[i]; return out }   # list -> "a,b,c"
function pairs(A, B, OUT,   i, j) { i = 0; while (++i <= A["n"]) { j = 0; while (++j <= B["n"]) add(OUT, A[i] "-" B[j]) } }   # every alpha x beta

$2 == "A1" { opt = 1; next }                                    # OptiType header row: the next row holds the 6 class I alleles
opt == 1 { i = 1; while (++i <= 7) if ($i != "") add(OPT, "HLA-" two($i)); opt = 0; next }   # columns 2-7 = A1 A2 B1 B2 C1 C2
$1 ~ /^HLA-/ && NF >= 8 {                                       # a T1K genotype row: gene, n, allele1, abundance, quality, allele2, abundance, quality
  g = $1; sub(/^HLA-/, "", g)                                   # gene name without prefix (A, DRB1, DQA1 ...)
  if ($3 != "." && $5 >= q) call[g, ++cnt[g]] = two($3)          # first allele, if T1K is confident enough (q = min quality, default 0)
  if ($6 != "." && $8 >= q) call[g, ++cnt[g]] = two($6)          # second allele (same as the first when homozygous)
  i = 0; while (++i <= cnt[g]) {                                # sort each allele into the list it belongs to
    a = call[g, i]                                              # current allele
    if (g ~ /^[ABC]$/) add(T1KI, "HLA-" a)                      # class I from T1K: only used if OptiType gave nothing
    else if (g ~ /^DRB[1345]$/) add(DR, a)                      # DR beta chains pair with the invariant DRA, so they stand alone
    else if (g == "DQA1") add(DQA1, a); else if (g == "DQB1") add(DQB1, a)   # DQ alpha / beta chains: paired up at the end
    else if (g == "DPA1") add(DPA1, a); else if (g == "DPB1") add(DPB1, a)   # DP alpha / beta chains: paired up at the end
  }
}
END {
  pairs(DQA1, DQB1, DR); pairs(DPA1, DPB1, DR)                  # DQ = DQA1 x DQB1, DP = DPA1 x DPB1 (all combinations)
  cI = OPT["n"] ? join(OPT) : join(T1KI); cII = join(DR)        # class I: OptiType first (it is the more accurate class I caller)
  n = split(override, o, /[;,]/); i = 0                         # clinical typing from the samplesheet / --hla, if any
  while (++i <= n) { x = o[i]; gsub(/ /, "", x); if (x == "") continue   # one override allele at a time, spaces removed
    if (x ~ /^(HLA-)?[ABC]\*/) add(UI, (x ~ /^HLA-/ ? x : "HLA-" x)); else add(UII, x) }   # class I gets the HLA- prefix, class II as given
  if (UI["n"]) cI = join(UI); if (UII["n"]) cII = join(UII)     # an override replaces the computed alleles of the same class only
  print cI (cI != "" && cII != "" ? "," : "") cII               # final comma-separated list for pvacseq run
}
