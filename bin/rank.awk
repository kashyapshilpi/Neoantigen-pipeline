# rank.awk - join pVACseq's class I and class II aggregated reports into one row per variant, with sort keys in front.
# Usage: awk -F'\t' -f rank.awk MHC_Class_I/<S>.all_epitopes.aggregated.tsv MHC_Class_II/<S>.all_epitopes.aggregated.tsv
# Output columns: 4 sort keys, then variant, gene, aa_change, dna_vaf, best_tier, both_classes, then tier/IC50/allele/peptide for I and II.
# Why: a variant that gives good class I (CD8 T cells) AND class II (CD4 T cells) epitopes is the most promising vaccine target.

FNR == 1 { k = 0; while (++k <= NF) col[$k] = k; c = (FILENAME ~ /Class_II/) ? 2 : 1; next }   # header: column name -> number; which class
{
  id = $col["ID"]; if (!(id in seen)) { seen[id] = 1; ids[++n] = id }   # remember each variant once, in first-seen order
  gene[id] = $col["Gene"]; aa[id] = $col["AA Change"]; vaf[id] = $col["DNA VAF"]   # variant details (same in both classes)
  tier[id, c] = $col["Tier"]; ic[id, c] = $col["IC50 MT"]; al[id, c] = $col["Allele"]; pep[id, c] = $col["Best Peptide"]   # per-class best epitope
}
function tr(t) { return t == "Pass" ? 0 : t == "Anchor" ? 1 : t == "Subclonal" ? 2 : t == "LowExpr" ? 3 : t == "NoExpr" ? 4 : t == "Poor" ? 5 : 6 }   # tier -> rank
END {
  i = 0; while (++i <= n) { id = ids[i]                         # every variant
    t1 = tr(tier[id, 1]); t2 = tr(tier[id, 2])                  # tier rank in each class (6 = no epitope in that class)
    best = t1 < t2 ? t1 : t2; both = (t1 <= 1 && t2 <= 1)       # best tier overall; both = good (Pass/Anchor) in class I AND class II
    b50 = ic[id, 1] != "" ? ic[id, 1] : (ic[id, 2] != "" ? ic[id, 2] : 1e9)   # strongest binding: class I IC50 first (better calibrated), else class II
    printf "%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", best, both, b50, (vaf[id] == "" ? 0 : vaf[id]), id, gene[id], aa[id], vaf[id], \
      (best == 6 ? "-" : (t1 <= t2 ? tier[id, 1] : tier[id, 2])), (both ? "yes" : "no")   # sort keys + variant details
    c = 0; while (++c <= 2) printf "\t%s\t%s\t%s\t%s", tier[id, c], ic[id, c], al[id, c], pep[id, c]   # class I, then class II columns
    print ""                                                    # end of the row
  }
}
