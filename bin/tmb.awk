# tmb.awk - tumor mutational burden (TMB) = qualifying nonsynonymous somatic mutations per megabase of exome.
# Input (tab-separated, one row per variant): CHROM POS REF ALT Consequence SYMBOL tumor_AD normal_AD   (AD = "ref_reads,alt_reads")
# Usage: awk -F'\t' -v mb=35.1 -v depth=20 -v alt=4 -v vaf=0.05 -v nvaf=0.02 -v high=10 -v list=counted.tsv -f tmb.awk variants.tsv
# Why these filters: low-depth / low-VAF calls are often noise, and a variant also seen in the normal is probably germline.

BEGIN { NS = "missense|frameshift|stop_gained|stop_lost|start_lost|inframe_insertion|inframe_deletion|protein_altering"   # changes the protein
        print "chrom\tpos\tref\talt\tgene\tconsequence\tdepth\tvaf" > list }   # header of the list of counted variants
{
  split($7, t, ","); split($8, n, ",")                          # tumor and normal read counts: [1] = ref reads, [2] = alt reads
  d = t[1] + t[2]; v = d ? t[2] / d : 0                          # tumor depth and variant allele fraction (VAF)
  nd = n[1] + n[2]; nv = nd ? n[2] / nd : 0                      # VAF in the normal (should be ~0 for a real somatic mutation)
  if (d < depth || t[2] < alt || v < vaf || nv > nvaf) next      # skip weak or germline-looking calls
  coding++                                                      # passes the quality filters (all coding/non-coding classes)
  if ($5 !~ NS) next                                            # TMB counts only mutations that change the protein
  ns++; printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.3f\n", $1, $2, $3, $4, $6, $5, d, v > list   # count it and list it
}
END {
  if (mb <= 0) { print "ERROR: target size (mb) must be > 0" > "/dev/stderr"; exit 1 }   # avoid dividing by zero
  tmb = ns / mb                                                 # mutations per Mb
  printf "nonsynonymous\tqualifying_all\ttarget_mb\ttmb_mut_per_mb\tclass\n%d\t%d\t%.3f\t%.2f\t%s\n", ns, coding, mb, tmb, (tmb >= high ? "TMB-high" : "TMB-low")
}
