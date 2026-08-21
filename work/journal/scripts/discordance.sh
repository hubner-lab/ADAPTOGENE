#!/bin/bash
set -euo pipefail
VCF=/mnt/data/eugene/ADAPTOGENE/Trifolium_results/_work/maf0.05_miss0.3_smiss0.7/Trifolium_chr1-6.vcf

check_pair() {
  local s1=$1 s2=$2 label=$3
  bcftools query -s "${s1},${s2}" -f '[%GT\t]\n' "$VCF" 2>/dev/null | \
  awk -v s1="$s1" -v s2="$s2" -v label="$label" '
    {
      gt1=$1; gt2=$2
      gsub(/\|/,"/",gt1); gsub(/\|/,"/",gt2)
      if (gt1=="./." || gt2=="./.") { miss++; next }
      if (gt1==gt2) { same++; next }
      # opposite homozygote: 0/0 vs 1/1
      if ((gt1=="0/0" && gt2=="1/1") || (gt1=="1/1" && gt2=="0/0")) { opp++ }
      diff++
    }
    END {
      tot = same+diff
      printf "%-12s %-10s %-10s n_compared=%-8d same=%-8d opp_hom=%-6d discordance=%.5f opp_rate=%.5f\n", \
        label, s1, s2, tot, same, opp, (tot>0?diff/tot:0), (tot>0?opp/tot:0)
    }'
}

echo "=== TOP (>=0.999) ==="
check_pair 28955-1 28955-4 "top1_0.9997"
check_pair 29024-10 29024-11 "top2_0.9993"
check_pair 28955-1 28955-24 "top3_0.9992"

echo "=== P90 (~0.73) ==="
check_pair 28955-3 28985-21 "p90_0.746"
check_pair 28967-2 28968-24 "p90_0.726"

echo "=== BASELINE (~0.55, lowest observed) ==="
check_pair 28960-15 29072-2 "low_0.549"
check_pair 28984-1 29017-7 "low_0.547"

echo "=== BOUNDARY (0.99-0.999) ==="
check_pair 28985-14 28985-5   "boundary_0.995"
check_pair 29072-2  29072-3   "boundary_0.996"
check_pair 28960-15 28960-17  "boundary_0.998"
check_pair 28955-1  28955-3   "boundary_0.999"
