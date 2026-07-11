#!/bin/bash
set -euo pipefail

# Download and process Arabidopsis 1001 Genomes full VCF
# Source: 1001genomes.org v3.1 release
# Output: data/Arabidopsis_1001G.vcf with S-prefixed sample names

DATA_DIR="data"
OUTVCF="${DATA_DIR}/Arabidopsis_1001G.vcf"
URL="https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
RAW_GZ="${DATA_DIR}/1001genomes_raw.vcf.gz"

echo "=== Step 1: Download full VCF (~18 GB) ==="
if [ ! -f "${RAW_GZ}" ]; then
    wget -O "${RAW_GZ}" "${URL}"
else
    echo "  Already downloaded: ${RAW_GZ}"
fi

echo "=== Step 2: Index VCF ==="
bcftools index "${RAW_GZ}"

echo "=== Step 3: Extract biallelic SNPs only ==="
# Filter to biallelic SNPs, remove indels, keep only chr 1-5
# Also add S prefix to sample names
echo "  Filtering with bcftools..."
bcftools view -m2 -M2 -v snps "${RAW_GZ}" \
    --regions 1,2,3,4,5 \
    -O v \
    | awk 'BEGIN{OFS="\t"} /^#CHROM/{for(i=10;i<=NF;i++) $i="S"$i} {print}' \
    > "${OUTVCF}"

# Count results
N_SAMPLES=$(grep "^#CHROM" "${OUTVCF}" | tr '\t' '\n' | tail -n +10 | wc -l)
N_SNPS=$(grep -vc "^#" "${OUTVCF}")
echo "  Result: ${N_SNPS} biallelic SNPs, ${N_SAMPLES} samples"
echo "  Output: ${OUTVCF}"

echo "=== Done ==="
