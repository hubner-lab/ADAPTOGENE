#!/bin/bash
set -euo pipefail

# ADAPTOGENE Pipeline — WBDC Dataset
# Runs all pipeline modes sequentially
# haplotype mode included since epsilon_selected is already set (0.8)

CONFIG="config_WBDC.yaml"
CORES=70
MEMORY="300g"
IMAGE="adaptogene:latest"

MODES=(
    processing
    structure
    structure_K
    association
    association_phenotypes
    overlapping
    regionplot
    maladaptation
    haplotype_scan
    haplotype
)

for MODE in "${MODES[@]}"; do
    echo "============================================="
    echo "  Starting mode: ${MODE}"
    echo "  $(date)"
    echo "============================================="

    docker run \
        --user "$(id -u):$(id -g)" \
        --rm \
        --memory="${MEMORY}" \
        --cpus="${CORES}" \
        -v "$PWD:/pipeline" \
        "${IMAGE}" \
        snakemake -c"${CORES}" -s Snakefile \
            --config mode="${MODE}" \
            --configfile "${CONFIG}" \
            --scheduler greedy \
            --rerun-incomplete

    echo "  Finished mode: ${MODE} at $(date)"
    echo ""
done

echo "============================================="
echo "  All modes complete!"
echo "============================================="
