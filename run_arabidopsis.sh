#!/bin/bash
set -euo pipefail

# ADAPTOGENE Pipeline — Arabidopsis 1001 Genomes
# Runs all pipeline modes sequentially (except haplotype, which needs epsilon_selected)

CONFIG="config_arabidopsis.yaml"
CORES=100
MEMORY="400g"
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
)

PROJ="Arabidopsis"

# Skip LD decay by touching outputs before structure_K
skip_ld_decay() {
    local RESDIR="${PROJ}_results"
    mkdir -p "${RESDIR}/_intermediate/ld_decay/"
    mkdir -p "${RESDIR}/structure_k/tables/"
    mkdir -p "${RESDIR}/structure_k/plots/ld_decay/"
    touch "${RESDIR}/_intermediate/ld_decay/prep_done.flag"
    touch "${RESDIR}/_intermediate/ld_decay/gw_done.flag"
    echo -e "group\tscope\thalf_decay_bp\tr2_02_bp\tn_samples\tn_pairs\tr2_intercept\tmethod" \
        > "${RESDIR}/structure_k/tables/ld_decay_half_distances.tsv"
    touch "${RESDIR}/structure_k/plots/ld_decay/ld_decay_genome_wide.png"
    touch "${RESDIR}/structure_k/plots/ld_decay/ld_decay_genome_wide.svg"
    echo "  [SKIP] LD decay outputs touched"
}

for MODE in "${MODES[@]}"; do
    echo "============================================="
    echo "  Starting mode: ${MODE}"
    echo "  $(date)"
    echo "============================================="

    # Touch LD decay outputs before structure_K
    if [ "${MODE}" = "structure_K" ]; then
        skip_ld_decay
    fi

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
echo "  NOTE: To run haplotype mode, review clustree"
echo "  plots and set haplotype.epsilon_selected in"
echo "  ${CONFIG}, then run mode=haplotype"
echo "============================================="
