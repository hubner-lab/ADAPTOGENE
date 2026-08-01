#!/usr/bin/env bash
# =============================================================================
# score_ladders.sh — score every hyperparameter ladder rung against the truth.
#
# mode=pregea selects rungs on CALIBRATION grounds (p-value histogram flatness
# for LFMM K, lambda_GC for EMMAX n_pcs, adjusted R2 for RDA condition_pcs). It
# never looks at the causal loci. This driver adds that missing view: each rung's
# persisted p-value table is scored against causal_loci.tsv, so rungs can be
# ranked by AUC-PR and best-F1 instead of by calibration proxies alone.
#
# Rungs consumed:
#   LFMM   {PROJECT}_results/_intermediate/pregea/lfmm/K{k}/pvalues.tsv
#   EMMAX  {PROJECT}_results/_intermediate/pregea/emmax/npc{n}/pvalues.tsv
#   RDA    benchmarks/laruson_sweep/{PROJECT}/rda/pc{n}/RDA_pvalues.tsv
#          (pregea's RDA block persists only diagnostics, so sweep_rda.sh
#           produces these separately)
#
# LFMM rungs are scored TWICE — raw, and with genomic control applied — because
# PreGEA runs LFMM with genomic_control = FALSE while production mode=gea uses
# TRUE, and at the full cohort that correction had lambda = 29.9.
#
# Usage: benchmarks/score_ladders.sh [PROJECT] [TRUTH_TSV]
# =============================================================================
set -euo pipefail

PROJECT="${1:-LARUSON1K}"
ROOT="${PIPELINE_ROOT:-/pipeline}"
TRUTH="${2:-${ROOT}/data/laruson_1k/causal_loci.tsv}"
RES="${ROOT}/${PROJECT}_results"
OUT="${ROOT}/benchmarks/laruson_eval/${PROJECT}_ladders.tsv"
EVAL="${ROOT}/benchmarks/eval_detection.R"

[[ -f "$TRUTH" ]] || { echo "ERROR: truth table not found: $TRUTH" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"
rm -f "$OUT"

score() {   # label, pvalues  — raw p-values, no genomic control
    local label="$1" pvals="$2"
    [[ -s "$pvals" ]] || { echo "SKIP $label — missing or empty: $pvals" >&2; return 0; }
    Rscript "$EVAL" --mode=rank --pvalues="$pvals" --truth="$TRUTH" \
        --label="$label" --out="$OUT" 2>&1 | grep -E "AUC-PR|TP@top" || true
    local adj val
    for spec in "qval 0.1" "bonf 0.05" "top 20"; do
        adj="${spec%% *}"; val="${spec##* }"
        Rscript "$EVAL" --mode=threshold --pvalues="$pvals" --truth="$TRUTH" \
            --label="$label" --adjust="$adj" --value="$val" --out="$OUT" \
            2>&1 | grep -E "^${label}\|" || true
    done
}

score_gc() {   # label, pvalues — same rung with genomic control applied
    local label="$1" pvals="$2"
    [[ -s "$pvals" ]] || { echo "SKIP $label (gc) — missing or empty: $pvals" >&2; return 0; }
    Rscript "$EVAL" --mode=rank --pvalues="$pvals" --truth="$TRUTH" \
        --label="$label" --genomic-control=yes --out="$OUT" 2>&1 | grep -E "AUC-PR" || true
    local adj val
    for spec in "qval 0.1" "bonf 0.05" "top 20"; do
        adj="${spec%% *}"; val="${spec##* }"
        Rscript "$EVAL" --mode=threshold --pvalues="$pvals" --truth="$TRUTH" \
            --label="$label" --adjust="$adj" --value="$val" --genomic-control=yes \
            --out="$OUT" 2>&1 | grep -E "^${label}\+gc\|" || true
    done
}

echo "=================== LFMM K ladder (raw p-values) ==================="
for d in "${RES}"/_intermediate/pregea/lfmm/K*/; do
    [[ -d "$d" ]] || continue
    score "LFMM_$(basename "$d")" "${d}pvalues.tsv"
done

echo "=================== LFMM K ladder (genomic control applied) ==================="
for d in "${RES}"/_intermediate/pregea/lfmm/K*/; do
    [[ -d "$d" ]] || continue
    score_gc "LFMM_$(basename "$d")" "${d}pvalues.tsv"
done

echo "=================== EMMAX n_pcs ladder ==================="
for d in "${RES}"/_intermediate/pregea/emmax/npc*/; do
    [[ -d "$d" ]] || continue
    score "EMMAX_$(basename "$d")" "${d}pvalues.tsv"
done

echo "=================== RDA condition_pcs ladder ==================="
for d in "${ROOT}"/benchmarks/laruson_sweep/"${PROJECT}"/rda/pc*/; do
    [[ -d "$d" ]] || continue
    score "RDA_$(basename "$d")" "${d}RDA_pvalues.tsv"
done

echo
echo "INFO: ladder scores written to $OUT"
