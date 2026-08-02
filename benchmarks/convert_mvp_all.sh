#!/usr/bin/env bash
# convert_mvp_all.sh -- run convert_mvp.R over every seed in the manifest whose three raw
# members are present, then roll the per-seed provenance rows up into one table.
#
# Seeds are independent, so conversions run in parallel. Re-running is cheap and idempotent:
# a seed whose outputs already exist is converted again from the same inputs.
#
# Usage:  benchmarks/convert_mvp_all.sh [JOBS]
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JOBS="${1:-8}"
SEEDS_TSV="$ROOT/benchmarks/mvp_seeds.tsv"
RAW="$ROOT/data/mvp/raw"
DOCKER=(nix shell nixpkgs#docker-client -c docker run --user "$(id -u):$(id -g)" --rm
        -e USER=pipeline -v "$ROOT:/pipeline" adaptogene:latest)

mapfile -t SEEDS < <(awk 'NR>1 {print $1}' "$SEEDS_TSV")
ready=(); skipped=()
for s in "${SEEDS[@]}"; do
  if [[ -s "$RAW/genotypes/${s}_Rout_Gmat_sample.txt.gz" &&
        -s "$RAW/mutations/${s}_Rout_muts_full.txt.gz" &&
        -s "$RAW/individuals/${s}_Rout_ind_subset.txt.gz" ]]; then
    ready+=("$s")
  else
    skipped+=("$s")
  fi
done
echo "ready: ${#ready[@]}/${#SEEDS[@]} -- ${ready[*]}"
# Never let an incomplete fetch look like a complete conversion.
[[ ${#skipped[@]} -gt 0 ]] && echo "SKIPPED (raw members missing): ${skipped[*]}" >&2

mkdir -p "$ROOT/data/mvp/logs"
run_one() {
  local s="$1"
  "${DOCKER[@]}" Rscript /pipeline/benchmarks/convert_mvp.R \
      --seed="$s" --raw-dir=/pipeline/data/mvp/raw \
      --out-dir="/pipeline/data/mvp/MVP${s}" > "$ROOT/data/mvp/logs/convert_${s}.log" 2>&1
  local rc=$?
  if [[ $rc -eq 0 ]]; then echo "  OK   $s"; else echo "  FAIL $s (rc=$rc, see data/mvp/logs/convert_${s}.log)" >&2; fi
  return $rc
}
export -f run_one; export ROOT
status=0; running=0
for s in "${ready[@]}"; do
  run_one "$s" & running=$((running+1))
  if [[ $running -ge $JOBS ]]; then wait -n || status=1; running=$((running-1)); fi
done
while [[ $running -gt 0 ]]; do wait -n || status=1; running=$((running-1)); done

echo "=== rolling up provenance ==="
"${DOCKER[@]}" Rscript -e '
suppressPackageStartupMessages(library(data.table))
f <- Sys.glob("/pipeline/data/mvp/MVP*/conversion_provenance.tsv")
if (!length(f)) { message("no provenance files found"); quit(status = 1) }
d <- rbindlist(lapply(f, fread, colClasses = c(seed = "character")), fill = TRUE)
setorder(d, seed)
fwrite(d, "/pipeline/data/mvp/conversion_summary.tsv", sep = "\t")
cat("wrote data/mvp/conversion_summary.tsv --", nrow(d), "seed(s)\n")
print(d[, .(seed, n_loci_raw, n_loci, n_duplicate_collapsed, af_max_deviation,
            causal_temp, causal_sal, causal_any, n_sites, n_individuals)])
'
exit $status
