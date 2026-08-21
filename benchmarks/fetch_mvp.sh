#!/usr/bin/env bash
# fetch_mvp.sh -- pull only the wanted seeds' members out of the MVP deposit's three tars.
#
# BCO-DMO doi:10.26008/1912/bco-dmo.889769.1 serves whole tars only -- there is no per-seed
# download. The tars are UNCOMPRESSED (members are individually gzipped), so we stream each one
# through tar and extract just the members matching our seeds. ~23.5 GB crosses the network;
# ~150 MB lands on disk.
#
# The datadocs URLs 307-redirect to presigned S3 links that expire after 1 h, so `curl -L`
# re-resolves on every invocation. Re-running is safe and cheap: a tar whose members are all
# already present is skipped.
#
# Usage:  benchmarks/fetch_mvp.sh [SEEDS_TSV] [RAW_DIR]
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SEEDS_TSV="${1:-$ROOT/benchmarks/mvp_seeds.tsv}"
RAW_DIR="${2:-$ROOT/data/mvp/raw}"
BASE="https://datadocs.bco-dmo.org/dataset/889769/file"

[[ -f "$SEEDS_TSV" ]] || { echo "ERROR: seed manifest not found: $SEEDS_TSV" >&2; exit 1; }

# Column 1 of the manifest is the seed; skip the header.
mapfile -t SEEDS < <(awk 'NR>1 {print $1}' "$SEEDS_TSV")
[[ ${#SEEDS[@]} -gt 0 ]] || { echo "ERROR: no seeds parsed from $SEEDS_TSV" >&2; exit 1; }
echo "Seeds (${#SEEDS[@]}): ${SEEDS[*]}"

# tar member dir | S3 object key | tar filename | per-seed member suffix
TARS=(
  "individuals|oAy9LYjtrl3kp4|individuals.tar|Rout_ind_subset.txt.gz"
  "genotypes|m7yVJwyCpRRJyy|genotypes.tar|Rout_Gmat_sample.txt.gz"
  "mutations|93g7482cYL8wzr|mutations.tar|Rout_muts_full.txt.gz"
)

fetch_one() {
  local dir key tar suffix
  IFS='|' read -r dir key tar suffix <<<"$1"

  mkdir -p "$RAW_DIR/$dir"

  local patterns=() missing=0
  for s in "${SEEDS[@]}"; do
    if [[ ! -s "$RAW_DIR/$dir/${s}_${suffix}" ]]; then
      patterns+=("$dir/${s}_${suffix}")
      missing=$((missing + 1))
    fi
  done

  if [[ $missing -eq 0 ]]; then
    echo "[$dir] all ${#SEEDS[@]} members already present -- skipping ($tar not downloaded)"
    return 0
  fi
  echo "[$dir] $missing/${#SEEDS[@]} members missing -- streaming $tar"

  local t0 rc
  t0=$(date +%s)
  # --ignore-command-error keeps a partial-match run from masking curl's exit status.
  curl -sSL --fail "$BASE/$key/$tar" \
    | tar -x --wildcards --no-same-owner -C "$RAW_DIR" -f - "${patterns[@]}"
  rc=$?
  local elapsed=$(( $(date +%s) - t0 ))

  # Verify by presence, not by exit code: tar returns non-zero when a stream ends after all
  # requested members were already extracted, which is a success for our purposes.
  local still=0
  for s in "${SEEDS[@]}"; do
    [[ -s "$RAW_DIR/$dir/${s}_${suffix}" ]] || still=$((still + 1))
  done
  if [[ $still -gt 0 ]]; then
    echo "[$dir] FAILED after ${elapsed}s (tar rc=$rc): $still member(s) still missing" >&2
    return 1
  fi
  echo "[$dir] OK in ${elapsed}s -- ${#SEEDS[@]}/${#SEEDS[@]} members present"
  return 0
}

status=0
pids=()
for spec in "${TARS[@]}"; do
  fetch_one "$spec" &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || status=1; done

echo "=== summary ==="
du -sh "$RAW_DIR"/* 2>/dev/null
ls -1 "$RAW_DIR"/*/ | wc -l
exit $status
