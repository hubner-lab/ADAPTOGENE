#!/usr/bin/env bash
# =============================================================================
# mvp_run_sweep.sh -- run the GEA parameter-cell sweep across MVP replicates.
#
#   benchmarks/mvp_run_sweep.sh [SEEDS_CSV|all] [SEED_JOBS] [CELLS]
#
# Defaults: all seeds in benchmarks/mvp_seeds.tsv, 7 seeds at a time, 5 cells.
#
# WHAT IT DOES, per seed, strictly serially:
#   upstream:  processing -> prestructure -> structure      (once; Snakemake skips if current)
#   pregea:                                                 (once; supplies the recommender arm)
#   cells:     mode=gea with config_MVP{seed}_c{i}.yaml, i = 1..CELLS
#
# After each cell the four p-value tables and the Manhattan plots are copied out to
#   benchmarks/mvp_eval/params/MVP{seed}/c{i}/
# because the pipeline writes every cell to the SAME path -- {method}_pvalues_K{k_best}.tsv
# encodes k_best, never the swept parameter (common.smk:1308). Harvesting is what turns a
# sequence of overwrites into a ladder.
#
# Cells are NOT forced. Snakemake's params rerun-trigger detects the changed rung on its
# own (verified: "Params have changed since last execution: before: 5 now: 0"), so only the
# method rules whose parameter actually moved are re-fitted. The HARVEST GATE below is what
# guards against that trigger silently failing: a cell whose harvested table is byte-identical
# to the previous cell's is a hard error, not a warning -- an identical table means the rung
# never took effect and the whole ladder would be a fiction.
#
# Seeds run in parallel; cells within a seed never do. Each {PROJECT}_results/ carries its own
# .snakemake/ (locks, metadata) -- verified: the repo-root .snakemake holds only log/ -- so
# cross-seed parallelism is safe, while two modes on ONE project would collide on its lock.
#
# Docker: /snap/bin/docker cannot start on this host (setuid blocked under nosuid+NoNewPrivs),
# so every invocation goes through a nix-provided client. See CLAUDE.local.md.
# -e USER=pipeline is required or snakemake's getpwuid() lookup aborts under --user.
# =============================================================================
set -uo pipefail

PIPELINE_ROOT="${PIPELINE_ROOT:-/mnt/data/eugene/ADAPTOGENE}"
IMAGE="${IMAGE:-adaptogene:latest}"
SEEDS_ARG="${1:-all}"
SEED_JOBS="${2:-6}"
CELLS="${3:-5}"

# BLAS_THREADS is not a tuning nicety, it is the difference between this sweep finishing
# overnight and not finishing at all. The container sees the HOST's nproc (192) regardless of
# --cpus, so OpenBLAS-pthread spawns 192 threads that thrash against the container's CPU
# share. Measured in-image on a 1500x1500 crossprod+svd, --cpus=32:
#     unset (192 threads)  52.65 s      1 thread   3.16 s
#     8 threads             1.42 s     16 threads  1.30 s   <- 40x faster than unset
# Unset, one RDA fit on seed 1232548 ran >27 min and climbed to 81 GB resident without
# finishing; the thread-local buffers are also where the memory went. Keep this pinned and
# IDENTICAL across every run in the sweep -- BLAS thread count can reorder floating-point
# summation, so mixing settings within one benchmark would be a reproducibility hole.
BLAS_THREADS="${BLAS_THREADS:-16}"
CPUS_PER_SEED="${CPUS_PER_SEED:-20}"
MEM_PER_SEED="${MEM_PER_SEED:-60g}"
SNAKE_CORES="${SNAKE_CORES:-4}"

MANIFEST="$PIPELINE_ROOT/benchmarks/mvp_seeds.tsv"
EXCLUSIONS="$PIPELINE_ROOT/benchmarks/mvp_method_exclusions.tsv"
PARAMS_DIR="$PIPELINE_ROOT/benchmarks/mvp_eval/params"
RUNLOG_DIR="$PIPELINE_ROOT/benchmarks/mvp_eval/runlogs"
METHODS=(EMMAX LFMM RDA BLINK)

mkdir -p "$PARAMS_DIR" "$RUNLOG_DIR"

DOCKER=(nix shell nixpkgs#docker-client -c docker)

if [[ "$SEEDS_ARG" == "all" ]]; then
    mapfile -t SEEDS < <(awk 'NR>1 && NF {print $1}' "$MANIFEST")
else
    IFS=',' read -r -a SEEDS <<< "$SEEDS_ARG"
fi

kbest_of() { awk -v s="$1" 'NR>1 && $1==s {print $11; exit}' "$MANIFEST"; }

snake() {   # snake <seed> <mode> <configfile> <logfile>
    "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm -e USER=pipeline \
        -e OPENBLAS_NUM_THREADS="$BLAS_THREADS" -e OMP_NUM_THREADS="$BLAS_THREADS" \
        --cpus="$CPUS_PER_SEED" --memory="$MEM_PER_SEED" \
        -v "$PIPELINE_ROOT:/pipeline" "$IMAGE" \
        snakemake -c"$SNAKE_CORES" -s Snakefile \
            --config mode="$2" --configfile "$3" --scheduler greedy \
            --rerun-incomplete \
        >> "$4" 2>&1
}

run_seed() {
    local seed="$1"
    local proj="MVP${seed}"
    local res="$PIPELINE_ROOT/${proj}_results"
    local k; k="$(kbest_of "$seed")"
    local log="$RUNLOG_DIR/${proj}.log"
    : > "$log"

    if [[ -z "$k" ]]; then echo "[$proj] FATAL: seed not in manifest" | tee -a "$log"; return 1; fi
    echo "[$proj] k_best=$k  cells=$CELLS  $(date -Is)" | tee -a "$log"

    # Clear a stale lock left by a killed run. Safe here specifically because this driver is
    # the only thing that ever runs snakemake against MVP* projects, and it serialises every
    # mode within a seed -- so a lock present at seed start is always a corpse, never a live
    # writer. Do NOT copy this into a context where a project can have a concurrent runner.
    if [[ -d "$res/.snakemake/locks" ]] && [[ -n "$(ls -A "$res/.snakemake/locks" 2>/dev/null)" ]]; then
        echo "[$proj] clearing stale snakemake lock" | tee -a "$log"
        "${DOCKER[@]}" run --user "$(id -u):$(id -g)" --rm -e USER=pipeline \
            -v "$PIPELINE_ROOT:/pipeline" "$IMAGE" \
            snakemake -s Snakefile --unlock --config mode=gea \
                --configfile "config_${proj}_c1.yaml" >> "$log" 2>&1
    fi

    # ---- upstream. Snakemake decides what is stale; MVP1231288 rebuilds by itself here
    # because the linkage-group-boundary fix changed its input VCF's mtime (MVP_README.md:183).
    # Upstream runs against the c1 CELL config, not the base config. Every cell config is
    # the base config with only the GEA: block and LDdecay.scope changed, and upstream reads
    # neither GEA params nor the per-chromosome LD-decay scope -- but the scope patch is what
    # keeps `ld_decay_run_chr` (which segfaults inside PopLDdecay on monomorphic deme x
    # chromosome subsets) out of the DAG. Using the base config here would reintroduce it.
    for mode in processing prestructure structure; do
        echo "[$proj] mode=$mode" | tee -a "$log"
        snake "$seed" "$mode" "config_${proj}_c1.yaml" "$log" \
            || { echo "[$proj] FAILED at mode=$mode" | tee -a "$log"; return 1; }
    done

    # ---- PreGEA: the hyperparameter recommender. Its ladders are fitted on the LD-pruned
    # marker set, where roughly half the causal loci are gone (seed 1232548: 30 -> 14 on the
    # temperature axis), so PreGEA NOMINATES a rung and the production cells below are what
    # actually get scored. Never mix a PreGEA recall number into the production surface.
    echo "[$proj] mode=pregea" | tee -a "$log"
    snake "$seed" pregea "config_${proj}_c1.yaml" "$log" \
        || echo "[$proj] WARNING: pregea failed, recommender arm unavailable" | tee -a "$log"

    # ---- cells
    local prev_dir=""
    for i in $(seq 1 "$CELLS"); do
        local cfg="config_${proj}_c${i}.yaml"
        local dest="$PARAMS_DIR/${proj}/c${i}"
        [[ -f "$PIPELINE_ROOT/$cfg" ]] || { echo "[$proj] FATAL: missing $cfg" | tee -a "$log"; return 1; }

        if [[ -f "$dest/BLINK_pvalues.tsv" ]]; then
            echo "[$proj] c$i already harvested, skipping" | tee -a "$log"
            prev_dir="$dest"; continue
        fi

        echo "[$proj] c$i  $(date -Is)" | tee -a "$log"
        local stamp="$res/.cell_start"
        touch "$stamp"
        snake "$seed" gea "$cfg" "$log" \
            || { echo "[$proj] FAILED at c$i" | tee -a "$log"; return 1; }

        mkdir -p "$dest/manhattan"
        for m in "${METHODS[@]}"; do
            src="$res/GEA/tables/methods/$m/${m}_pvalues_K${k}.tsv"
            if [[ ! -s "$src" ]]; then
                # A missing table is FATAL unless this (seed, method) is a declared exclusion
                # in benchmarks/mvp_method_exclusions.tsv. Without that distinction a genuine
                # silent failure would look identical to a known gap, which is exactly the
                # class of bug this whole harness exists to rule out.
                if awk -F'\t' -v s="$seed" -v m="$m" 'NR>1 && $1==s && $2==m {found=1}
                                                      END {exit !found}' "$EXCLUSIONS" 2>/dev/null; then
                    echo "[$proj] c$i $m EXCLUDED (declared in mvp_method_exclusions.tsv)" | tee -a "$log"
                    continue
                fi
                echo "[$proj] FATAL: c$i missing/empty $src" | tee -a "$log"; return 1
            fi
            cp "$src" "$dest/${m}_pvalues.tsv"
        done
        cp -r "$res/GEA/plots/manhattan/." "$dest/manhattan/" 2>/dev/null
        cp "$res/GEA/tables/selected_snps.tsv" "$dest/" 2>/dev/null

        # The Manhattans are what a user reads an operating point off, so a stale figure
        # paired with a fresh p-value table would be a silently wrong plot. Any PNG older
        # than this cell's start means Snakemake decided it did not need regenerating.
        local stale
        stale=$(find "$dest/manhattan" -name '*.png' ! -newer "$stamp" | wc -l)
        if (( stale > 0 )); then
            echo "[$proj] FATAL: c$i harvested $stale Manhattan PNG(s) older than the cell" \
                 "start — figures do not match this cell's p-values." | tee -a "$log"
            return 1
        fi

        # ---- HARVEST GATE (see header). An unchanged table means the rung did not apply.
        if [[ -n "$prev_dir" ]]; then
            for m in "${METHODS[@]}"; do
                [[ -s "$dest/${m}_pvalues.tsv" && -s "$prev_dir/${m}_pvalues.tsv" ]] || continue
                if cmp -s "$dest/${m}_pvalues.tsv" "$prev_dir/${m}_pvalues.tsv"; then
                    echo "[$proj] FATAL: c$i $m is byte-identical to the previous cell — the" \
                         "parameter did not take effect. Ladder would be a fiction; aborting." | tee -a "$log"
                    return 1
                fi
            done
        fi
        prev_dir="$dest"
    done

    echo "[$proj] DONE $(date -Is)" | tee -a "$log"
}

# ------------------------------------------------------------- bounded worker pool
echo "INFO: ${#SEEDS[@]} seeds, $SEED_JOBS at a time, $CELLS cells each"
echo "INFO: per seed --cpus=$CPUS_PER_SEED --memory=$MEM_PER_SEED, snakemake -c$SNAKE_CORES"
fail=0
for seed in "${SEEDS[@]}"; do
    while (( $(jobs -rp | wc -l) >= SEED_JOBS )); do wait -n; done
    run_seed "$seed" &
done
wait
for seed in "${SEEDS[@]}"; do
    grep -q "^\[MVP${seed}\] DONE" "$RUNLOG_DIR/MVP${seed}.log" 2>/dev/null \
        || { echo "INCOMPLETE: MVP${seed} (see $RUNLOG_DIR/MVP${seed}.log)"; fail=1; }
done
echo "INFO: sweep finished, $( [[ $fail -eq 0 ]] && echo 'all seeds complete' || echo 'SOME SEEDS INCOMPLETE' )"
exit $fail
