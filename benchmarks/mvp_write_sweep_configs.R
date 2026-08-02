#!/usr/bin/env Rscript
# mvp_write_sweep_configs.R -- emit the per-cell GEA parameter-sweep configs.
#
# Usage:  Rscript benchmarks/mvp_write_sweep_configs.R [--seeds=1232548,...] [--cells=5]
#
# WHY THIS EXISTS
# ---------------
# The pipeline rejects duplicate method names inside one GEA.configs list
# (common.smk:404), and the p-value filename is keyed on {method} + sNMF.k_best only
# (common.smk:1308) -- the swept parameter appears NOWHERE in the output path. So a
# ladder over EMMAX n_pcs / LFMM K / RDA condition_pcs / BLINK n_pcs cannot be produced
# by one run. It needs one `mode=gea` run per rung, with the p-value tables harvested
# out between runs (benchmarks/mvp_run_sweep.sh does the harvesting).
#
# Each emitted config is the seed's BASE config with ONLY the `GEA:` block replaced.
# `project_name` stays MVP{seed}, so upstream (processing/prestructure/structure) is
# computed once per seed and reused by all cells. Input.* is byte-identical across a
# seed's cells, so this does not trip the CLAUDE.md hazard about two configs sharing a
# project_name while pointing at different raw data -- here they point at the same data
# on purpose.
#
# THE LADDERS
# -----------
# k = sNMF.k_best for the seed (from the deposit's own K column, floored at 2).
#
#   n_pcs   (EMMAX n_pcs, RDA condition_pcs, BLINK n_pcs -> GAPIT PCA.total):
#           {0, round(k/2), k, 2k, min(20, 4k)}
#           c3 == k is the pipeline default (the @k_best sentinel, gea.py:43).
#           c1 == 0 is the load-bearing rung on this dataset: MVP seeds were selected
#           for R2(PC1~temp) in [0.30, 0.60], so structure correction has real work to
#           do. It also reproduces Lotterhos 2023's uncorrected RDA arm, because
#           condition_pcs=0 skips the second fit and makes rda.R's pmax a no-op.
#
#   LFMM K: {k-2, k-1, k, k+2, k+4}, clamped to [1, 20], deduplicated, padded upward
#           until 5 distinct rungs. Mirrors PreGEA's own k_offset=2 convention
#           (common.smk:614) but extends further up, because LFMM latent factors are
#           not the same quantity as PCA covariates and the pruned-set PreGEA ladder
#           tops out at k+2.
#
# Cells are a compute-amortization bundle, NOT a joint parameter setting under test:
# cell i simply sets every method to its own i-th rung so that one run yields one rung
# for each of the four methods. Each method's ladder is scored independently
# downstream; benchmarks/mvp_score_sweep.sh additionally builds mixed bundles
# (each method at its own best rung) for the combine rules.
#
# adjust/threshold are left at the base config's bonf 0.05. They only decide which
# {method}_pvalues_K{k}_sig_snps_{adjust}.tsv the pipeline writes; the sweep scores the
# raw p-value table post hoc across the whole threshold grid, so this choice does not
# constrain any reported number.

suppressPackageStartupMessages(library(data.table))

A <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
    a <- sub("^--", "", a)
    A[[gsub("-", "_", sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
}

ROOT      <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
SEEDS_TSV <- file.path(ROOT, "benchmarks/mvp_seeds.tsv")
N_CELLS   <- as.integer(if (is.null(A$cells)) 5L else A$cells)
METHODS   <- c("EMMAX", "LFMM", "RDA", "BLINK")

# anova.cca permutations, for BOTH the GEA RDA params and PreGEA.Advanced.
#
# Measured with benchmarks/probe_rda_cost.R (1000 individuals x 4000 SNPs, 16 BLAS threads):
#     rda() ordination               4.5 s
#     anova.cca full   (parallel=8)  5.4 s
#     anova.cca by=axis   (SERIAL) 112.1 s
#     anova.cca by=margin (SERIAL)  88.2 s
# The permutation tests are 98% of the cost, and the two expensive ones pass no `parallel=`
# argument (rda.R:323,328), so they stay single-threaded no matter how many cores the
# container gets. At 999 permutations that block is ~35 min on 4k SNPs, ~90 min per fit on
# the full 10-13k SNP set, and cells c2-c5 fit TWICE -- roughly 3 h per cell, which is
# exactly the stall observed.
#
# Why 99 does not change a single scored number here: the anova only feeds `RDA_anova.tsv`
# (a diagnostics output) and the `axes: auto` selection. With two predictors the constrained
# rank is 2, and rda.R caps the selection at the rank and floors it at 2 -- so K_sel is 2
# whatever the anova returns. That is the same fixed K Lotterhos 2023 used
# (`rdadapt(rdaout, 2)`, MVP_README.md:278). The p-values we score come from rdadapt over
# those axes, never from the permutation test.
#
# 99 is the pipeline's own documented minimum (gea.py: permutations min=99). The claim that
# the p-value table is unchanged is VERIFIED by byte-comparing RDA_pvalues against a
# 999-permutation run, not assumed -- see the journal.
PERMUTATIONS <- as.integer(if (is.null(A$permutations)) 99L else A$permutations)

man <- fread(SEEDS_TSV, colClasses = c(seed = "character"))
if (!is.null(A$seeds)) {
    want <- strsplit(A$seeds, ",")[[1]]
    man  <- man[seed %in% want]
    if (!nrow(man)) stop("No manifest rows matched --seeds=", A$seeds)
}

# ---------------------------------------------------------------- the two ladders
# Both return exactly N_CELLS distinct integers in [lo, 20]: dedupe first, then pad
# upward. Without the pad, k=2 would collapse c1/c2 of the LFMM ladder onto K=1 and
# that cell would silently re-run an identical fit.
pad_to <- function(v, n, lo, hi = 20L) {
    v <- sort(unique(pmin(pmax(as.integer(round(v)), lo), hi)))
    step <- 2L
    while (length(v) < n) {
        nxt <- max(v) + step
        if (nxt > hi) { step <- 1L; nxt <- max(v) + 1L }
        if (nxt > hi) stop("Cannot pad ladder to ", n, " rungs below ", hi)
        v <- sort(unique(c(v, nxt)))
    }
    if (length(v) > n) v <- v[seq_len(n)]
    v
}

npc_ladder  <- function(k) pad_to(c(0, round(k / 2), k, 2 * k, min(20, 4 * k)), N_CELLS, lo = 0L)
lfmm_ladder <- function(k) pad_to(c(k - 2, k - 1, k, k + 2, k + 4),            N_CELLS, lo = 1L)

# ------------------------------------------------------------- GEA block emitter
# Methods declared unusable for a given replicate in mvp_method_exclusions.tsv are omitted
# from that seed's GEA.configs entirely. Leaving one in and letting it fail would abort the
# whole mode=gea run for the seed; omitting it lets the other three methods through and keeps
# the gap explicit and auditable in one file rather than implicit in a stack trace.
EXCL_F <- file.path(ROOT, "benchmarks/mvp_method_exclusions.tsv")
# Read by hand rather than with fread: the file is meant to be self-documenting, so it
# carries '#' comment lines, and fread on a header-plus-comments-only body fails to resolve
# any columns at all. Comment-stripping first makes the empty case (no exclusions) behave
# like the populated one.
EXCL <- local({
    empty <- data.table(seed = character(), method = character())
    if (!file.exists(EXCL_F)) return(empty)
    ln <- grep("^\\s*(#|$)", readLines(EXCL_F, warn = FALSE), invert = TRUE, value = TRUE)
    if (length(ln) < 2L) return(empty)          # header only -> nothing excluded
    p <- tstrsplit(ln[-1], "\t", fixed = TRUE)
    data.table(seed = as.character(p[[1]]), method = as.character(p[[2]]))
})
# Argument is `s`, not `seed`: inside a data.table `i` expression a bare `seed` resolves to
# the COLUMN, so naming the argument `seed` makes the filter compare the column to itself and
# silently match every row.
methods_for <- function(s) setdiff(METHODS, EXCL[seed == s, method])

gea_block <- function(seed, cell, npc, lfk, k_best) {
    keep <- methods_for(seed)
    entry <- function(method, param, value) sprintf(
        '    - method: "%s"\n      adjust: "bonf"\n      threshold: \'0.05\'\n      params:\n        %s: %d',
        method, param, value)
    # RDA carries one extra param: see the permutations note below.
    rda_entry <- function(value) paste0(entry("RDA", "condition_pcs", value),
                                        "\n        permutations: ", PERMUTATIONS)
    paste0(
        sprintf('#-----------------------------------------------------------------------------
# GEA -- parameter-sweep cell c%d of %d for seed %s (k_best = %d).
# Generated by benchmarks/mvp_write_sweep_configs.R -- do not hand-edit; regenerate.
#
# This file differs from config_MVP%s.yaml in the GEA: block ONLY. project_name,
# Input.*, Filter.*, LD.*, sNMF.*, Climate.*, PreGEA.* and GWAS.* are untouched, so all
# cells share one {PROJECT}_results/ tree and one upstream computation.
#
# Rungs in this cell:  EMMAX n_pcs = %d | LFMM K = %d | RDA condition_pcs = %d | BLINK n_pcs = %d
# (c%d of the n_pcs ladder; the pipeline default rung is the one where n_pcs == k_best == %d.)
#
# snp_clumping_distance stays 5000: the 100 kb default would collapse an entire 50 kb
# linkage group into a single region on this simulated genome.
#-----------------------------------------------------------------------------
GEA:
  configs:
', cell, N_CELLS, seed, k_best, seed, npc, lfk, npc, npc, cell, k_best),
        paste(Filter(Negate(is.null), list(
                if ("EMMAX" %in% keep) entry("EMMAX", "n_pcs", npc),
                if ("LFMM"  %in% keep) entry("LFMM",  "K",     lfk),
                if ("RDA"   %in% keep) rda_entry(npc),
                if ("BLINK" %in% keep) entry("BLINK", "n_pcs", npc))), collapse = "\n"),
        "\n  snp_clumping_distance: 5000\n  promoter_length: 1000\n")
}

# --------------------------------------------------------------- LDdecay scope patch
# PopLDdecay SEGFAULTS (rc 139) when a subpopulation has zero usable sites on a
# chromosome. Reproduced deterministically on seed 1232218: `-SubPop deme_008` on linkage
# group 1 skips all 491 sites as non-bi-allelic/singleton, leaving 0 SNPs, and the binary
# dies instead of reporting an empty result. That is unavoidable on this dataset by
# construction -- 100 demes x 10 individuals against 50 kb linkage groups carrying only
# ~250-500 SNPs, so some deme x chromosome combinations are entirely monomorphic.
#
# `ld_decay_run_chr` is on the mode=gea critical path (common.smk:1862 puts
# O['ld_decay_table'] in the GEA targets), so the crash blocks the sweep outright.
#
# Setting scope to genome_wide removes that rule from the DAG entirely
# (structure.smk:327 gates it on the scope) while KEEPING `ld_decay_table`, which is the
# only LD-decay product any GEA rule consumes (_assoc_downstream.smk:89-111, group 'All').
# The genome-wide run pools all 20 linkage groups (~5-13k SNPs), so no 10-sample deme is
# monomorphic there -- it already completes on every seed. Per-chromosome LD-decay curves
# are a side product of mode=structure that the GEA benchmark never reads.
#
# Side benefit: it drops ~2,000 PopLDdecay invocations per seed (20 chromosomes x 101
# groups), which were pure cost on the critical path.
# PreGEA runs the same three anova.cca calls on its own RDA ladder (mode=pregea's
# pregea_rda_setup), so it inherits the same 98%-of-cost problem -- all six upstream seeds
# were observed parked there. Same reasoning applies: two predictors, rank 2, axes pinned.
patch_pregea_perms <- function(lines) {
    start <- grep("^PreGEA:", lines)
    if (!length(start)) return(lines)
    nxt <- grep("^[A-Za-z][A-Za-z0-9_]*:", lines)
    nxt <- nxt[nxt > start[1]]
    stop_at <- if (length(nxt)) nxt[1] - 1L else length(lines)
    i <- grep("^\\s*permutations:", lines[start[1]:stop_at])
    if (!length(i)) return(lines)
    j <- start[1] + i[1] - 1L
    lines[j] <- sprintf("    permutations: %d                        # was 999: anova.cca is 98%% of RDA cost and cannot move K_sel (rank 2)", PERMUTATIONS)
    lines
}

patch_lddecay <- function(lines) {
    i <- grep("^\\s*scope:", lines)
    i <- i[i > grep("^LDdecay:", lines)[1]]
    if (!length(i)) return(lines)
    lines[i[1]] <- paste0(
        '  scope: "genome_wide"                     # NOT "both": PopLDdecay segfaults on ',
        'monomorphic deme x chromosome subsets (see benchmarks/mvp_write_sweep_configs.R)')
    lines
}

# ------------------------------------------------------------------- patch a config
# Replace the [GEA: .. line before the next top-level key] span. Working on the base
# config's own text rather than re-templating the whole file keeps config_MVP{seed}.yaml
# the single source of truth for everything that is not being swept.
patch_config <- function(base_lines, new_block) {
    starts <- grep("^GEA:", base_lines)
    if (length(starts) != 1L) stop("Expected exactly one '^GEA:' line, found ", length(starts))
    s <- starts[1]
    tops <- grep("^[A-Za-z][A-Za-z0-9_]*:", base_lines)
    nxt  <- tops[tops > s]
    if (!length(nxt)) stop("No top-level key follows GEA:")
    e <- nxt[1] - 1L
    # Walk back over the comment header that belongs to the NEXT block, and over the
    # blank line separating the two, so the replacement lands cleanly.
    while (e > s && grepl("^(#|\\s*$)", base_lines[e])) e <- e - 1L
    # Same for the comment header that introduces GEA: itself.
    h <- s - 1L
    while (h > 0L && grepl("^#", base_lines[h])) h <- h - 1L
    c(base_lines[seq_len(h)],
      strsplit(new_block, "\n", fixed = TRUE)[[1]],
      base_lines[(e + 1L):length(base_lines)])
}

# ------------------------------------------------------------------------- emit
rows <- list()
for (i in seq_len(nrow(man))) {
    r      <- man[i]
    seed   <- r$seed
    k_best <- as.integer(r$k_best)
    npcs   <- npc_ladder(k_best)
    lfks   <- lfmm_ladder(k_best)

    base_f <- file.path(ROOT, sprintf("config_MVP%s.yaml", seed))
    if (!file.exists(base_f)) stop("Base config missing: ", base_f)
    base_lines <- readLines(base_f, warn = FALSE)

    for (cell in seq_len(N_CELLS)) {
        out_f <- file.path(ROOT, sprintf("config_MVP%s_c%d.yaml", seed, cell))
        writeLines(patch_pregea_perms(patch_lddecay(patch_config(base_lines,
                       gea_block(seed, cell, npcs[cell], lfks[cell], k_best)))),
                   out_f)
        for (m in methods_for(seed)) {
            param <- switch(m, EMMAX = "n_pcs", LFMM = "K", RDA = "condition_pcs", BLINK = "n_pcs")
            value <- if (m == "LFMM") lfks[cell] else npcs[cell]
            rows[[length(rows) + 1L]] <- data.table(
                seed = seed, cell = paste0("c", cell), k_best = k_best,
                method = m, param = param, value = value,
                is_default = (m == "LFMM" && value == k_best) || (m != "LFMM" && value == k_best),
                config = basename(out_f))
        }
    }
    message(sprintf("INFO: MVP%s k_best=%d  n_pcs=[%s]  LFMM_K=[%s]",
                    seed, k_best, paste(npcs, collapse = ","), paste(lfks, collapse = ",")))
}

cells <- rbindlist(rows)
out_manifest <- file.path(ROOT, "benchmarks/mvp_sweep_cells.tsv")
fwrite(cells, out_manifest, sep = "\t")
message("INFO: wrote ", out_manifest, " (", nrow(cells), " rows = ",
        uniqueN(cells$seed), " seeds x ", N_CELLS, " cells x ", length(METHODS), " methods)")
