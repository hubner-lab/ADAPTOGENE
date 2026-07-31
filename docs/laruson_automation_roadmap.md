# Headless Automation Roadmap

This document specifies every pipeline change needed to run ADAPTOGENE end-to-end — from raw simulated inputs through evaluation — **without the Shiny GUI**. The primary target is the Láruson 2022 simulation benchmark; the changes are designed to be general enough to support any future benchmark dataset.

> **Status update.** This was originally written as a pure specification (see note below) — it has since been partially implemented, in two passes:
> - **Phases 1-3 implemented** (commit `e2218b7`): `Climate.source: custom`, `scripts/stage_custom_climate.R`, `scripts/promote_snp_set.R`, `benchmarks/run_benchmark.sh`. Two of this section's own assumptions turned out **stale** relative to what was actually built — see the callouts inline below (Phase 1c's per-garden directory, Phase 2's `promote_snp_set` argument order, the config template's renamed keys).
> - **Phase 0 (Gate G1) confirmed and Phase 4 (`benchmarks/convert_laruson.R`) implemented**, verified against a real archive replicate (Case 1, seed `2889863491989`) — see `docs/laruson_dataset_notes.md` for confirmed facts and `docs/laruson_benchmark.md` for the current run guide. Several Phase 4 assumptions below (fitness matrix shape, causal-loci schema) turned out different from what's marked "FILL" here.
> - **Phase 5 (eval harness) and Shiny UI exposure of `Climate.source`/`Climate.custom.*` remain unimplemented** — explicitly deferred, not started.
>
> The phase-by-phase spec below is left as originally written (historical record of the design); do not treat its "FILL"/placeholder values as current — `docs/laruson_dataset_notes.md` supersedes them.

> **This is a specification, not an implementation.** File paths, function signatures, acceptance criteria, and config keys are fully specified here. No pipeline or eval code has been written yet. Implement phases sequentially: Phase 0 → 1 → 2 → 3 unlock a headless run; Phase 4 feeds inputs; Phase 5 evaluates outputs.

---

## Background: Why Headless Is Blocked Today

Two structural blockers prevent running the pipeline without the GUI:

### Blocker 1 — Climate hardwired to WorldClim/CMIP6

`download_climate_present.R` builds a `SpatVector` from `latitude`/`longitude` columns and calls `terra::extract()` against WorldClim rasters; it issues a hard `stop()` if coordinates produce NA values (L166-183). The `download_climate_future_model.R` and `merge_climate_future.R` scripts follow the same pattern, querying CMIP6 rasters. Simulated datasets (including Láruson 2022) use **abstract environments** — no real geography — so real-valued lat/lon that happen to be in the ocean, or negative coordinates outside the WorldClim extent, produce NA extractions and abort.

Affected Snakemake rules:
- `download_climate_present` (`workflow/rules/prestructure.smk` or `structure.smk`)
- `download_climate_future_model` (`workflow/rules/maladaptation.smk:15-30`)
- `merge_climate_future` (`workflow/rules/maladaptation.smk:34-56`)

### Blocker 2 — Maladaptation SNP set is a Shiny artifact

Both `gradient_forest_adaptive` and `geometric_offset` rules (`maladaptation.smk:76-99, 155-190`) consume `sigsnps = lambda wc: snp_set_file(wc.run_label)` → `_intermediate/snp_sets/<name>/selected_snps.tsv`. This file is written by the Shiny GEA tab (user clicks "Save SNP set" after reviewing results). `resolve_active_snp_sets()` in `workflow/rules/common.smk:786-822` raises an error if no snp_sets directory exists. No Snakemake rule produces this file.

### Gap 3 — No evaluation-vs-truth code

`scripts/compare_offsets.R` computes model-vs-model concordance (Spearman/Kendall/ExDet/Dutilleul/Kendall-W). Its Shiny banner explicitly states "No fitness ground truth exists for these offsets." Both axes of the benchmark — offset accuracy vs fitness, and causal-SNP detection rate — are net-new science code.

### Gap 4 — No sequential run-all driver

The pipeline has separate Docker invocations per mode. No script chains them sequentially, handles the reciprocal-transplant multi-garden loop (each garden = one future-climate table), or enforces the "never run modes in parallel" constraint.

---

## Phase 0 — Gate G1: Deep Research Pass (prerequisite)

**Goal.** Confirm dataset details before writing the converter and evaluation scripts. Attempting to implement Phases 4-5 without verifying the Láruson archive format risks building for the wrong structure.

**Action.** Run `/deep-research` on the Láruson 2022 Dryad archive (DOI: `10.5061/dryad.x95x69pkk`) and the paper itself (DOI: `10.1111/eva.13354`). Specifically confirm:

1. The exact file manifest (names, formats, whether VCF is per-chromosome or whole-genome)
2. Fitness structure: is it a reciprocal-transplant home × garden matrix (rows = populations, columns = gardens) or a scalar per population? If matrix, how many gardens?
3. Causal/adaptive loci: is there an explicit truth file? If yes, what are its columns? Are loci categorized (causal, linked-neutral, background-neutral)?
4. Simulation architecture: number of linkage groups / chromosomes, sites per LG, and approximate LD block size — this sets the TP-matching window for `eval_detection.R`
5. Number of populations/demes and samples per deme

**Output.** A short grounded note in `docs/laruson_dataset_notes.md` with the confirmed values. Use this to fill the manifest table in `docs/laruson_benchmark.md` and to parameterize the converter and eval scripts.

**Gate.** Do not proceed to Phases 4 or 5 until the manifest, fitness structure, and causal-loci format are confirmed.

---

## Phase 1 — Custom-Environment Source (Critical Path)

**Goal.** Let the pipeline accept user-supplied environment tables instead of downloading from WorldClim/CMIP6. With `Climate.source: custom`, all downstream rules (GEA, maladaptation) receive correctly formatted environment files without any geographic coordinate query.

### 1a. Config changes — `config_LARUSON.yaml` + `workflow/rules/common.smk`

Add three new config keys under the `Climate` group (validated in `common.smk`):

```yaml
Climate:
  source: custom          # worldclim (default) | custom
  present_file: data/laruson/environments_present.tsv   # abs or relative to pipeline root
  future_dir:   data/laruson/environments_future/       # directory of per-garden TSVs
```

Parse in `common.smk` near other `Climate.*` reads. `CLIMATE_SOURCE = _cfg('Climate', 'source', 'worldclim')`. Guard all existing WorldClim logic behind `if CLIMATE_SOURCE == 'worldclim':`. Emit a clear error if `source: custom` is set but `present_file` or `future_dir` are missing.

### 1b. New script — `scripts/stage_custom_climate.R`

**Purpose.** Validate and reformat user-supplied environment tables into the exact file layout that downstream rules expect. Runs once at the start of `gea` and `maladaptation` modes when `Climate.source: custom`.

**Args (positional):**
```
1  PRESENT_FILE    — user-supplied TSV: site, bio_1..bio_N (unscaled)
2  FUTURE_DIR      — directory of per-garden TSVs: garden_{id}.tsv, same cols as PRESENT_FILE
3  SAMPLES         — metadata.tsv (col 1 = site, col 3 = latitude, col 4 = longitude)
4  PREDICTORS      — comma-separated predictor names (subset of bio_N cols to use)
5  OUT_PRES_SITE   — output: climate_present_site.tsv (unscaled, sites only)
6  OUT_PRES_SCALED — output: climate_present_site_scaled.tsv (z-scaled)
7  OUT_PRES_ALL    — output: climate_present_all.tsv (with ID col, for raster cells)
8  OUT_FUT_SITE    — output: climate_future_site.tsv (default/first garden, or averaged)
9  OUT_FUT_ALL     — output: climate_future_all.tsv
10 OUT_RASTER      — output: dummy present raster .tif (n_sites × 1 grid, for spatial template)
11 GARDEN_ID       — optional; if supplied, selects a specific garden from FUTURE_DIR
```

**Spec — what the script must do:**
- Read `PRESENT_FILE`. Validate columns: `site` + all `PREDICTORS` present, no NA in predictor columns.
- Join to `SAMPLES` by `site` to get `latitude`/`longitude` for abstract grid (used only for piemap rendering, never for raster extraction).
- Emit `OUT_PRES_SITE`: `site, bio_1..bio_N` rows ordered by sample order in SAMPLES.
- Emit `OUT_PRES_SCALED`: same structure, each predictor z-scaled (center+scale over sites).
- Emit `OUT_PRES_ALL`: synthesize an abstract "landscape" grid — a regular lat/lon grid covering `[min_lon-1, max_lon+1] × [min_lat-1, max_lat+1]` at 0.5-degree resolution, extract environment values by nearest-site assignment. Columns: `ID` (1-based integer), `bio_1..bio_N`. Required by `geometric_offset.R` args 7-8 and `gradient_forest_offset.R`.
- Emit `OUT_RASTER`: a dummy single-band `SpatRaster` (`terra::rast()`) covering the abstract grid extent at 0.5-degree resolution, with NA values. Required by `geometric_offset.R` arg 9 and `gradient_forest_offset.R` arg 5 as a spatial template; the raster values are never extracted for site offsets.
- If `GARDEN_ID` is supplied: read `FUTURE_DIR/garden_{GARDEN_ID}.tsv`, emit `OUT_FUT_SITE` + `OUT_FUT_ALL` for that garden. If not supplied: average all garden TSVs in `FUTURE_DIR` and emit the average.

**colClasses safety.** Use `fread(..., colClasses = c("site" = "character"))` on all TSVs (per `fread` rule in CLAUDE.md).

### 1c. Snakemake rule branching

In the rules that currently download WorldClim (`download_climate_present`, `download_climate_future_model`, `merge_climate_future`), add a conditional branch:

**Pattern (pseudocode for each rule):**
```python
if CLIMATE_SOURCE == 'custom':
    rule <original_name>:
        input:  [user-supplied files]
        output: [same output paths as the WorldClim version]
        shell: "Rscript /pipeline/scripts/stage_custom_climate.R ..."
else:
    rule <original_name>:
        # existing WorldClim download logic unchanged
```

Because Snakemake chooses rules at parse time, use `if/else` at module level (same pattern as the existing `if MODE ==` blocks). Alternatively, use a single rule with a `params.source` branch in the shell — whichever is cleaner.

For the reciprocal-transplant loop (N gardens = N offset runs), the future site/all tables must be swappable per garden. The approach: `stage_custom_climate.R` accepts `GARDEN_ID`; `run_benchmark.sh` (Phase 3) invokes the pipeline N times with different `GARDEN_ID` injected via `--config garden_id=<id>`. The Snakemake rules pass `{config[garden_id]}` as an arg.

> **As implemented (stale note above):** `Climate.custom.future_table` (`common.smk:143`) is a **single file path**, not a directory, and the real `scripts/stage_custom_climate.R` has no `GARDEN_ID` argument. There is no per-garden loop today — the merged implementation stages exactly one future scenario. A real reciprocal-transplant Axis-1 evaluation (Phase 5) will need to add that looping mechanism when it's built; this pass (`docs/laruson_benchmark.md`) works around it with a single placeholder future table.

**Acceptance criteria:**
- `mode=gea` and `mode=maladaptation` run to completion with `Climate.source: custom` on abstract environments, with zero network activity and no NA-coordinate abort.
- Present env tables have the same column structure as the WorldClim-produced versions (validated by reading one WorldClim output and one custom output side-by-side for column types).
- All site offsets in `genetic_offset_site.tsv` are non-NA.

---

## Phase 2 — Headless SNP-Set Promotion

**Goal.** Let pipeline-produced GEA sig SNPs flow into the maladaptation module without Shiny intervention.

### 2a. New script — `scripts/promote_snp_set.R`

**Purpose.** Copy the combined `selected_snps.tsv` from the GEA output into the `_intermediate/snp_sets/<name>/` directory that `resolve_active_snp_sets()` expects.

**Args (positional):**
```
1  GEA_SELECTED_SNPS  — GEA/tables/selected_snps.tsv (pipeline output, has SNPID chr:pos col)
2  OUT_DIR            — _intermediate/snp_sets/<set_name>/  (parent dir for the snp set)
3  SET_NAME           — name of the snp set (matching Maladaptation snp_sets config entry)
```

**Spec:**
- Read `GEA_SELECTED_SNPS` with `fread(..., colClasses = c("SNPID" = "character"))`.
- Validate `SNPID` column exists and is chr:pos format (regex `^[^:]+:[0-9]+$`).
- Write to `OUT_DIR/selected_snps.tsv` (same schema; the colClasses rule above ensures character SNPID).
- Exit 0; log number of SNPs promoted.

**Args note.** `geometric_offset.R:102-104` reads `CANDIDATE_SNPS` as `fread(CANDIDATE_SNPS, header = TRUE, colClasses = c('SNPID' = 'character'))` and accesses `sig_dt$SNPID`. The promoted file must have a header with column `SNPID`.

### 2b. Snakemake rule — `promote_snp_set`

Add to `workflow/rules/gea.smk` (or a new `workflow/rules/headless.smk`):

```python
# Guard: only active when headless_snp_set is set in config
if config.get('headless_snp_set'):
    rule promote_snp_set:
        input:  O['gea_selected_snps']     # GEA/tables/selected_snps.tsv
        output: W['snp_set_promoted']      # _intermediate/snp_sets/{set_name}/selected_snps.tsv
        params:
            set_name = config['headless_snp_set'],
            out_dir  = lambda wc: str(INTER / 'snp_sets' / config['headless_snp_set'])
        shell:
            "Rscript /pipeline/scripts/promote_snp_set.R "
            "{input} {params.out_dir} {params.set_name} > {log} 2>&1"
```

`W['snp_set_promoted']` = `INTER / 'snp_sets' / config['headless_snp_set'] / 'selected_snps.tsv'` — this is exactly the path that `snp_set_file(name)` returns in `common.smk`.

### 2c. Config entry in `config_LARUSON.yaml`

```yaml
headless_snp_set: laruson_gea_combined   # triggers the promote_snp_set rule

Maladaptation:
  methods:
    gradient_forest:
      run_label: laruson_gea_combined    # must match headless_snp_set name
      ...
    geometric_offset:
      run_label: laruson_gea_combined    # must match headless_snp_set name
      ...
  snp_sets:
    - laruson_gea_combined
```

**Acceptance criteria:**
- `mode=maladaptation` resolves the SNP set from `_intermediate/snp_sets/laruson_gea_combined/selected_snps.tsv` without Shiny.
- `resolve_active_snp_sets()` finds the set and does not error.
- Both GF and geometric offset rules receive the promoted SNPID list and produce non-empty offset tables.

> **As implemented (stale note above):** the real `scripts/promote_snp_set.R` signature is `SOURCE_SELECTED_SNPS SET_NAME SNP_SETS_DIR [SOURCE_MODULE]` (argument order differs from §2a above) and it also writes/upserts `manifest.json` (Shiny-compatible), not just the TSV. `benchmarks/run_benchmark.sh` already calls it correctly — no changes needed.

---

## Phase 3 — Run-All Driver

**Goal.** A single shell script that runs all pipeline modes sequentially (fail-fast), including the reciprocal-transplant garden loop for maladaptation.

### Script — `benchmarks/run_benchmark.sh`

**Full spec:**

```bash
#!/usr/bin/env bash
# run_benchmark.sh — headless ADAPTOGENE benchmark on simulated data
# Usage: bash benchmarks/run_benchmark.sh [--config config_LARUSON.yaml] [--cores 4]
#
# Environment:
#   PIPELINE_DIR   — absolute path to the pipeline checkout (default: $PWD)
#   GARDEN_IDS     — space-separated list of garden IDs for the fitness matrix loop
#                    (read from config_LARUSON.yaml if not set)
```

**Execution sequence:**
1. `mode=processing` — filter, impute, normalize
2. `mode=prestructure` — PCA, cross-entropy
3. `mode=structure` — sNMF K-best, piemaps
4. `mode=gea` — GEA methods, sig SNPs, regions, genes (environments → GEA, no GWAS)
5. `mode=gea` post-step: `promote_snp_set` rule (headless SNP-set promotion, triggered by `headless_snp_set` config key)
6. For each garden in `GARDEN_IDS` (reciprocal-transplant loop):
   - Inject `--config garden_id=<id>` so `merge_climate_future` / `stage_custom_climate.R` selects that garden's future table
   - `mode=maladaptation`
   - Rename/copy offset output to a per-garden subdirectory: `Maladaptation/tables/<method>/<run_label>_nospatial/garden_{id}/`
7. `benchmarks/eval_offset.R` — Axis 1 (offset vs fitness)
8. `benchmarks/eval_detection.R` — Axis 2 (causal-SNP detection)
9. `benchmarks/report.R` — produce `benchmark_eval.tsv` + plots

**Docker command template** (reuse exactly from CLAUDE.md):
```bash
docker run --user $(id -u):$(id -g) --rm --memory=20g \
  -v ${PIPELINE_DIR}:/pipeline adaptogene:latest \
  snakemake -c${CORES:-4} -s Snakefile \
  --config mode=${MODE} garden_id=${GARDEN_ID:-""} headless_snp_set=laruson_gea_combined \
  --configfile ${CONFIG:-config_LARUSON.yaml} \
  --scheduler greedy
```

**Fail-fast.** Use `set -euo pipefail` at the top. Each mode must complete before the next starts. If any mode fails, the script prints the mode name + log path and exits non-zero.

**Acceptance criteria:**
- One command on a machine with Docker takes raw Láruson inputs (after conversion) through all modes to evaluation outputs, with no manual steps and no GUI.
- Garden loop produces N per-garden offset tables, one per garden in `fitness.tsv`.
- Failure in any mode surfaces clearly (mode name + log path in stderr).

---

## Phase 4 — Láruson Converter (Dataset-Specific)

**Goal.** Transform the raw Láruson Dryad archive into the exact ADAPTOGENE input set.

> **Prerequisite: Gate G1 (Phase 0) must be complete.** The exact file names, fitness structure (matrix vs scalar), and causal-loci format must be confirmed before implementing this script.

### Script — `benchmarks/convert_laruson.R`

**Input.** All files from the Dryad archive extracted into `data/laruson/` (read-only raw data). The archive has 4 Case scenarios × 10 replicates; the converter targets one Case+replicate at a time.

**Args (positional):**
```
1  RAW_DIR        — data/laruson/ (archive root, read-only)
2  OUT_DIR        — data/laruson_converted/ (output root)
3  CASE_N         — integer Case number (e.g. 3 for polygenic — confirm from README.txt)
4  SEED           — replicate seed string (from seeds_source_R.txt)
5  TP_WINDOW_KB   — half-window for linked-neutral classification (fill from Gate G1 LG architecture)
```

The input files for Case N, seed S are:
- `{RAW_DIR}/Case_{N}_{S}.vcf.gz` — genotypes
- `{RAW_DIR}/Case_{N}_{S}_causal_mutations_pos_filtered.txt` — QTN positions (explicitly listed)
- `{RAW_DIR}/Case_{N}_{S}_ML_WF_CG_sum_Gen.txt` — summary stats by generation including Common Garden (CG) fitness data
- `{RAW_DIR}/README.txt` — column dictionary (read first; extract env var names + fitness column layout)

**Outputs** (written to `{OUT_DIR}/`; never overwrite `data/laruson/`):

| Output file | Description |
|-------------|-------------|
| `laruson.vcf` | Uncompressed VCF with normalized chr names (no "chr" prefix, integer chr IDs) |
| `metadata.tsv` | `site, sample, latitude, longitude, bio_1..bio_N` — lat/lon = abstract grid (evenly spaced, valid numeric); bio cols = sim environmental variables renamed `bio_1..bio_N` |
| `environments_present.tsv` | `site, bio_1..bio_N` — present environment at each deme |
| `environments_future/garden_{id}.tsv` | Per garden: `site, bio_1..bio_N` — target environment for that garden (CG location) |
| `causal_loci.tsv` | `chr, pos, effect_size, category` — category ∈ `causal|linked_neutral|background_neutral` |
| `fitness.tsv` | `site, garden, fitness` (long format) — one row per (home_deme, garden) pair, extracted from `_ML_WF_CG_sum_Gen.txt` |
| `laruson_minimal.gff3` | One synthetic gene record per causal locus: `chr \t . \t gene \t (pos-500) \t (pos+500) \t . \t . \t . \t ID=gene_{chr}_{pos};Name=causal_{chr}_{pos}` — ensures gene/enrichment rules resolve without errors |

**Spec — key normalization steps:**
1. **VCF.** Decompress `.vcf.gz` (`system("bgzip -d ...")`). Strip "chr" prefix from CHROM column if present (`gsub("^chr", "", CHROM)`). Validate all POS are integers (not scientific notation).
2. **metadata.tsv.** One row per sample. Extract sample-to-deme mapping from VCF header or `ind.txt`. Use abstract lat/lon (assign each deme a grid cell: `lat = deme_id * 2.0, lon = deme_id * 3.0` or similar uniform grid). Rename simulation env variables to `bio_1..bio_N` matching the `Climate.predictors` list in `config_LARUSON.yaml` (column names from `README.txt`).
3. **causal_loci.tsv.** Parse `causal_mutations_pos_filtered.txt` (column names from `README.txt`). Assign category: SNPs explicitly listed as QTNs = `causal`; SNPs within ±[TP_WINDOW_KB] of a causal locus = `linked_neutral`; all others = `background_neutral`. TP_WINDOW_KB is parameterizable (default = LD-block size from Gate G1).
4. **fitness.tsv.** Extract from `_ML_WF_CG_sum_Gen.txt`. The "CG" component (Common Garden summary) contains fitness by generation for each transplant scenario. Column layout must be verified from `README.txt`. Pivot to long format: `site = home_deme, garden = garden_id, fitness = value`. Validate no NA fitness values.
5. **environments_future/garden_{id}.tsv.** Each Common Garden location has a specific environment (the "future" scenario for offset calculation). Extract these from the CG metadata (in `README.txt`) or from simulation parameters. One TSV per garden, same column layout as `environments_present.tsv`.
6. **GFF.** Write with TAB separators, standard GFF3 header, integer coordinates (no scientific notation — per CLAUDE.md GFF coord rule).

**`colClasses` safety.** Use `fread(..., colClasses = c("site" = "character", "sample" = "character"))` everywhere.

**Acceptance criteria:**
- `laruson.vcf` passes `vcfR::read.vcfR()` without errors.
- `metadata.tsv` has no NA lat/lon and no NA bio columns.
- `causal_loci.tsv` has at least 1 row per category; all positions are integers.
- `fitness.tsv` has N_demes × N_gardens rows in long format.
- The converter outputs pass ADAPTOGENE's early normalization (processing mode completes on them without error).

---

## Phase 5 — Two-Axis Evaluation Harness

**Goal.** Quantify both axes of benchmark quality: (1) does genetic offset predict fitness? (2) does the pipeline detect causal loci, and does combining methods help?

### Script — `benchmarks/eval_offset.R` (Axis 1)

**Purpose.** Per-garden, per-method correlation between pipeline-predicted genetic offset and simulated fitness.

**Args (positional):**
```
1  OFFSET_DIR      — Maladaptation/tables/ (pipeline root for offset outputs)
2  FITNESS_FILE    — data/laruson_converted/fitness.tsv
3  METADATA        — data/laruson_converted/metadata.tsv (maps sample → site)
4  METHODS         — comma-separated: gradient_forest,geometric_offset
5  RUN_LABEL       — snp_set name used as run_label in pipeline (e.g. laruson_gea_combined)
6  OUT_TSV         — benchmark_eval.tsv (appended to, long format)
7  OUT_PLOTS_DIR   — directory for per-garden offset-vs-fitness scatter plots
```

**Spec:**
- For each method in METHODS:
  - For each garden in FITNESS_FILE:
    - Load `OFFSET_DIR/<method>/<run_label>_nospatial/garden_{garden_id}/genetic_offset_site.tsv` (site-level offset).
    - Load fitness for that garden from FITNESS_FILE (filter `garden == garden_id`).
    - Join on `site`. Assert no NA after join.
    - Compute Spearman ρ and Kendall τ: `cor.test(offset, fitness, method = "spearman")` and `"kendall"`.
    - Expected direction: negative (higher offset → lower fitness).
  - Aggregate across gardens: median ρ, median τ, N_gardens where ρ < 0.
- Write to `OUT_TSV`: one row per (method, garden) + one summary row per method. Columns: `axis, method, garden, metric, value`.
- For each garden: produce a scatter plot `OUT_PLOTS_DIR/<method>_garden_{id}_offset_vs_fitness.png` (x = offset, y = fitness, loess curve, Spearman ρ in subtitle).

**Reuse.** The loess scatter style mirrors `scripts/compare_offsets.R`'s plotting patterns. Use `ggplot2` with the pipeline's standard theme (or a minimal base-R equivalent for the benchmark scripts — no bslib or Shiny dependency).

**Acceptance criterion.** At least one method shows Spearman ρ < 0 (p < 0.05) for the majority of gardens (>50%). If not, flag in `benchmark_eval.tsv` with `interpretation: offset_does_not_predict_fitness` and check the garden-id loop alignment (the most likely source of error is a mis-joined fitness garden column).

---

### Script — `benchmarks/eval_detection.R` (Axis 2)

**Purpose.** For each GEA method and each combine strategy, compute precision / recall / F1 against the causal-loci truth table.

**Args (positional):**
```
1  SIG_SNPS_DIR    — GEA/tables/methods/ (per-method significant SNP files)
2  COMBINED_SNPS   — GEA/tables/selected_snps.tsv (combined across methods)
3  CAUSAL_LOCI     — data/laruson_converted/causal_loci.tsv
4  METHODS         — comma-separated GEA method names (e.g. EMMAX,LFMM,GLM)
5  TP_WINDOW_KB    — half-window in kb for TP matching (filled from Gate G1 LG architecture)
6  OUT_TSV         — benchmark_eval.tsv (appended to, same long-format file as eval_offset.R)
7  OUT_PLOTS_DIR   — directory for precision-recall bar plots and confusion heatmaps
```

**TP matching rule (window/LD-based, non-exact).** A called SNP `c` is a True Positive for causal locus `q` if `|c.pos - q.pos| ≤ TP_WINDOW_KB × 1000` AND `c.chr == q.chr`. One called SNP can match at most one causal locus (first-match priority). All causal loci not matched by any called SNP = False Negatives.

**Locus categories** (from `causal_loci.tsv.category`). Report separately:
- `causal` — true QTNs; TP only from these
- `linked_neutral` — physically linked to QTNs; hits here are **expected** (not pipeline bugs); report as "expected-linked hits" not as FP
- `background_neutral` — unlinked neutral; true FP (these inflate FP count)

**Per-configuration output:**
- Configurations = each method individually + each combine strategy (single, Sum, Overlap)
- For each config: TP, FP (background-neutral hits), FN, expected-linked hits, precision = TP/(TP+FP), recall = TP/(TP+FN), F1 = 2·P·R/(P+R)
- Number of called SNPs total

**Multi-method comparison.** The goal is to show that combining methods (Sum or Overlap) Pareto-dominates the best single method on the FP/FN trade-off. Table: `config, precision, recall, F1, FP_background, expected_linked`. Plot: a 2D precision-recall scatter with each config as a point, connected to the Pareto frontier.

**Write to `OUT_TSV`.** Columns: `axis, config, category, metric, value`. Append to the same `benchmark_eval.tsv` produced by `eval_offset.R`.

**Acceptance criterion.** Combined strategies (Sum or Overlap) should not be dominated by any single method on both precision and recall simultaneously. If a single method Pareto-dominates all combinations, flag in the TSV and note it in `docs/laruson_benchmark.md` — this would be a valid scientific finding (the combine step adds no value on this simulation architecture).

---

### Script — `benchmarks/report.R`

**Purpose.** Assemble final `benchmark_eval.tsv` + summary plots for the manuscript.

**Args:** `EVAL_TSV, OUT_PLOTS_DIR, OUT_SUMMARY_TSV`

**Spec:**
- Read `EVAL_TSV`. Produce:
  - Figure 1: Offset vs fitness — one panel per method, scatter + loess + ρ annotation, gardens averaged.
  - Figure 2: Precision-recall — bar chart or scatter, one group per config (methods + combine strategies).
  - Figure 3 (optional): Confusion breakdown — stacked bar showing TP / expected-linked / FP / FN per config.
- Write `OUT_SUMMARY_TSV` (one row per axis/method, human-readable): `axis, method_or_config, metric, value, interpretation`.

---

## Configuration Template — `config_LARUSON.yaml`

Full config template (fill concrete values after Gate G1 and Phase 4):

```yaml
project_name: LARUSON
cpu: 4

Input:
  dir:      data/laruson_converted/
  vcf:      laruson.vcf
  metadata: metadata.tsv
  gff:      laruson_minimal.gff3

Filter:
  maf:         0.05
  snp_miss:    0.2
  sample_miss: 0.5

LD:
  window: 200
  step:   50
  r2:     0.2

sNMF:
  k_start:  1
  k_end:    5
  k_best:   <FILL from structure run>
  ploidy:   2
  repeats:  10

Climate:
  source:       custom
  present_file: data/laruson_converted/environments_present.tsv
  future_dir:   data/laruson_converted/environments_future/
  predictors:   bio_1,bio_2         # FILL after Gate G1 — number of sim env vars

Map:
  climate_extent: [-180, 180, -90, 90]   # global bbox for abstract coords
  gap:            0.1
  resolution:     10
  zoom_extent:    []

GFF:
  feature:    gene
  gene_name:  Name
  biotype:    ""
  go_field:   ""             # no GO terms in synthetic GFF

GEA:
  configs:
    - method:     EMMAX
      adjust:     bonferroni
      threshold:  0.05
    - method:     LFMM
      adjust:     bonferroni
      threshold:  0.05
  combine_method:      Sum
  combine_gap:         50000
  sig_snp_distance:    50000
  region_distance:     100000
  top_regions:         20
  promoter_length:     2000
  scattermore_threshold: 30000

GWAS:
  configs: []              # no phenotype GWAS — environments run through GEA

Enrichment:
  top_terms:       10
  plot_width:      8
  plot_height:     6
  cnet_label:      gene_id
  top_plot_regions: 0

Future:
  ssp:    245              # placeholder; not used when Climate.source: custom
  year:   2070
  models: [ACCESS-CM2]    # placeholder; not used when Climate.source: custom

Maladaptation:
  methods:
    gradient_forest:
      run_label:         laruson_gea_combined
      ntree:             500
      cor_threshold:     0.5
      spatial_correction: false
      random_model:      false
    geometric_offset:
      run_label:         laruson_gea_combined
      scale:             true
  snp_sets:
    - laruson_gea_combined

headless_snp_set: laruson_gea_combined    # triggers promote_snp_set rule
```

---

## File Map Summary

| New / Changed file | Change type | Phase |
|--------------------|-------------|-------|
| `workflow/rules/common.smk` | Add `CLIMATE_SOURCE`, `CLIMATE_PRESENT_FILE`, `CLIMATE_FUTURE_DIR` parsing; guard WorldClim logic | 1a |
| `workflow/rules/prestructure.smk` (or structure.smk) | Branch `download_climate_present` on `CLIMATE_SOURCE` | 1c |
| `workflow/rules/maladaptation.smk` | Branch `download_climate_future_model` + `merge_climate_future` on `CLIMATE_SOURCE`; pass `garden_id` | 1c |
| `scripts/stage_custom_climate.R` | **New** — validate + reformat user env tables → expected output paths | 1b |
| `workflow/rules/gea.smk` (or headless.smk) | Add `promote_snp_set` rule, guarded by `config.get('headless_snp_set')` | 2b |
| `scripts/promote_snp_set.R` | **New** — copy GEA selected_snps.tsv into snp_sets dir with correct SNPID column | 2a |
| `config_LARUSON.yaml` | **New** — full benchmark config | 1a, 2c |
| `benchmarks/run_benchmark.sh` | **New** — sequential multi-mode driver with garden loop | 3 |
| `benchmarks/convert_laruson.R` | **New** — Dryad → ADAPTOGENE inputs + truth tables | 4 |
| `benchmarks/eval_offset.R` | **New** — Axis 1: offset vs fitness per garden per method | 5 |
| `benchmarks/eval_detection.R` | **New** — Axis 2: causal-SNP TP/FP/FN per method and combine strategy | 5 |
| `benchmarks/report.R` | **New** — assemble benchmark_eval.tsv + summary figures | 5 |
| `docs/laruson_dataset_notes.md` | **New** — Gate G1 confirmed dataset facts | 0 |
| `docs/laruson_benchmark.md` | **New** — end-to-end benchmarking guide (see companion doc) | — |
| `benchmarks/TESTING_PLAN.md` | Mode name corrections (stale: structure_K/association/overlapping → structure/gea/gea_x_gwas) | — |

---

## Sequencing Diagram

```
Phase 0 (G1 research)
    │
    ├──→ Phase 1 (custom env source) ──→ Phase 3 (run-all driver) ──→ END-TO-END RUN
    │
    └──→ Phase 2 (SNP-set promotion) ──→ (feeds into Phase 3's maladaptation step)
    │
    └──→ Phase 4 (converter) ──→ (produces inputs consumed by Phase 3)
    │
    └──→ Phase 5 (eval harness) ──→ (consumes pipeline outputs from Phase 3)
```

Phase 1 is the critical path — Phases 2, 3, and 5 all depend on it. Phase 4 (converter) can be developed in parallel with Phases 1-3 as long as Gate G1 (Phase 0) is satisfied.
