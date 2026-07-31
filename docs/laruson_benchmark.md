# Simulation Benchmark: Láruson et al. 2022

This document is the end-to-end guide for running ADAPTOGENE on a simulated dataset with known ground truth. It covers data download, format conversion, pipeline configuration, headless execution, and evaluation on two independent axes: offset accuracy and adaptive-locus detection.

> **Status.** Phases 1-3 of `docs/laruson_automation_roadmap.md` (custom-environment source, headless SNP-set promotion, run-all driver) are implemented (commit `e2218b7`). Phase 0 (Gate G1) is confirmed — see `docs/laruson_dataset_notes.md`. Phase 4 (this dataset's converter, `benchmarks/convert_laruson.R`) and `config_LARUSON.yaml` are implemented and verified against a real archive replicate (Case 1, seed 2889863491989). **Phase 5 (the two-axis eval harness) and Shiny UI exposure of `Climate.source: custom` are explicitly deferred** — this guide covers running the pipeline end-to-end on simulated input, not yet scoring it against truth.

---

## 1. Why Simulated Data?

Empirical datasets (Arabidopsis, Populus) are the primary validation targets for ADAPTOGENE — they test the full pipeline on real genomic data. But they have a fundamental limitation: **no ground truth**. We can observe which SNPs are associated with traits or climate, but we cannot know which are truly causal, and we cannot perform controlled common-garden fitness experiments at the scale needed to validate offset predictions.

Simulated datasets provide exactly what empirical datasets cannot:

- **Known causal loci.** Every adaptive SNP (QTN) is recorded at the time of simulation. We can compute precision and recall directly.
- **Fitness measurements.** Reciprocal-transplant common garden fitness is a direct output of the SLiM model. We can test whether predicted offset negatively correlates with measured fitness.
- **Controlled architecture.** We know the number of QTNs, their effect sizes, the demographic model, and the environmental gradient. This lets us distinguish true failures (pipeline misses real adaptive signal) from expected behavior (linked-neutral loci trigger false positives — this is the population genetics, not a bug).

### The two benchmark axes

**Axis 1 — Offset accuracy.** Does the pipeline's genetic offset predict fitness? A valid offset method should produce higher offset values for populations transplanted to environments far from their home climate, and lower fitness should follow. We measure Spearman and Kendall rank correlations between per-site offset and per-garden fitness across the reciprocal-transplant matrix.

**Axis 2 — Adaptive-SNP detection.** What fraction of true QTNs does the pipeline recover (recall)? How many non-adaptive hits are produced (precision)? **Critically: does combining multiple GEA methods reduce false positives and false negatives?** The pipeline supports combining EMMAX, LFMM, and GAPIT models — this benchmark tests whether that combination is scientifically justified (Pareto-improves on any single method in the precision-recall plane).

---

## 2. Dataset: Láruson et al. 2022

### Citation

Láruson ÁJ, Lotterhos KE, Chamberland VF, Kelley JL, Sundaram M, Lind BM (2022). Comparing Measures of Adaptive Genetic Variation for Use in Conservation — A New Tool for Conservationists. *Evolutionary Applications* 15(7):1048-1063. DOI: [10.1111/eva.13354](https://doi.org/10.1111/eva.13354)

### Background

The dataset uses SLiM3 (Haller & Messer 2019) to simulate diploid individuals on a **100-deme, 10x10 grid landscape** with two orthogonal environmental optima axes and two quantitative traits (confirmed — see `docs/laruson_dataset_notes.md`). The archive ships four "Case" replicated scenarios, all multilocus/polygenic, distinguished by **landscape shape** rather than genetic architecture:

- **Case 1 (`ML_WF_Cline.slim`)** — simple monotonic cline (closest to a "1D gradient" reading)
- **Case 2 (`ML_WF_MountainRange.slim`)** — non-monotonic ("mountain") gradient
- **Case 3 / Case 4** — undocumented beyond script name in the shipped README

`seeds_Neutral.txt` / `seeds_SingleLocus.txt` are seed lists only — **the true Neutral and Single-Locus scenarios have no corresponding data files in this archive.**

Per-replicate fitness data (`_ML_WF_CG_sum_Gen.txt`) is a **population-average generation time series** (`sympatry`, `allopatry`, `local_adaptation` per ~100-generation snapshot) — **not** a home-deme x garden-deme reciprocal-transplant matrix. Building a real Axis-1 (offset-vs-fitness) evaluation will need a different construction (e.g. per-deme `P#_fit` columns from `_Freq_ML_WF.txt`) — deferred to the eval-harness phase, not attempted by the converter.

The paper's key finding: GF offset correlates negatively with fitness across most scenarios, validating the method — but with important confounders (deme size effects, nonlinear environments). This benchmark aims to test whether ADAPTOGENE reproduces that correlation **and** whether the pipeline's multi-method GEA approach identifies the underlying QTNs — once the eval harness exists.

### Why Láruson 2022 (not other benchmarks)?

- Ships in `.vcf.gz` format → direct VCF input (no allele-count-matrix conversion)
- Explicit QTN list → Axis 2 (causal-SNP detection) is directly testable
- Reciprocal-transplant fitness matrix → Axis 1 (offset vs fitness) across multiple garden scenarios
- Pilot scale (pilot-fast) — validates the harness end-to-end before the larger Lind & Lotterhos 2025 MVP benchmark

> **Future benchmark:** Lind & Lotterhos 2025 (BCO-DMO DOI: `10.26008/1912/bco-dmo.889769.1`, Zenodo: `10.5281/zenodo.7622893`) covers 2,250 simulation replicates and is the publication-grade gold standard. It requires a separate converter (allele-count-matrix format, not VCF). Not in scope for this pass.

---

## 3. Download Instructions

### 3a. Prerequisites

```bash
# Required tools (all available via nix)
nix shell nixpkgs#wget nixpkgs#bcftools nixpkgs#htslib

# Docker image must be built
docker build -t adaptogene .

# Free disk space needed: ~2 GB for archive + converted files
```

### 3b. Fetch the Dryad archive

```bash
# Create raw data directory (read-only after extraction — never modify data/laruson/)
mkdir -p data/laruson

# Download the Dryad archive
wget -O data/laruson/laruson2022_dryad.zip \
  "https://datadryad.org/stash/downloads/file_stream/doi:10.5061/dryad.x95x69pkk"

# Verify checksum (fill SHA256 after Gate G1 confirms it)
# sha256sum data/laruson/laruson2022_dryad.zip
# expected: <FILL_FROM_GATE_G1>

# Extract
cd data/laruson && unzip laruson2022_dryad.zip && cd ../..
```

> **Checkpoint.** After extraction, `data/laruson/` should contain the files listed in Section 3c. If the download URL above is stale (Dryad periodically rotates direct-download links), go to [https://datadryad.org/dataset/doi:10.5061/dryad.x95x69pkk](https://datadryad.org/dataset/doi:10.5061/dryad.x95x69pkk) and click the download button manually, or use the Dryad API: `wget "https://api.datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.x95x69pkk/download" -O laruson2022_dryad.zip`.

### 3c. File manifest (confirmed against the real archive)

The archive (1.14 GB, `doi_10_5061_dryad_x95x69pkk__v20220525.zip`) contains four simulation Case scenarios x 10 replicates per case. Per replicate, five file types are produced (filenames use `Case{N}_{SEED}`, no underscore between "Case" and the number):

| File pattern | Format | Size (approx) | Role |
|--------------|--------|---------------|------|
| `Case{N}_{SEED}.vcf.gz` | VCF.gz | ~27.7 MB | Genotype input → `laruson.vcf` (decompress only — already chr-normalized, single contig "1") |
| `Case{N}_{SEED}_ind.txt` | space-delimited | ~635 KB | 10,000 rows x 8 cols: `indID indSubpopIndex subpop phen0 phen1 opt0 opt1 fitness`. **VCF join is by row order**, not by any ID column — see `docs/laruson_dataset_notes.md` §3 |
| `Case{N}_{SEED}_causal_mutations_pos_filtered.txt` | headerless, 1 column | ~700 B, 102 rows | Exact VCF `POS` values of causal loci (verified 100% match) → `causal_loci.tsv`. No effect-size column |
| `Case{N}_{SEED}_Freq_ML_WF.txt` | space-delimited | ~157 KB, 33 rows x 613 cols | Per-generation, per-deme (`P1..P100`) fitness/frequency/phenotype/optima — not used by the converter |
| `Case{N}_{SEED}_ML_WF_CG_sum_Gen.txt` | space-delimited | ~2.5 KB, 33 rows | Population-average generation time series (`sympatry`, `allopatry`, `local_adaptation`, ...) — **not** a per-deme×garden fitness matrix |

**Plus documentation:** `README.txt` (column dictionary, no `src/` scripts shipped), `seeds_source_R.txt` (seed→Case mapping), `seeds_Neutral.txt`/`seeds_SingleLocus.txt` (seed lists only, no data files).

**Simulation cases (corrected — all four are multilocus/polygenic, distinguished by landscape shape, not architecture):**

| Case | SLiM script | Landscape |
|------|-------------|-----------|
| 1 | `ML_WF_Cline.slim` | Simple monotonic cline |
| 2 | `ML_WF_MountainRange.slim` | Non-monotonic "mountain" gradient |
| 3 | `ML_WF_Case3.slim` | Undocumented beyond script name |
| 4 | `ML_WF_Case4.slim` | Undocumented beyond script name |

**Chosen pilot replicate: Case 1, seed `2889863491989`** — simplest landscape (monotonic cline), matches this doc's original "1D gradient" framing most closely. Full details in `docs/laruson_dataset_notes.md`.

---

## 4. Data Model and Conversion

Láruson simulations use abstract environments on a **100-deme, 10x10 grid** — no real geography. `benchmarks/convert_laruson.R` bridges this to ADAPTOGENE's expected input set:

- **Valid, well-spaced lat/lon** — `latitude = row * 2.0`, `longitude = col * 3.0`, where row/col are derived per-replicate from the actual unique environmental-optima grid (not a hardcoded formula). Never used for raster extraction under `Climate.source: custom`.
- **Environment variables** — the sim's two optima axes (`opt0`, `opt1`) renamed `bio_1`, `bio_2`.
- **VCF** — decompress only (native R `gzfile()`, no external `zcat`/bash dependency); already chr-normalized (single contig `"1"`, integer `POS`).
- **Minimal GFF** — one synthetic `gene` record per causal locus (±500 bp) so gene-annotation/enrichment rules resolve without a real GFF. `GFF.go_field` left `"NULL"` — no GO terms, enrichment cleanly no-ops.

### Run the converter

```bash
docker run --user $(id -u):$(id -g) --rm --memory=20g -v $PWD:/pipeline adaptogene:latest \
  Rscript /pipeline/benchmarks/convert_laruson.R \
    data/laruson/            `# raw archive dir (read-only)` \
    data/laruson_converted/  `# output dir` \
    1                         `# CASE_N = 1 (Cline)` \
    2889863491989             `# SEED (first Case-1 seed)` \
    5                         `# TP_WINDOW_KB (placeholder — see dataset notes §4)`
```

### Converter outputs (confirmed against a real run)

| File | Description |
|------|-------------|
| `laruson.vcf` | Uncompressed VCF, unmodified otherwise |
| `metadata.tsv` | `site, sample, latitude, longitude, bio_1, bio_2` — 10,000 rows |
| `environments_present.tsv` | `site, bio_1, bio_2` — 100 rows (one per deme) |
| `environments_future.tsv` | `site, bio_1, bio_2` — present + fixed +0.5 synthetic shift. **Placeholder only** — no real "future"/garden concept exists in this static landscape; see dataset notes §4 |
| `causal_loci.tsv` | `chr, pos, category` (causal / linked_neutral / background_neutral) — **no `effect_size` column** (not present in source data), one row per VCF variant row |
| `fitness_timeseries.tsv` | Raw pass-through of `_ML_WF_CG_sum_Gen.txt` — population-average time series, archived for a future eval harness, **not consumed by the adapter-only pipeline run** |
| `laruson_minimal.gff3` | One synthetic gene record per causal locus |

---

## 5. Configuration

`config_LARUSON.yaml` (repo root) — current schema, not the stale placeholder previously drafted here. Key settings that differ from a typical run:

```yaml
Climate:
  source: custom
  predictors: "bio_1,bio_2"
  custom:
    present_table:   environments_present.tsv
    future_table:    environments_future.tsv    # single table, no per-garden directory (no garden loop today)
    columns:         "bio_1,bio_2"
    grid_resolution: 0.5
    key:             site

GEA:
  configs:
    - { method: EMMAX, adjust: bonf, threshold: '0.05' }
    - { method: LFMM,  adjust: bonf, threshold: '0.05' }

GWAS:
  configs: []                        # no phenotype data — environments run through GEA only

Maladaptation:
  methods:
    gradient_forest: { spatial_correction: without, random_model: false }
    geometric_offset: {}
  snp_sets: [laruson_gea]
```

`sNMF.k_best` is intentionally left `null` — same "review the cross-entropy plot before setting K" convention as any dataset; see `config_LARUSON.yaml`'s inline comment.

---

## 6. Headless Run

```bash
docker build -t adaptogene .   # confirm image current

# Pass 1: runs processing, prestructure, structure, gea, promote_snp_set, maladaptation.
# Will fail at mode=structure if sNMF.k_best is still null — expected.
bash benchmarks/run_benchmark.sh config_LARUSON.yaml laruson_gea

# Inspect LARUSON_results/PreStructure/plots/cross_entropy_K*.png (view it yourself —
# per project rules, image files are not opened programmatically), set sNMF.k_best in
# config_LARUSON.yaml, then re-run — Snakemake skips the already-completed steps:
bash benchmarks/run_benchmark.sh config_LARUSON.yaml laruson_gea
```

`benchmarks/run_benchmark.sh` already implements exactly this chain and already calls `promote_snp_set.R` with its real argument order — no script changes were needed for this pass. There is no reciprocal-transplant garden loop (the roadmap's original multi-garden design assumed a per-garden directory that doesn't exist in the merged `Climate.custom.future_table` implementation — it's a single table) — Phase 5 (eval harness) will need to design that separately once it exists.

---

## 7-9. Evaluation, Expected Results, Caveats — deferred to Phase 5

**Not implemented in this pass.** `benchmarks/eval_offset.R`, `benchmarks/eval_detection.R`, and `benchmarks/report.R` do not exist yet. Gate G1 findings relevant to designing them later (full detail in `docs/laruson_dataset_notes.md`):

- **No per-deme×garden fitness matrix exists in this dataset.** `_ML_WF_CG_sum_Gen.txt` is a population-average `sympatry`/`allopatry` time series, not itemized per garden. A real Axis-1 (offset-vs-fitness) design will need `_Freq_ML_WF.txt`'s per-generation, per-deme `P#_fit` columns, or a different construction entirely.
- **No `effect_size` for causal loci** — Axis-2 (precision/recall) is still buildable from `causal_loci.tsv`'s `category` column, just without effect-size-weighted analysis.
- `TP_WINDOW_KB` (5 kb default in the converter) is an unvalidated placeholder pending a real LD-decay estimate once `mode=processing`/`mode=structure` have run on the converted data.
- Abstract environments have no spatial correlation structure — `spatial_correction: without` in `config_LARUSON.yaml` is intentional, not a stopgap.
- Linked-neutral hits (within `TP_WINDOW_KB` of a causal locus) are **expected**, not false positives — report them as a separate category once the eval harness exists.
- **MVP/Lind & Lotterhos 2025 benchmark** (BCO-DMO/Zenodo, 2,250 seeds, no VCF — allele-count matrix requiring a from-scratch VCF synthesizer) remains a distinct future follow-on, not in scope here.

---

## References

- Láruson ÁJ et al. (2022). Comparing Measures of Adaptive Genetic Variation for Use in Conservation — A New Tool for Conservationists. *Evolutionary Applications* 15(7):1048-1063. [DOI: 10.1111/eva.13354](https://doi.org/10.1111/eva.13354)
- Haller BC & Messer PW (2019). SLiM 3: Forward Genetic Simulations Beyond the Wright-Fisher Model. *Molecular Biology and Evolution* 36(3):632-637. [DOI: 10.1093/molbev/msy228](https://doi.org/10.1093/molbev/msy228)
- Fitzpatrick MC & Keller SR (2015). Ecological genomics meets community-level modelling of biodiversity: mapping the genomic landscape of current and future environmental adaptation. *Ecology Letters* 18:1–11. [DOI: 10.1111/ele.12376](https://doi.org/10.1111/ele.12376)
- Gain C et al. (2023). A quantitative theory for genomic offset statistics. *Molecular Biology and Evolution* 40(6):msad140. [DOI: 10.1093/molbev/msad140](https://doi.org/10.1093/molbev/msad140)
- Lind BM & Lotterhos KE (2025). Evaluation of genomic offset predictions using individual fitness measured in common gardens. *Molecular Ecology Resources*. [DOI: 10.1111/1755-0998.14008](https://doi.org/10.1111/1755-0998.14008) *(The gold-standard framework; this benchmark follows its evaluation design)*
- Lotterhos KE (2024). Model Validation Program dataset (BCO-DMO: 10.26008/1912/bco-dmo.889769.1; Zenodo: 10.5281/zenodo.7622893) *(Future extension)*
