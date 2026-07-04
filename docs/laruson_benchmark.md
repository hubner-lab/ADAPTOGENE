# Simulation Benchmark: Láruson et al. 2022

This document is the end-to-end guide for running ADAPTOGENE on a simulated dataset with known ground truth. It covers data download, format conversion, pipeline configuration, headless execution, and evaluation on two independent axes: offset accuracy and adaptive-locus detection.

> **This benchmark requires pipeline changes** not yet implemented. Follow `docs/laruson_automation_roadmap.md` to implement Phases 0-5 before attempting the run. The two critical changes are: (1) custom-environment source bypassing WorldClim, and (2) headless SNP-set promotion bypassing Shiny. The guide below assumes those changes are in place.

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

The dataset uses SLiM3 (Haller & Messer 2019) to simulate a population of diploid individuals distributed across a 1D environmental gradient. Multiple architectural scenarios are tested:

- **Single-locus (oligogenic)** — one or few QTNs of large effect
- **Polygenic** — many QTNs of small effect (more realistic for complex traits)

Individuals are then subjected to simulated reciprocal-transplant common garden experiments, producing fitness values across all home-deme × garden-deme combinations. This is the direct analog of a common garden experiment used to validate genetic offset predictions in empirical systems.

The key finding of the paper: GF offset correlates negatively with fitness across most scenarios, validating the method — but with important confounders (deme size effects, nonlinear environments). This benchmark tests whether ADAPTOGENE reproduces that correlation **and** whether the pipeline's multi-method GEA approach identifies the underlying QTNs.

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

### 3c. File manifest

The archive (1.14 GB) contains four simulation Case scenarios × 10 replicates per case (different random seeds). Per replicate, five file types are produced:

| File pattern | Format | Size (approx) | Role |
|--------------|--------|---------------|------|
| `Case_N_[SEED].vcf.gz` | VCF.gz | ~27.7 MB | Genotype input → `Input.vcf` (after decompress + normalize) |
| `Case_N_[SEED]_ind.txt` | TSV | ~635 KB | Individual genotype/metadata matrix (may duplicate VCF) |
| `Case_N_[SEED]_causal_mutations_pos_filtered.txt` | TSV | ~714 B | **QTN positions** (explicitly listed) → `causal_loci.tsv` |
| `Case_N_[SEED]_Freq_ML_WF.txt` | TSV | ~157 KB | Allele frequencies from WF/ML models (used for some offset methods) |
| `Case_N_[SEED]_ML_WF_CG_sum_Gen.txt` | TSV | ~2.5 KB | **Summary statistics by generation for Common Garden** — primary source for fitness values |

**Plus documentation:**

| File | Description |
|------|-------------|
| `README.txt` | Column-level data dictionary — **read this first** |
| `seeds_source_R.txt`, `seeds_Neutral.txt`, `seeds_SingleLocus.txt` | R seeds used per Case |

**Simulation cases:**

| Case | Architecture | Notes |
|------|-------------|-------|
| Neutral | No adaptive loci | Baseline; `causal_mutations_pos_filtered.txt` should be empty |
| Monogenic/SingleLocus | 1 QTN of large effect | Simplest detection test |
| Polygenic | Multiple QTNs of small effect | More realistic; harder to detect |
| Additional (N=4) | TBD from README | Confirm from README.txt |

**Recommended benchmark target.** Start with the **Polygenic case, one replicate** (e.g., the first seed in `seeds_source_R.txt`). This gives the hardest, most realistic detection test. After the harness validates, optionally loop all 10 replicates and average metrics across seeds.

> **Still unresolved (need README.txt or paper methods).** Column names in `causal_mutations_pos_filtered.txt` and `ML_WF_CG_sum_Gen.txt`, exact fitness matrix dimensions (number of demes × gardens), and the sim's chromosome/LG count and LD block size. Download `README.txt` from the Dryad page first; see `docs/laruson_dataset_notes.md` (filled during Phase 0) for confirmed values.

---

## 4. Data Model and Conversion

Láruson simulations use abstract environments — no real geography. ADAPTOGENE normally expects WorldClim bioclimatic variables extracted at real lat/lon coordinates. The conversion step bridges this gap:

- **Abstract lat/lon** — each deme is assigned a unique grid coordinate (e.g., `lat = deme_id × 2.0`, `lon = deme_id × 3.0`) so the metadata is valid for piemap rendering. These coordinates are never used for raster extraction when `Climate.source: custom`.
- **Environment variables** — the simulation's environment gradient is renamed `bio_1..bio_N` to match ADAPTOGENE's predictor naming convention.
- **VCF normalization** — chromosome prefix stripped (`chr1` → `1`), coordinates verified as integers.
- **Minimal GFF** — one synthetic gene record per QTN spans ±500 bp so gene/enrichment rules resolve without errors. No real GO terms (GFF's `go_field` is left empty in `config_LARUSON.yaml`).

### Run the converter

```bash
# Requires Gate G1 complete and benchmarks/convert_laruson.R implemented (Phase 4)
# Choose a case and seed (example: Polygenic case, first seed from seeds_source_R.txt)
CASE=3         # FILL: case number for polygenic (confirm from README.txt)
SEED=1442299973452   # FILL: first seed from seeds_source_R.txt (10 available)

docker run --user $(id -u):$(id -g) --rm --memory=20g -v $PWD:/pipeline adaptogene:latest \
  Rscript /pipeline/benchmarks/convert_laruson.R \
    data/laruson/                  `# raw archive dir (read-only)` \
    data/laruson_converted/        `# output dir` \
    ${CASE}                        `# case number` \
    ${SEED}                        `# replicate seed` \
    50                             `# TP_WINDOW_KB — fill from Gate G1 LG architecture`
```

To run all 10 replicates (after the harness validates on one):
```bash
while IFS= read -r SEED; do
  docker run ... Rscript /pipeline/benchmarks/convert_laruson.R \
    data/laruson/ data/laruson_converted_${SEED}/ ${CASE} ${SEED} 50
done < data/laruson/seeds_source_R.txt
```

### Converter outputs

| File | Description |
|------|-------------|
| `laruson.vcf` | Uncompressed, chr-normalized VCF |
| `metadata.tsv` | `site, sample, latitude, longitude, bio_1..bio_N` |
| `environments_present.tsv` | `site, bio_1..bio_N` — present environment per deme |
| `environments_future/garden_{id}.tsv` | Per garden: `site, bio_1..bio_N` — target environment |
| `causal_loci.tsv` | `chr, pos, effect_size, category` (causal / linked_neutral / background_neutral) |
| `fitness.tsv` | `site, garden, fitness` (long format) |
| `laruson_minimal.gff3` | Synthetic GFF covering QTN positions |

---

## 5. Configuration

The benchmark uses `config_LARUSON.yaml` (defined in full in `docs/laruson_automation_roadmap.md`, Phase 1a). Key settings that differ from a typical run:

```yaml
Climate:
  source:       custom                  # bypass WorldClim; use user-supplied env tables
  present_file: data/laruson_converted/environments_present.tsv
  future_dir:   data/laruson_converted/environments_future/

GEA:
  configs:
    - { method: EMMAX, adjust: bonferroni, threshold: 0.05 }
    - { method: LFMM,  adjust: bonferroni, threshold: 0.05 }
  combine_method: Sum                   # tested against Overlap and single in eval_detection.R

GWAS:
  configs: []                           # no phenotype GWAS — sim environments → GEA only

headless_snp_set: laruson_gea_combined  # triggers headless SNP-set promotion (Phase 2)

Maladaptation:
  methods:
    gradient_forest:
      run_label: laruson_gea_combined
      spatial_correction: false         # no real geography, no spatial correction
      random_model: false
    geometric_offset:
      run_label: laruson_gea_combined
```

---

## 6. Headless Run

After Phases 0-4 are complete (Gate G1, pipeline changes, converter), the full run is a single command:

```bash
bash benchmarks/run_benchmark.sh \
  --config config_LARUSON.yaml \
  --cores 4
```

The script runs pipeline modes sequentially (fail-fast), loops over gardens for the fitness matrix, promotes the GEA SNP set headlessly, and invokes the evaluation harness. No Shiny interaction is required. See `docs/laruson_automation_roadmap.md` Phase 3 for the complete spec.

### Manual step-by-step (if `run_benchmark.sh` not yet implemented)

```bash
CFG="config_LARUSON.yaml"
DOCKER="docker run --user $(id -u):$(id -g) --rm --memory=20g -v $PWD:/pipeline adaptogene:latest"
SNAKE="snakemake -c4 -s Snakefile --scheduler greedy --configfile ${CFG}"

# 1. Processing
$DOCKER $SNAKE --config mode=processing

# 2. Pre-structure
$DOCKER $SNAKE --config mode=prestructure

# 3. Structure
$DOCKER $SNAKE --config mode=structure

# 4. GEA (uses custom environment tables)
$DOCKER $SNAKE --config mode=gea

# 5. Headless SNP-set promotion
$DOCKER $SNAKE --config mode=gea headless_snp_set=laruson_gea_combined

# 6. Maladaptation (repeat per garden, varying garden_id)
for GARDEN_ID in 1 2 3; do   # FILL: actual garden IDs from fitness.tsv
  $DOCKER $SNAKE --config mode=maladaptation garden_id=${GARDEN_ID} headless_snp_set=laruson_gea_combined
done

# 7. Evaluate
$DOCKER Rscript /pipeline/benchmarks/eval_offset.R \
  LARUSON_results/Maladaptation/tables/ \
  data/laruson_converted/fitness.tsv \
  data/laruson_converted/metadata.tsv \
  gradient_forest,geometric_offset \
  laruson_gea_combined \
  LARUSON_results/benchmark_eval.tsv \
  LARUSON_results/benchmark_plots/

$DOCKER Rscript /pipeline/benchmarks/eval_detection.R \
  LARUSON_results/GEA/tables/methods/ \
  LARUSON_results/GEA/tables/selected_snps.tsv \
  data/laruson_converted/causal_loci.tsv \
  EMMAX,LFMM \
  50 \
  LARUSON_results/benchmark_eval.tsv \
  LARUSON_results/benchmark_plots/
```

---

## 7. Evaluation Outputs

Both eval scripts append to `LARUSON_results/benchmark_eval.tsv` (long format: `axis, method_or_config, metric, value`).

### Axis 1 — Offset accuracy

| Metric | Expected direction | Interpretation if wrong |
|--------|--------------------|------------------------|
| Spearman ρ (offset vs fitness) | Negative | Check garden-ID loop alignment; verify `site` join between offset and fitness tables |
| % gardens with ρ < 0 | > 50% | Both methods should track each other; if GF works but geometric doesn't (or vice versa), check the LFMM K alignment |
| Median |ρ| across gardens | > 0.3 is a reasonable pilot bar | Low but negative = correct direction, weak signal (expected for polygenic sims) |

### Axis 2 — Causal-SNP detection

| Metric | Interpretation |
|--------|----------------|
| Recall (causal) | Fraction of true QTNs recovered. Very low recall = GEA misses adaptive signal; check K value, method thresholds |
| Precision (excl. linked) | Fraction of called SNPs that are true QTNs (excluding expected-linked hits). Very low = mostly noise; check MAF filter, combine threshold |
| Expected-linked hits | SNPs near QTNs (within TP window) called as significant. These are **expected** — do not treat as FP |
| Combined > best-single | Combined (Sum or Overlap) should not be Pareto-dominated by any single method |

### Plots produced

| Plot | What to look for |
|------|------------------|
| `<method>_garden_<id>_offset_vs_fitness.png` | Negative trend + loess curve; no visible horizontal cloud |
| `precision_recall.png` | Combined config point should be up-right of all single-method points |
| `confusion_breakdown.png` | Causal TP bar grows with combination; background-neutral FP shrinks |

---

## 8. Expected Results and Interpretation

Based on Láruson et al. 2022 (Tables 2-4, Figures 3-4):

**Axis 1.** GF offset should negatively correlate with fitness across most (not all) gardens and both architectures. Effect size is typically moderate (Spearman |ρ| ~ 0.3–0.6 for oligogenic, weaker for polygenic). Geometric offset (LEA::genetic.gap) was not evaluated in the original paper — this benchmark provides new validation data for that method. A positive or near-zero ρ after controlling for the garden-ID join is a pipeline failure, not a dataset issue.

**Axis 2.** Recall for QTNs is expected to be imperfect — adaptive SNPs at intermediate frequency and in polygenic architectures are the hardest to detect. What matters is that:
(a) the pipeline produces recall above a random-draw null (recall > QTN_count / total_SNPs)
(b) combining methods shifts the precision-recall curve upward or to the right relative to any single method

**Caveats.**
- Abstract environments have no spatial correlation structure; the spatial correction (`spatial_correction: true`) is turned off in `config_LARUSON.yaml` because the spatial PCNM approach assumes geographic autocorrelation.
- Linked-neutral hits (within TP_WINDOW_KB of a QTN) are expected to inflate the called-SNP list — this reflects LD, not a pipeline bug. Report them as a separate category, not as false positives.
- The Láruson simulations were designed with relatively simple demographic models. More complex realistic demography (admixture, bottlenecks, isolation-by-distance) may produce different confounder patterns — see Láruson 2022 Section 4.2 for a discussion.

---

## 9. Caveats, Gates, and Known Limitations

| Caveat | Detail |
|--------|--------|
| Gate G1 required | The manifest table (Section 3c), TP window (Section 7), and garden IDs (Section 6) must be filled from the Phase 0 deep-research pass before any conversion or evaluation |
| Pipeline changes required | `Climate.source: custom` + headless SNP-set promotion are not yet implemented; see `docs/laruson_automation_roadmap.md` Phases 1-2 |
| No GO enrichment | The minimal GFF has no GO annotation; enrichment plots will be empty. This is expected and not a failure |
| Abstract map quality | Piemaps use abstract lat/lon grid coords; the geographic output is meaningless. Ignore map-based plots for this benchmark; focus on tabular outputs |
| Fitness matrix vs scalar | If Gate G1 confirms a scalar (one fitness per pop, not a matrix), the garden loop in `run_benchmark.sh` is simplified to a single maladaptation run; update `run_benchmark.sh` accordingly |
| MVP benchmark | Lind & Lotterhos 2025 (BCO-DMO/Zenodo) is the publication-grade follow-up but requires a separate converter and larger compute; not in scope for this pass |

---

## References

- Láruson ÁJ et al. (2022). Comparing Measures of Adaptive Genetic Variation for Use in Conservation — A New Tool for Conservationists. *Evolutionary Applications* 15(7):1048-1063. [DOI: 10.1111/eva.13354](https://doi.org/10.1111/eva.13354)
- Haller BC & Messer PW (2019). SLiM 3: Forward Genetic Simulations Beyond the Wright-Fisher Model. *Molecular Biology and Evolution* 36(3):632-637. [DOI: 10.1093/molbev/msy228](https://doi.org/10.1093/molbev/msy228)
- Fitzpatrick MC & Keller SR (2015). Ecological genomics meets community-level modelling of biodiversity: mapping the genomic landscape of current and future environmental adaptation. *Ecology Letters* 18:1–11. [DOI: 10.1111/ele.12376](https://doi.org/10.1111/ele.12376)
- Gain C et al. (2023). A quantitative theory for genomic offset statistics. *Molecular Biology and Evolution* 40(6):msad140. [DOI: 10.1093/molbev/msad140](https://doi.org/10.1093/molbev/msad140)
- Lind BM & Lotterhos KE (2025). Evaluation of genomic offset predictions using individual fitness measured in common gardens. *Molecular Ecology Resources*. [DOI: 10.1111/1755-0998.14008](https://doi.org/10.1111/1755-0998.14008) *(The gold-standard framework; this benchmark follows its evaluation design)*
- Lotterhos KE (2024). Model Validation Program dataset (BCO-DMO: 10.26008/1912/bco-dmo.889769.1; Zenodo: 10.5281/zenodo.7622893) *(Future extension)*
