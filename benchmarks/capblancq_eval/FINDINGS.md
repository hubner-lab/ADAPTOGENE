# Capblancq deposit — screen results (2026-08-04)

Data: Dryad doi:10.5061/dryad.jdfn2z3m6, CC0, downloaded 2026-08-04.
`data/capblancq/raw/` (archives + `SHA256SUMS.txt` + `PROVENANCE.txt`), unpacked to
`data/capblancq/unpacked/`. Screen script: `benchmarks/capblancq_screen.py`.
Per-run table: `benchmarks/capblancq_eval/screen_all.tsv` (400 rows, no failures).

## Deposit shape

**400 runs, not 100**: 4 selection-gradient shapes (`alpha0.5`, `alpha1`, `alpha2`, `alpha3`)
× 100 replicates. α is the exponent in the fitness function and is the paper's central variable:

```
fitness = exp(-(1/2) * ( (|opt_var1 - pheno_var1| / sigma_K)^α
                       + (|opt_var2 - pheno_var2| / sigma_K)^α ))
```

α=2 is Gaussian, α=1 Laplace/exponential — matching the paper's theory that geometric offset
predicts fitness under Gaussian selection while RONA implicitly assumes a Laplace gradient.

Per run: `genome_step1.vcf` (~98 MB), `fitness_step1`, `fitness_pred`, `mutationm2_step1.txt`,
`mutationm3_step1.txt`, `position_ind_step1`, `trait1`, `trait2`, `var1_step1`, `var1_pred`,
`var2_step1`, `var2_pred`. Unpacked total 32 GB.

Genome: **single chromosome, ~110 kb**, `initializeMutationRate(3e-6)`,
`initializeRecombinationRate(1e-2)`. Spatial: 12×12 continuous landscape, nonWF, `K=11`.
`trait1 = sumOfMutationsOfType(m2)`, `trait2 = sumOfMutationsOfType(m3)` — **m2 is causal for
var1 only, m3 for var2 only**, giving clean per-axis truth.

## The five checks — all resolved

| # | Question | Answer |
|---|---|---|
| 1 | 120 or 240 causal? | **~240 total, ~120 per axis.** Range 227–240 (median 239); m2 median 119, m3 median 120. Varies per run — must be counted, never assumed |
| 2 | Fitness per-individual or per-deme? | **Per-individual.** 1,108–1,788 individuals per run (median 1,639). `position_ind_step1` ships x,y so we can bin demes ourselves. Env values are also per-individual |
| 3 | Causal positions in VCF coordinates? | **`VCF_POS = mutfile_pos + 1`**, uniform across all 119/120 (standard SLiM 0-based → VCF 1-based). Better: causal SNPs are labelled directly in the VCF as `MT=2` / `MT=3`, so truth can be built from the VCF alone |
| 4 | Unlinked-neutral set? | **No QTL-free linkage groups** — one chromosome. But `r = 1e-2/bp` over 110 kb ≈ 1,100 Morgans, so **every locus is effectively unlinked**. There is no `linked_neutral` class; all `MT=1` SNPs are `background_neutral`, and `precision_strict` is computable cleanly. Cost: no LD-based detection, region clustering and WZA are meaningless here |
| 5 | Re-adaptation before post-shift fitness? | **None, by construction.** `fitness_step1` and `fitness_pred` are computed in the *same* `12000 late()` block from the *same* phenotypes (`pheno1_var1`, `pheno1_var2`); only the optimum map changes (`map1`/`map3` → `map2`/`map4`). Zero generations elapse |

**MAF trap: absent** — the opposite of the retired Láruson dataset. Recall ceiling at MAF 0.01 has
median 1.000 (337/400 runs lose zero causal loci); at MAF 0.05 median 0.983, worst 0.923. Causal
MAF median 0.367.

## Structure — the blocking finding for the GEA arm

Measured on **neutral SNPs only** (`MT==1`, MAF ≥ 0.05), which is the honest diagnostic here: causal
loci are ~17% of the MAF-0.05 set, so an all-SNP PCA is dominated by the adaptive signal itself.
Both PCAs agree, so this is not an artefact.

| α | `pc1_var_neutral` median | `r2_pc1_env_neutral` median (min–max) |
|---|---|---|
| 0.5 | 0.144 | **0.747** (0.549–0.915) |
| 1 | 0.138 | **0.856** (0.613–0.956) |
| 2 | 0.140 | **0.906** (0.652–0.961) |
| 3 | 0.145 | **0.885** (0.703–0.950) |

- **In the target band (0.3–0.6): 1 / 400.**
- **Degenerate (≥ 0.8): 279 / 400.**
- `r2_env1_env2` median 0.0006 — the two gradients are orthogonal exactly as designed.

This is the **MVP `Est-Clines` failure mode**, not the Láruson one: structure magnitude is genuine
(PC1 ≈ 14%, versus Láruson's 0.38%), but it is *aligned with* the environment. Cause is structural
to the design — a continuous landscape with local dispersal plus an environment that is a smooth
function of x and y makes isolation-by-distance a near-perfect proxy for the environment.

**Consequence.** Structure correction can only remove signal on this dataset, so it cannot be used
to find the optimal LFMM `K` / EMMAX `n_pcs` / RDA `condition_pcs` for real data. The collapse
curve is still measurable and reportable, and α provides a partial lever (α=0.5 is least
degenerate, and holds the only in-band run).

## Offset arm — usable, with a demanding null

`r2_null_dist_vs_logfit` = R² of Euclidean climate-transfer distance against realized log-fitness
decline. Any genomic-offset method must beat this to add value.

| α | median null R² | median log-fitness decline |
|---|---|---|
| 0.5 | **0.327** | −0.538 |
| 1 | 0.769 | −0.668 |
| 2 | 0.750 | −0.749 |
| 3 | 0.513 | −0.983 |

Overall min 0.094, median 0.635, max 0.870. For reference, Gain et al. 2023 reported ≈0.45 for
their equivalent null.

**α is therefore a difficulty knob for the offset arm**: α=0.5 leaves the most headroom (null 0.33),
α=1 and α=2 leave very little (null ≈0.75–0.77). Reporting offset performance without this null,
per α, would be misleading.

## Track B extensibility (SLiM scripts)

4 scripts, 321 lines each, **identical except the α exponent and a hardcoded `setwd()`**. The
2-axis design is hardcoded, not parameterized: ~42 axis-specific references, two mutation types
(`m2`, `m3`), four spatial maps (`map1`–`map4` = var1 present/future, var2 present/future) built by
explicit `mapVar1`/`mapVar2` blocks, and the fitness expression written out term-by-term in three
places (lines 249, 309, 315).

Extending to ≥5 axes is mechanical but touches ~10 sites: add mutation types `m4`…, add map-pair
construction blocks, extend the phenotype and fitness expressions. There is no loop over axes to
generalize. Note also that the maps are built as independent per-column / per-row uniforms, so
**adding correlated axes requires replacing the map construction**, not just repeating it — which
is the substantive part of the work, per the plan's §0.2/§0.3.

---

# Conclusions and decision (2026-08-04)

## Capblancq vs MVP, measured

The screen was run to decide whether Capblancq should replace or supplement MVP. Compared against
our 14 selected MVP seeds, using the MVP authors' own shipped statistics
(`data/mvp/selection/summary_20220428_20220726.csv`, `cor_PC1_temp` squared to R²):

| | **MVP (our 14 seeds)** | **Capblancq (400 runs)** |
|---|---|---|
| environmental predictors | 2 | 2 — **no improvement** |
| R²(PC1 ~ env) | median **0.442** (0.310–0.940) — **in the 0.3–0.6 band** | median 0.75–0.91; **1/400 in band, 279/400 degenerate** |
| orthogonal control axis | **yes** — salinity at R² = 0.006 | none; both axes equally degenerate |
| mean F<sub>ST</sub> | median 0.192 — in the 0.05–0.20 target range | — |
| LD | usable — `linked_neutral` class exists, region/window methods meaningful | **none** — `r = 1e-2` leaves every locus unlinked |
| architectures | oligogenic (2–15 causal) → highly polygenic (268–395); 1-trait and 2-trait arms | one architecture, ~240 causal, always 2 traits |
| replicates | 2,250 available, 14 converted | 400 |
| prior work | converted, scored, journals 04–07, published baseline to compare against | none |

**MVP is better on every axis that matters for GEA parameter tuning.** The seed-selection exercise
did its job: MVP sits in the informative band that Capblancq misses in 399 of 400 runs.

Capblancq's only genuine advantages:

1. **Post-shift fitness and future environment ship directly** — no need to reimplement the fitness
   function, so no risk of getting it wrong.
2. **α — four selection-gradient shapes.** MVP has no equivalent. It tests the paper's own theory
   (geometric offset ↔ Gaussian, RONA ↔ Laplace) and empirically controls offset difficulty: the
   climate-distance null moves 0.33 → 0.77 across α.

It does **not** fix the multivariate problem, which is what motivated the search in the first place.

## MVP does support a genetic-offset arm

We never ran one, but nothing blocks it. Our converted MVP already ships per-individual phenotypes
and per-deme optima:

```
site  sample  latitude  longitude  phen_temp  phen_sal  fitness     # per individual
site  bio_1  bio_2                                                  # per deme = the optima
```

Transplant fitness of individual *i* in deme *j* is a published function of
(phenotype_i − optimum_j) — Lotterhos 2023 Eq. 3 — which is exactly what the MVP-offsets pipeline
computes to build its 100 reciprocal-transplant and 11 climate-novelty gardens.

## Correction to an earlier claim

The cross-arm question — *does feeding GEA-selected SNPs into the offset model beat using all
SNPs?* — was previously recorded as answerable only on a deposit carrying both truth and fitness.
That is backwards: it requires structure in the informative band, so it is answerable on **MVP**
and **not** on Capblancq.

## Decision

- **MVP is the primary dataset for both arms.** GEA parameter tuning (structure in band, orthogonal
  control axis available) and genetic offset (transplant fitness computed from data already held).
- **Capblancq is the validation set.** Because its post-shift fitness is shipped rather than
  derived, it validates our transplant-fitness generator and `benchmarks/eval_offset.R` before
  either is trusted on MVP. The α result stands on its own as a reportable figure MVP cannot
  produce.
- **Track B (extending the simulation to ≥5 predictors) is parked.** It does not fix the
  multivariate gap on either dataset, and the map construction — not the axis count — is the real
  work.

## Standing field-level result

Three independent literature searches (2026-08-04) found **no published deposit with ≥5
environmental predictors + a causal truth set + realized fitness**. The gap is structural:
truth+fitness deposits are all 2-causal-axis; ≥5-predictor deposits are empirical and so carry no
truth set. Anchor for the paper: **Forester et al. 2018 — the benchmark that established RDA for
multivariate GEA — uses a single causal environmental axis** (verified from the author-hosted PDF
with SI, Methods §2.3–2.6: *"For simulation data, only one predictor variable is used"*). Full
survey in the session plan file.

## Next steps (not started)

1. `benchmarks/eval_offset.R` — score predicted offset against realized log-fitness (Spearman ρ and
   R²) against the Euclidean climate-distance null. Does not exist; `benchmarks/` has
   `eval_detection.R` for p-values-vs-truth only.
2. MVP transplant-fitness generator — Lotterhos 2023 Eq. 3 over our converted phenotypes + optima.
3. Validate both on Capblancq (fitness shipped), then run on MVP.
4. Offset parameter sweeps: `gradient_forest` (`ntree`, `cor_threshold`, `spatial_correction`),
   `geometric_offset` (`scale`, `k`), `rda_offset` (`axes`, `axis_alpha`, `condition_pcs`), and
   `Maladaptation.snp_sets` = all vs GEA-selected — the cross-arm question.
5. Retarget `benchmarks/score_ladders.sh` and `sweep_rda.sh`, which still default to the retired
   `LARUSON1K` project.
