# GEA hyperparameter optimization on the Láruson benchmark — findings

Dataset: Láruson et al. 2022, Case 1 (cline), replicate seed 2889863491989.
Truth: `data/laruson_converted/causal_loci.tsv` — 102 causal / 27,186 linked_neutral /
5,134 background_neutral positions.

Two projects:
- **LARUSON** — full cohort, 10,000 samples (100 sites x 100), 5,583 SNPs post-filter, 28 causal testable.
- **LARUSON1K** — 1,000-sample subsample (10 per site, seed 3), 5,601 SNPs post-filter,
  30 causal testable on the full set / 24 on the LD-pruned set used by `mode=pregea`.
  Seed 3 was selected because it retains **all 28** of the 10k-passing causal loci, so the
  two stages score against the same truth (`data/laruson_1k/subsample_provenance.tsv`).

Scoring: `benchmarks/eval_detection.R` (+ `lib_detection.R`), which reuses the pipeline's own
`compute_pval_threshold()` so threshold calls reproduce `find_sig_snps.R` exactly. Verified
against the existing 10k outputs: EMMAX 0 hits, LFMM bonf 10/10, LFMM qval 15/16, RDA qval
15/15, ceiling 28 — all reproduced exactly.

---

## 1. The dominant constraint: MAF, not any GEA parameter

Only **28 of 102** causal loci survive `Filter.maf: 0.05`. The causal MAF distribution at the
full cohort is median **0.026**, q75 0.060 — the 74 lost loci sit far below the filter, not at
its margin. No GEA hyperparameter can recover them.

Recall is therefore reported against two denominators throughout: `recall_testable` (/28 at 10k,
/30 at 1k full set, /24 on the LD-pruned set pregea uses) and `recall_simulated` (/102). MAF was
held fixed per instruction; a MAF arm remains available and unstarted.

**Because the testable denominators differ between stages, 1k and 10k results are compared in raw
TP counts, never in `recall_testable`.** The best 10k configuration reaches `recall_simulated`
= 0.167 (17/102) — the ceiling costs more than every parameter tuned in this document combined.

## 2. EMMAX — root cause of zero hits found

`EMMAX.params.n_pcs` defaults to `sNMF.k_best` (= 3 here). On this cline landscape the PCA
covariates are collinear with the environment and absorb the signal.

| n_pcs | 0 | 1 | 2 | 3 (default) | >=4 |
|---|---|---|---|---|---|
| lambda_gc | 0.990 | 1.007 | 0.981 | 0.979 | 0.973-0.979 |
| hits_qval | 5 | 8 | 2 | **0** | 0 |
| AUC-PR (vs truth) | **0.580** | 0.383 | 0.062 | **0.052** | 0.049-0.055 |
| best F1 | **0.619** | 0.471 | 0.176 | 0.150 | 0.120-0.143 |

lambda is ~1.0 at *every* rung — calibration diagnostics alone would never flag this. Only the
truth-based score exposes it: **AUC-PR collapses 11x from n_pcs=0 to the default n_pcs=3.**
This is why the 10k baseline produced zero EMMAX calls at both bonf 0.05 and qval 0.1.

pregea's own selection rule (smallest n_pcs with lambda in [0.90, 1.15]) independently picks
n_pcs=0 — heuristic and truth agree here.

**Recommendation: `EMMAX.params.n_pcs: 0`.**

## 3. LFMM — K is not a lever on this dataset

Sweeping K = 1..10 (the pipeline default +/-2 would have been far too narrow):

| | K=1 | K=2 | K=3 (k_best) | K=5 | K=9 | K=10 |
|---|---|---|---|---|---|---|
| lambda_gc (raw) | 3.824 | 3.814 | 3.787 | 3.788 | 3.772 | 3.787 |
| hist_flatness_ks | 0.3137 | 0.3151 | 0.3151 | 0.3135 | **0.3124** | 0.3129 |
| AUC-PR (vs truth) | **0.770** | 0.767 | 0.765 | 0.764 | 0.764 | 0.763 |
| best F1 | **0.791** | 0.780 | 0.780 | 0.780 | 0.769 | 0.769 |

**lambda is flat across the entire ladder** (3.77-3.82) — adding latent factors does not reduce
inflation at all. The inflation is genuine genome-wide association (polygenic architecture plus
IBD on a cline), not unmodelled structure: the environment *is* the dominant axis of genetic
variation, so no number of latent factors can separate confounder from signal.

Truth-based scoring agrees the knob is inert: AUC-PR spans 0.763-0.770 across K=1..10.

Note pregea's rule (min `hist_flatness_ks` among rungs retaining >=50% of max hits) selects
K=9, which is marginally the *worst* rung on truth. The gap is negligible (0.764 vs 0.770), so
this is a tie-breaking artefact rather than a bad heuristic — but it does show that flatness
alone does not track detection power here.

**Recommendation: leave LFMM K at `sNMF.k_best`; spend the effort on the threshold instead.**
(For reference: at the 10k cohort LFMM reported gif = 29.9, consistent with lambda scaling with
n for a real polygenic signal rather than indicating miscalibration.)

## 4. RDA — `condition_pcs: 0`, confirmed across the full 0..10 ladder

Per-SNP p-values per rung from `benchmarks/sweep_rda.sh` (full marker set, n=1000, 5,601 SNPs,
30 testable causal). gif is stable at 1.34-1.37 across the ladder, so this is not a calibration
effect:

| condition_pcs | 0 | 1 | 2 | 3 | 4 | 5 | 6-10 |
|---|---|---|---|---|---|---|---|
| AUC-PR | **0.699** | 0.573 | 0.187 | 0.218 | 0.218 | 0.217 | 0.190-0.211 |
| best F1 | **0.714** | 0.593 | 0.258 | 0.340 | 0.326 | 0.326 | 0.279-0.318 |
| TP at best F1 | **20** | 16 | 8 | 8 | 7 | 7 | 6-7 |

Conditioning on structure removes the signal, and the drop is compounded by `rda.R`'s B6 double
fit, which takes `pmax()` of the partial and unconstrained p-values whenever `condition_pcs > 0`
and is therefore strictly conservative. The existing `config_LARUSON.yaml` already sets 0; this
confirms it empirically rather than by assumption.

## 5. Thresholds are the real lever (free, no refits)

Sweeping adjust x threshold on the **existing 10k p-value tables** — no model refitting —
already beats the shipped configuration:

| config | called | TP | background FP | precision | recall (/28) |
|---|---|---|---|---|---|
| RDA qval 0.1 *(current default)* | 15 | 15 | 0 | 1.00 | 0.536 |
| RDA qval 0.2 | 18 | 16 | 0 | 1.00 | 0.571 |
| RDA custom 1e-3 | 22 | 17 | 1 | 0.94 | 0.607 |
| LFMM top 20 | 34 | 18 | 5 | 0.78 | 0.643 |
| union of 3 methods @ qval 0.1 | 17 | 16 | 0 | 1.00 | 0.571 |
| RDA top 100 | 99 | 20 | 9 | 0.69 | 0.714 |

Pareto frontier (10k): `benchmarks/laruson_eval/LARUSON10k_baseline_pareto.tsv`,
plot `LARUSON10k_baseline_precision_recall.png`.

## 5b. End-to-end confirmation at 1k, where combining does beat every single method
<!-- Scope note: the combine-wins result below holds at 1k. At 10k it only ties — see section 6. -->


`mode=gea` re-run on LARUSON1K with the tuned parameters (EMMAX n_pcs=0, LFMM K=k_best,
RDA condition_pcs=0). At the shipped `qval 0.1` threshold, EMMAX went from **0 hits to 5** —
the fix works end to end, not just in the ladder.

Full threshold x combine sweep (30 testable causal, exact-position matching), best recall at
precision >= 0.9:

| config | called | TP | background FP | precision | recall (/30) |
|---|---|---|---|---|---|
| **≥2-of-3 methods, top 20** | 22 | **17** | **0** | **1.000** | **0.567** |
| LFMM top 20 | 35 | 17 | 1 | 0.944 | 0.567 |
| union top 10 | 26 | 17 | 1 | 0.944 | 0.567 |
| RDA top 20 | 19 | 16 | 0 | 1.000 | 0.533 |
| rank-sum top 30 | 30 | 14 | 1 | 0.933 | 0.467 |
| EMMAX top 10 | 17 | 10 | 1 | 0.909 | 0.333 |

**The `at_least_2` combine rule Pareto-dominates every single method**: it matches the best
single-method recall (17 TP, tied with LFMM) while being the only configuration at that recall
with zero background false positives, and it does so calling 22 SNPs rather than LFMM's 35.
This satisfies the acceptance criterion in `docs/laruson_automation_roadmap.md` Phase 5 — the
combine step adds real value on this architecture rather than being dominated by a single method.

Note EMMAX only becomes a useful combine member *after* the n_pcs fix; at the default it
contributes nothing, which would have made any ≥2-of-3 rule collapse to an LFMM∩RDA intersection.

Pareto frontier: `LARUSON1K_tuned_pareto.tsv`, plot `LARUSON1K_tuned_precision_recall.png`.

## 6. Transfer to the full 10,000-sample cohort

`mode=gea` re-run on LARUSON with the tuned parameters. Only EMMAX and RDA re-ran; LFMM was
correctly reused from cache, since its parameters did not change.

**EMMAX, before vs after the `n_pcs` fix, at every FDR level (28 testable causal):**

| qval | 0.05 | 0.1 | 0.2 | 0.3 |
|---|---|---|---|---|
| TP before (`n_pcs`=3, the default) | 0 | 0 | 0 | 0 |
| TP after (`n_pcs`=0) | 4 | **5** | **8** | 8 |
| background FP after | 0 | 0 | 0 | 0 |

EMMAX goes from contributing nothing at any threshold to 8 true positives with zero background
false positives. The structural choice transferred cleanly from 1k to 10k, as expected — it
tracks the landscape, not the sample size.

**Best configurations at 10k (exact-position matching, 28 testable causal):**

| config | called | TP | background FP | precision | recall (/28) | **recall (/102)** | F1 |
|---|---|---|---|---|---|---|---|
| RDA custom 1e-3 | 22 | **17** | 1 | 0.944 | 0.607 | **0.167** | **0.739** |
| RDA qval 0.2 | 18 | 16 | **0** | **1.000** | 0.571 | 0.157 | 0.727 |
| union of 3 @ qval 0.1 | 17 | 16 | **0** | **1.000** | 0.571 | 0.157 | 0.727 |
| ≥2-of-3 @ top 20 | 20 | 16 | 1 | 0.941 | 0.571 | 0.157 | 0.711 |
| *shipped default (RDA qval 0.1)* | 15 | 15 | 0 | 1.000 | 0.536 | 0.147 | 0.698 |

**The `/102` column is the point of section 1.** The best configuration recovers 17 of 102
simulated causal loci — 16.7%. Every parameter tuned in this document moves the `/28` column;
the gap between 0.607 and 0.167 is the MAF filter, and it is larger than the total effect of all
the tuning combined.

**Why the union row did not move when EMMAX was fixed.** EMMAX went from 0 to 5 calls at
qval 0.1, yet `union @ qval 0.1` reports the same 17 called / 16 TP as the baseline. Verified,
not assumed: EMMAX's 5 tuned calls are a **strict subset** of LFMM ∪ RDA (`EMMAX ⊆ LFMM ∪ RDA`
is TRUE, zero positions outside). All 5 are causal — EMMAX at n_pcs=0 is precise but recovers
only the loci the other two already find. So the union is genuinely unchanged; the fixed EMMAX
adds confirmation, not new discoveries, at this threshold.

**Honest caveat on the combine claim.** At 1k, `≥2-of-3` strictly Pareto-dominated every single
method (the only zero-FP configuration at maximum recall). At 10k it does **not**: RDA alone at
`qval 0.2` matches the best combine rule exactly (16 TP, 0 FP). With 10x the samples RDA is
strong enough on its own that combining adds nothing on this architecture.

Compare in **TP counts, not recall** — the two stages have different denominators (30 testable
at 1k, 28 at 10k), so the near-identical 0.567 vs 0.571 recall figures are an artefact of that,
not evidence of similarity. The raw counts are 17 TP (1k) vs 16 TP (10k) for `≥2-of-3`.

The roadmap's Phase-5 acceptance criterion ("combined strategies should not be dominated by any
single method") is satisfied at both sizes — combining is never worse — but the stronger claim
that combining *wins* holds only at the smaller sample size.

**Window sensitivity.** All tables above use exact-position matching. Re-scored at the measured
LD-decay scale (±0.6 kb, r²=0.2 at 591 bp) and at the ±5 kb window the `causal_loci.tsv`
categories were built with:

| config | window | called | TP (testable) | causal matched (any) | background FP | precision |
|---|---|---|---|---|---|---|
| RDA qval 0.2 | 0 | 18 | 16 | 16 | 0 | 1.000 |
| RDA qval 0.2 | ±0.6 kb | 18 | 16 | 17 | 0 | 1.000 |
| RDA qval 0.2 | ±5 kb | 18 | 16 | 18 | 0 | 1.000 |
| RDA custom 1e-3 | 0 | 22 | 17 | 17 | 1 | 0.944 |
| RDA custom 1e-3 | ±5 kb | 22 | 17 | 18 | 1 | 0.955 |

The conclusion is robust to window choice: precision never falls, and the hits labelled
`linked_neutral` at window 0 resolve to real QTNs detected through linked markers as the window
widens (RDA qval 0.2's 2 linked hits become 0 at ±5 kb, with `causal matched` rising 16→18).
Widening only ever helps here, so nothing in the recommendation depends on the window rule.

## 7. Recommended configuration

Written into `config_LARUSON.yaml` (and `config_LARUSON1K.yaml`), with the measured numbers in
comments at each key:

```yaml
GEA:
  configs:
    - method: "EMMAX"
      adjust: "qval"
      threshold: '0.1'
      params: {n_pcs: 0}            # <- the one change that matters
    - method: "LFMM"
      adjust: "qval"
      threshold: '0.1'              # K stays at sNMF.k_best (inert knob)
    - method: "RDA"
      adjust: "qval"
      threshold: '0.1'
      params: {condition_pcs: 0, axes: 2, permutations: 99}
```

Threshold choice depends on which point of the curve is wanted, and all three are free to switch
(no refit):
- **Zero false positives:** `RDA qval 0.2` or `union @ qval 0.1` — 16 TP, precision 1.000.
- **Maximum recall at >=0.9 precision:** `RDA custom 1e-3` — 17 TP, 1 background FP.
- The shipped `RDA qval 0.1` gives 15 TP; relaxing to `qval 0.2` costs nothing and gains one.

## 8. Search design

Model fitting and significance calling were decoupled. Each method has exactly one
structure-correction knob (LFMM K, EMMAX n_pcs, RDA condition_pcs), so the ladders are a **sum**
(53 model fits: 20 LFMM + 22 EMMAX + 11 RDA) rather than a product with the 15-point threshold
grid and 4 combine rules (~3,200 fits). Thresholds and combine rules are pure post-processing on
p-value tables that already exist.

---

## What was NOT produced

- **`PreGEA/tables/pregea_recommendations.tsv` does not exist for either project.** `mode=pregea`
  completed its LFMM and EMMAX ladders (both scored above), but its RDA block was stopped: with
  the 999 permutations originally configured, `anova.cca` had not finished rung 0 after ~35 min,
  extrapolating to many hours for 11 rungs. Since `pregea_recommend` takes the RDA ladder as
  input, the recommendations table was never written. This is a deliberate stop, not a failed
  run — `benchmarks/sweep_rda.sh` produced strictly more information (per-SNP p-values per rung,
  which pregea does not persist at all, and which is what made the truth-based RDA ladder in
  section 4 possible).
  **Consequence:** the Shiny GEA tab's method editor reads that file for parameter pre-fill and
  will find nothing. To generate it, re-run `mode=pregea` now that the config carries
  `permutations: 99` — the LFMM/EMMAX rungs are cached, so only the RDA block re-runs. Cap BLAS
  threads when doing so (see below).
- **The MAF arm (`Filter.maf: 0.01`) was not run**, per instruction to leave dataset filters
  fixed. Section 1 quantifies what it would be worth.

## Operational notes (not results)

- **BLAS thread oversubscription, not algorithmic cost, was the RDA bottleneck.** Unbounded, a
  single `rda()` rung consumed 3.5 CPU-hours in 13 minutes wall-clock and did not finish; with
  `OMP_NUM_THREADS=4` / `OPENBLAS_NUM_THREADS=4` the same rung takes **70 seconds**. This box is
  shared (15-minute load average was ~181 before any of this work started). Any future RDA run
  here should cap BLAS threads.
- `axes=<fixed int>` does **not** skip `anova.cca` — `scripts/rda.R` lines 318-329 run the full,
  by-axis and by-margin anovas unconditionally, and the latter two are unparallelized. Permutation
  count therefore drives RDA wall-clock regardless of the `axes` setting. It does not affect
  per-SNP p-values, which come from `rdadapt`'s robust Mahalanobis distance.
- A **third pyslim VCF-writer quirk** beyond the two documented in
  `docs/laruson_dataset_notes.md` section 3: three rows (1:84913, 1:109742, 1:167837) carry a GT
  allele index of 2 while declaring a single ALT, which aborts `bcftools` AC/AN recalculation.
  None are causal and none survive the MAF filter. `benchmarks/subsample_laruson.R` uses
  `bcftools view -I` (no INFO update) to pass records through verbatim.
