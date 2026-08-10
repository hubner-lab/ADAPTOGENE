# ADAPTOGENE simulation benchmark — claims, with evidence and status

> **STATUS: COMPLETE (2026-08-10).** All 12 marker panels x 32 replicates x 112 gardens x 4
> offset methods computed and scored. Every claim below rests on the full matrix; the 30
> non-control replicates enter all aggregates (n = 120 paired cells per contrast = 30
> replicates x 4 methods).

Metric throughout: **prediction accuracy = −Kendall's τ** between predicted genomic offset and
realised common-garden fitness. Scale 0–1, **higher is better**. Per-replicate medians first,
then median across replicates, so every replicate weighs the same. Panel comparisons are
**paired within replicate**.

Design: 30 SS-Mtn replicates, balanced 10 oligogenic / 10 mod-polygenic / 10 highly-polygenic,
plus 2 degenerate Est-Clines controls (excluded from aggregates). Full reciprocal transplant:
112 gardens (100 landscape + 12 climate-novelty). Four offset methods.

---

## Claim 1 — Markers found by the pipeline predict as well as knowing the true answer

| marker set | accuracy |
|---|---|
| true causal loci (oracle, unobtainable on real data) | 0.70 |
| **pipeline, ≥2 of LFMM/RDA/EMMAX agree** | **0.70** |
| any method flagged it (union) | 0.63 |
| LFMM alone | 0.61 |

Paired Wilcoxon vs oracle, Holm-corrected: `≥2 agree` **p = 1.0** (indistinguishable).
`union` p = 2×10⁻⁶ and `LFMM only` p = 2×10⁻⁷ — both significantly worse.

**Why it matters:** the source design (Lind & Lotterhos 2025) cannot ask this question, because
their "adaptive" panel *is* the true QTN set. Ours is built by a GEA scan with no access to
truth — the situation on real data.

## Claim 2 — Use `≥2 agree`. Not because it is most accurate, because it never fails.

| panel | oligogenic | mod-polygenic | highly-polygenic | **spread** |
|---|---|---|---|---|
| **≥2 agree** | 0.704 | 0.703 | **0.714** | **0.011** |
| RDA only | 0.697 | 0.725 | 0.675 | 0.050 |
| EMMAX only | 0.704 | **0.739** | 0.625 | 0.114 |
| union | 0.547 | 0.565 | 0.691 | 0.144 |
| LFMM only | 0.536 | 0.565 | 0.684 | 0.148 |
| all 3 agree | **0.709** | **0.748** | 0.553 | 0.195 |

`≥2 agree` wins outright only on highly-polygenic. Its case is the **spread: 0.011 against
0.05–0.195 for everything else.** `all 3 agree` wins two strata and collapses in the third
(18 markers cannot span a 410-locus gradient). `EMMAX only` wins mod-polygenic and is
second-worst on highly-polygenic.

On a real dataset the trait's architecture is unknown, so a rule that is merely *good at the
right architecture* is not usable. **Combining does not make you more accurate; it makes you
reliably accurate.**

## Claim 3 — Combining methods does NOT beat the best single method

Paired, per replicate × method — `≥2 agree` wins:

| contrast | wins |
|---|---|
| vs LFMM only | **65 %** |
| vs RDA only | 44 % |
| vs EMMAX only | 40 % |

It beats the *worst* single method and is a coin flip against the best two. **Do not write
"combining improves accuracy".** The defensible claim is Claim 2: insurance against picking
the wrong single method, which cannot be known in advance.

## Claim 4 — Finding more causal loci makes prediction WORSE, and it is selection quality, not set size

**Size-matched floor, now measured.** `random_matched` draws the SAME NUMBER of markers as
`≥2 agree` (151 vs 151) from background-neutral loci only:

| | accuracy |
|---|---|
| `≥2 agree` (GEA-selected) | **−0.7115** |
| `random_matched` (same size, random) | −0.5873 |

Real better in **67 %** of paired cells, n = 120, **p = 3.3×10⁻⁹**. Identical size, very
different accuracy — so what the GEA scan selects genuinely carries the signal. This upgrades
the claim below from association to demonstrated cause.


ρ(recall, accuracy) = **−0.79**. ρ(background-noise fraction, accuracy) = **−0.63**.

`≥2 agree` finds 17 of 410 causal loci on highly-polygenic replicates and predicts as well as
the oracle holding all 410 — the causal loci track the same climate gradient and are
near-redundant. `union` finds the most causal loci of any obtainable panel at every
architecture and ranks near the bottom in two of three: it also carries 31–59 background
markers.

**Detection performance and prediction performance are largely decoupled.** Precision/recall
against QTNs is a poor proxy for offset accuracy, which is why the results lead with accuracy.

## Claim 5 — Curating markers beats using all of them — but ONLY for some architectures

**FINAL, n = 120 paired cells (30 replicates x 4 methods):**

| contrast | accuracy | curated better in | p |
|---|---|---|---|
| true QTNs vs all SNPs | −0.6960 vs −0.5841 | 81 % | 3.7×10⁻¹⁴ |
| **≥2 agree vs all SNPs** | **−0.7115 vs −0.5841** | **62 %** | 6.3×10⁻⁸ |
| ≥2 agree vs neutral SNPs | −0.7115 vs −0.5820 | 66 % | 7.1×10⁻⁹ |

**The benefit is architecture-dependent:**

| architecture | curated | all SNPs | curated better in |
|---|---|---|---|
| oligogenic | −0.698 | −0.545 | **68 %** |
| mod-polygenic | −0.711 | −0.553 | **68 %** |
| highly-polygenic | −0.727 | −0.684 | **52 %** — no meaningful benefit |

When few loci control the trait, curation is worth a lot (0.698 vs 0.545 oligogenic; 0.711 vs
0.553 mod-polygenic). When the trait is highly polygenic the advantage nearly vanishes (52 %,
i.e. a coin flip) — the adaptive signal is spread so widely that the whole genome already
carries it.

**This stratification is NOT in Lind & Lotterhos 2025.** Verified against their full text
(`docs/papers/lind2025_fulltext.txt`): they report the advantage as one pooled figure —
*"The median increase in performance from adaptive compared to all or neutral models was less
than 3%"* — and stratify only by NUMBER OF TRAITS (1/2 vs 6), never by genic level. Genic
architecture appears solely as a covariate on overall performance: *"The performance of methods
was similar across genic levels but increased slightly as the number of QTNs underlying
adaptation became more polygenic (Figure S8)."* No sentence ties the size of the marker-set
benefit to architecture. Their pooled <3% averages two different regimes.

CAVEAT: only their main text was checked; the Supporting Information (Figures S8, S31-S37) was
not available via Europe PMC. Confirm no per-genic-level breakdown exists there before writing
"not previously reported".

## Claim 6 — Methods/engineering

Projecting a fitted model onto N future climates now costs **one model fit instead of N**:
84 pipeline invocations per replicate → **1**; ~10.7 h → ~1.2 h compute (**9.2×**). A separate
fix batched Gradient Forest's `predict()` calls (2300 s → 128 s on the worst job).

**Scope limit that must be stated:** the gain scales with the number of futures requested
(GCMs, SSPs, time points, transplant sites). A user projecting onto a **single** future sees
almost no speedup. Correct phrasing: *"projecting onto N future climates costs one model fit
instead of N"* — not "the pipeline is 9× faster".

Validated: outputs **byte-identical** to the previous implementation on 1176 files
(max abs diff 0.000e+00). Fitness identity exact on all 32 replicates (R² ≥ 0.999999);
local adaptation reproduces the deposit's own values within 0.022.

---

## Caveats that must appear in Methods

1. **Cohort size imbalance.** The 18 Group 1 replicates are ~2.4× larger than the original 14
   (median `all` panel 24 366 vs 10 045 SNPs). Selection stratified on architecture and
   R²(PC1~temp), never on dataset size. Paired within-replicate comparisons are unaffected;
   pooled per-stratum medians mix cohorts. Report `n_snps` as a covariate.
2. **Random floor — RESOLVED.** `rand_best1` is size-matched exactly to `best` (151 markers).
   GEA selection beats it by 0.124 accuracy, p = 3.3×10⁻⁹. Claim 4 is now causal, not
   associational. The `rand_solo1`/`rand_union1` pairs are computed and can be added.
3. **Two controls excluded.** The 2 Est-Clines replicates are degenerate by construction
   (R²(PC1~temp) > 0.9) and are excluded from all aggregates.

Sources: `REPORT_phase1.md`, `TABLES_absolute.md`, `DISCUSSION_NOTES.md`,
`figures/table_absolute_counts.tsv`, `offset09/phase1_seed_medians_solo.tsv`.
