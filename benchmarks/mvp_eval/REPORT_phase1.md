# Group 1 interim report — GEA-detected marker panels vs the QTN oracle

**Date:** 2026-08-08 · **Panel:** 32 replicates (30 SS-Mtn primary + 2 Est-Clines controls)
**Scope of this report:** phase-1 gardens (30 landscape + 12 novelty) × the *adaptive* and
*GEA-detected* marker panels only. The reference/floor panels and the remaining 70 gardens
are still running — see **Not yet answered** at the bottom. Nothing here depends on them.

---

## 1. What was measured

Kendall's τ between predicted genomic offset and realised common-garden fitness.
**More negative is better**: a good offset predicts low fitness. Sign is left raw so the
numbers match the tables and the source paper.

Four marker panels, all derived without touching the truth set except `adaptive`:

| panel | what it is | obtainable on real data? |
|---|---|---|
| `adaptive` | the true QTN loci — an **oracle** | no |
| `gea_best` | agreement of ≥2 of LFMM/RDA/EMMAX, 5 kb window | yes |
| `gea_strict` | agreement of all three | yes |
| `gea_union` | any of the three | yes |

**Balancing.** Only (replicate × garden × method) cells where **all four panels** are present
were used, so a panel difference cannot be an artefact of which replicates happen to carry
which panel. Per-replicate medians are taken first, then the median across replicates, so a
replicate with more gardens does not outweigh one with fewer.

**3600 paired cells: 30 replicates × 30 landscape gardens × 4 methods.**
Architecture strata are exactly balanced: 10 oligogenic / 10 mod-polygenic / 10 highly-polygenic.

---

## 2. Headline result

**GEA-identified panels match the QTN oracle.**

| panel | median τ | % of oracle | beats oracle in |
|---|---|---|---|
| `gea_best` (≥2 methods) | **−0.7041** | **100.8 %** | 40 % of replicate × method |
| `gea_strict` (all 3) | −0.6896 | 98.7 % | 44 % |
| `gea_union` (any) | −0.6376 | 91.3 % | 17 % |
| `adaptive` (oracle) | −0.6985 | — | — |

The agreement-based panel is **statistically indistinguishable from the true QTN set**, and
median-wise slightly ahead of it. Requiring agreement is what matters: dropping to `union`
(any single method) costs ~9 points of oracle performance.

**This question cannot be posed in the source design.** Lind & Lotterhos's "adaptive" panel
*is* the true QTN set, so their comparison is oracle-vs-neutral. Here the candidate panels are
built by a GEA scan with no access to truth, which is the situation on real data.

## 3. It holds across genetic architecture

`gea_best` vs oracle, per stratum (all four offset methods pooled, n = 40 cells each):

| architecture | % of oracle | beats oracle |
|---|---|---|
| mod-polygenic | 101.5 % | 35 % |
| oligogenic | 100.7 % | 40 % |
| highly-polygenic | 99.5 % | 45 % |

Architecture does not decide *whether GEA panels work* — they track the oracle in all three
strata. (Architecture did strongly govern the *curated-vs-uncurated* contrast in the 14-seed
run; that contrast needs the `all`/`neutral` panels and is therefore still pending.)

## 4. Method ranking is stable

Median τ by method (per-replicate medians, `adaptive` panel):

| method | median τ |
|---|---|
| GFoffset | −0.730 |
| LFMM2offset | −0.719 |
| RDA-uncorrected | −0.668 |
| RDA-corrected | −0.664 |

RDA-corrected remains the weakest, reproducing the 14-seed finding that
`condition_pcs: 0` is the right default for offset.

## 5. Cohort check — dataset size does not drive the result

The 18 Group 1 replicates are ~2.5× larger than the original 14 (median 24 720 vs 10 040 SNPs);
the selection rule stratified on architecture and R², not on dataset size. Effect on the
headline metric is small:

| cohort | median SNPs | median τ (`gea_best`, GF) | n |
|---|---|---|---|
| group1_expansion | 24 720 | −0.7271 | 18 |
| initial | 10 040 | −0.7446 | 12 |

Because the panel comparisons are **paired within replicate**, cross-replicate size differences
cannot bias them. Size remains a covariate worth reporting, not a confound to correct for.

---

## 6. Not yet answered (running)

1. **The floor.** `neutral` and `random_matched` (size-matched random panels) are the controls
   that prove curation beats an arbitrary set of the same size. Coverage is currently 27 and 26
   replicates and not balanced, so no floor comparison is made here.
2. **`all` vs curated.** Whether marker curation pays off at all — the strongest claim of the
   14-seed run, and the one that depends most on architecture.
3. **The remaining 70 gardens.** This report uses 30 of 100 landscape gardens. The full
   reciprocal-transplant design is phase 2.

Reasons 1 and 2 are the same reason: `all` (10–26 k SNPs) and `neutral_all` (~5 k) dominate
Gradient Forest build cost, so they were deferred to get this result hours sooner.

## 7. Provenance

- Every offset value produced by the scenario-aware pipeline, validated byte-for-byte against
  the previous per-garden driver on 1176 files (max abs diff 0.000e+00).
- Fitness identity recovered exactly on all 32 replicates (R² ≥ 0.999999).
- Local adaptation reproduces the deposit's own values to within 0.022.
- Numbers behind this report: `benchmarks/mvp_eval/offset09/phase1_seed_medians.tsv`,
  `garden_performance.tsv`.

---

# Addendum — single-method panels: does combining methods help?

Added at request: one panel per GEA method (`solo_lfmm`, `solo_rda`, `solo_emmax`), each built
from the **same scan** at the **same verbatim operating point** as the combined panels. The only
thing that differs is the combination rule, so any gap is attributable to combining.

Balanced and paired exactly as above: **3600 cells, 30 replicates x 30 gardens x 4 methods**,
seven panels present in every cell.

## A1. The result — combining does NOT beat the best single method

| panel | median tau | % of oracle | beats oracle |
|---|---|---|---|
| `gea_best` (>=2 agree) | **-0.7041** | 100.8 % | 40 % |
| **`gea_emmax_only`** | **-0.6982** | **100.0 %** | 48 % |
| `gea_strict` (all 3) | -0.6896 | 98.7 % | 44 % |
| `gea_rda_only` | -0.6889 | 98.6 % | 42 % |
| `gea_union` (any) | -0.6376 | 91.3 % | 17 % |
| `gea_lfmm_only` | -0.6217 | 89.0 % | 13 % |

Head-to-head, combining vs each single method (paired per replicate x method):

| contrast | combo better in |
|---|---|
| `gea_best` vs `gea_lfmm_only` | **65 %** |
| `gea_best` vs `gea_rda_only` | 44 % |
| `gea_best` vs `gea_emmax_only` | 40 % |
| `gea_strict` vs `gea_emmax_only` | 34 % |

**Read plainly: combining beats LFMM alone, and does not beat EMMAX alone or RDA alone.**
`gea_best` edges `gea_emmax_only` on the median (-0.7041 vs -0.6982) but wins only 40 % of
paired cells -- that is a coin flip, not an improvement. The honest statement is that
combining **protects against picking the wrong single method**, not that it beats the right one.
Since which method is "right" is unknowable without truth, that protection still has value --
but it is a weaker claim than "combining improves accuracy" and should not be written as one.

## A2. What actually predicts offset accuracy: precision, not recall

Ranking the panels by size against accuracy:

| panel | median SNPs | median tau |
|---|---|---|
| `gea_emmax_only` | 12 | -0.6982 |
| `gea_strict` | 39 | -0.6896 |
| `adaptive` (oracle) | 43 | -0.6985 |
| `gea_rda_only` | 99 | -0.6889 |
| `gea_best` | 152 | -0.7041 |
| `gea_lfmm_only` | 197 | -0.6217 |
| `gea_union` | 225 | -0.6376 |

**Spearman rho(panel size, tau) = +0.50 — bigger panels are systematically WORSE.**
The two largest panels are the two worst; a 12-SNP EMMAX panel matches the QTN oracle.
Maximising recall (`union`, `solo_lfmm`) actively hurts: diluting a panel with false positives
degrades the offset model more than the extra true positives help it.

This reframes the pipeline's practical recommendation. The useful knob is not "combine more
methods" but "use a high-precision operating point" -- and EMMAX at p < 1e-04 supplies one
essentially for free.

## A3. Caveat that must travel with A2

Panel size and method identity are **confounded** in this comparison: `solo_emmax` is small
*because* EMMAX at 1e-04 calls few SNPs. The size-matched random floor (`rand_solo1` and
friends) is the control that separates "small and precise" from "small", and it is part of the
deferred reference set. Until it lands, A2 is a strong association, not a demonstrated cause.

Note the confound does NOT flatter the result: the naive expectation is that more markers help,
and the data show the opposite, so this is not a case of a favourable panel being handed an
advantage by its size.
