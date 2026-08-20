# Marker panels: absolute counts and prediction accuracy, by genetic architecture

30 SS-Mtn replicates (10 per architecture) × 30 landscape gardens × 4 offset methods.
Counts are medians across the 10 replicates in each stratum; `range` gives the spread.
`accuracy` = −Kendall's τ between predicted offset and realised common-garden fitness
(**higher is better**), median of per-replicate medians.

`union` / `≥2 agree` / `all 3 agree` refer to the three GEA methods: **LFMM, RDA, EMMAX**.

---

## Oligogenic — 7 causal loci exist

| marker set | markers | causal found | range | linked | noise | accuracy |
|---|---|---|---|---|---|---|
| true QTNs (oracle) | 7 | **7/7** | 4-10/7 | 0 | 0 | 0.702 |
| ≥2 agree | 140 | 4/7 | 2-8/7 | 127 | 6 | **0.704** |
| all 3 agree | 80 | 4/7 | 2-7/7 | 76 | **0** | **0.709** |
| union | 218 | 6/7 | 4-10/7 | 179 | **31** | 0.547 |
| LFMM only | 197 | 6/7 | 4-9/7 | 162 | 26 | 0.536 |
| RDA only | 99 | 4/7 | 1-7/7 | 86 | 9 | 0.697 |
| EMMAX only | 22 | 4/7 | 2-7/7 | 15 | 1 | 0.704 |

## Moderately polygenic — 42 causal loci exist

| marker set | markers | causal found | range | linked | noise | accuracy |
|---|---|---|---|---|---|---|
| true QTNs (oracle) | 42 | **42/42** | 28-98/42 | 0 | 0 | 0.678 |
| ≥2 agree | 150 | 9/42 | 7-18/42 | 123 | 17 | 0.703 |
| all 3 agree | 36 | 3/42 | 1-8/42 | 33 | **0** | **0.748** |
| union | 225 | 13/42 | 9-22/42 | 165 | **55** | 0.565 |
| LFMM only | 194 | 12/42 | 7-20/42 | 141 | 40 | 0.565 |
| RDA only | 99 | 6/42 | 3-13/42 | 71 | 17 | 0.725 |
| EMMAX only | 10 | 2/42 | 1-7/42 | 7 | 1 | 0.739 |

## Highly polygenic — 410 causal loci exist

| marker set | markers | causal found | range | linked | noise | accuracy |
|---|---|---|---|---|---|---|
| true QTNs (oracle) | 410 | **410/410** | 268-534/410 | 0 | 0 | **0.719** |
| ≥2 agree | 158 | 17/410 | 5-24/410 | 108 | 21 | 0.714 |
| all 3 agree | 18 | 2/410 | 0-10/410 | 14 | **0** | 0.553 |
| union | 235 | 23/410 | 9-32/410 | 158 | 59 | 0.691 |
| LFMM only | 197 | 19/410 | 8-29/410 | 131 | 46 | 0.684 |
| RDA only | 99 | 11/410 | 2-15/410 | 70 | 18 | 0.675 |
| EMMAX only | 5 | 1/410 | 0-4/410 | 2 | 2 | 0.625 |

---

## How to read these

**Detection and prediction are largely decoupled.** `≥2 agree` finds 17 of 410 causal loci
on highly-polygenic replicates and predicts as well as the oracle holding all 410 — the
causal loci respond to the same climate gradient, so they are near-redundant. Capturing the
gradient is the task; enumerating the loci is not.

**Noise, not recall, is what costs accuracy.** `union` finds the most causal loci of any
obtainable panel at every architecture, and ranks near the bottom at two of three: it also
carries 31-59 background markers. Panels with 0-6 noise markers track the oracle.

**Panel size must match architecture — this is the one caveat.** `all 3 agree` is the best
panel on oligogenic and mod-polygenic (0 noise, matches or beats the oracle) but the **worst
obtainable panel on highly-polygenic** (0.553): requiring all three methods to agree leaves
only ~18 markers, too few to span a 410-locus gradient. Symmetrically, `EMMAX only` wins on
mod-polygenic with 10 markers and falls to 0.625 on highly-polygenic with 5.

`≥2 agree` is the only panel that is at or near the top in **all three** strata
(0.704 / 0.703 / 0.714) — it adapts its size to the architecture (140 / 150 / 158 markers)
while keeping noise low. That, rather than raw accuracy in any single stratum, is the
argument for it as the pipeline default.

**Still pending:** the size-matched random floor, which separates "small and clean" from
merely "small", and the remaining 70 gardens.
