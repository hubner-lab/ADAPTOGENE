# Discussion notes — observations to carry into the manuscript

**Status: PROVISIONAL.** All of this rests on 30 SS-Mtn replicates × 30 landscape gardens
(phase-1 garden set) with the adaptive + GEA-detected panels. Each observation below must be
re-checked against (a) the remaining 70 gardens and (b) the size-matched random floor before
it goes into a manuscript. Recorded now so the reasoning is not reconstructed from memory later.

Date: 2026-08-08. Numbers: `figures/table_absolute_counts.tsv`, `offset09/phase1_seed_medians_solo.tsv`.

---

## Note 1 — Panel size must scale with genetic architecture; "cleanest" is not "best"

The intuitive reading of the pooled result is "less noise wins". Per-architecture accuracy
contradicts it:

| panel | markers (oligo / mod / highly) | oligogenic | mod-polygenic | highly-polygenic |
|---|---|---|---|---|
| all 3 agree | 80 / 36 / 18 | **0.709** | **0.748** | **0.553** ← worst obtainable |
| EMMAX only | 22 / 10 / 5 | 0.704 | 0.739 | 0.625 |
| ≥2 agree | 140 / 150 / 158 | 0.704 | 0.703 | **0.714** |
| union | 218 / 225 / 235 | 0.547 | 0.565 | 0.691 |

`all 3 agree` carries **0 background markers at every architecture** and is still the worst
obtainable panel on highly-polygenic replicates. Requiring all three methods to agree leaves
~18 markers to span a ~410-locus gradient; below some size, cleanliness stops compensating.
`EMMAX only` shows the same failure at 5 markers.

Note also that `union` moves the *opposite* way — it is worst on oligogenic (0.547) and
respectable on highly-polygenic (0.691). Its 31-59 noise markers are fatal when only 7 causal
loci exist and tolerable when 410 do.

**Manuscript claim this supports:** the operating point of a GEA scan should be chosen relative
to the expected polygenicity of the trait, not fixed. `≥2 agree` is the only panel at or near
the top in all three strata (0.704 / 0.703 / 0.714) because its size tracks the architecture
(140 → 158 markers) while noise stays low. That is an argument from **robustness across
architectures**, which is stronger than any single pooled median and is what should be written.

**Check before use:** does the ordering hold on the full 112-garden design? `all 3 agree` on
highly-polygenic is the specific number most likely to move, since it rests on ~18 markers.

---

## Note 2 — GEA panels can beat the QTN oracle; linked markers are not merely a substitute

On two of three architectures the obtainable panels **exceed** the true-causal-loci panel:

| architecture | oracle | best obtainable | margin |
|---|---|---|---|
| oligogenic | 0.702 | 0.709 (`all 3 agree`) | +0.007 |
| mod-polygenic | 0.678 | **0.748** (`all 3 agree`) | **+0.070** |
| highly-polygenic | 0.719 | 0.714 (`≥2 agree`) | −0.005 |

The mod-polygenic margin is not marginal. Two candidate explanations, not yet separated:

1. **Linked markers can be better predictors than the causal loci themselves.** A causal variant
   may sit at low frequency or be poorly tagged; a linked marker at intermediate frequency can
   track the environmental gradient with less sampling noise. If so, the oracle is not an upper
   bound on offset accuracy — only on detection.
2. **The oracle panel is heterogeneous in size** (28-98 causal loci across the 10 mod-polygenic
   replicates) and may be over- or under-sized for the offset model on individual replicates,
   whereas the GEA panels are size-stabilised by construction (fixed operating points).

**Manuscript claim this supports:** "recover the causal loci" is the wrong objective function
for genomic offset. The objective is to assemble markers that span the adaptive gradient, and
detection metrics (precision/recall against QTNs) are a poor proxy for it — see the decoupling
in Note 3.

**Check before use:** separate the two explanations. Explanation 2 is testable now by
regressing the oracle's per-replicate accuracy on its own size; explanation 1 needs the
allele-frequency spectrum of causal vs linked markers.

---

## Note 3 — Detection and prediction are largely decoupled (supporting context for Notes 1-2)

`≥2 agree` finds **17 of 410** causal loci on highly-polygenic replicates and predicts as well
as the oracle holding all 410. `EMMAX only` finds **1 of 410** and reaches 0.625. Meanwhile
`union` finds the most causal loci of any obtainable panel at every architecture and ranks
near the bottom at two of three.

Across panels: rho(background-noise fraction, accuracy) = **-0.63**;
rho(recall, accuracy) = **-0.79**. Recall is *negatively* associated with prediction accuracy.

This is why the report leads with prediction accuracy rather than precision/recall, and why the
precision/recall tables are presented as composition (causal / linked / background) rather than
as a detection score.

---

## Open items that gate all three notes

1. **Size-matched random floor** (`rand_*` panels) — separates "small and clean" from merely
   "small". Without it, Note 1's size argument is an association.
2. **Remaining 70 gardens** — the phase-1 set is 30 of 100 landscape gardens.
3. **`all` and `neutral_all` panels** — needed for the curated-vs-uncurated contrast, which the
   14-seed run found to be strongly architecture-dependent.
