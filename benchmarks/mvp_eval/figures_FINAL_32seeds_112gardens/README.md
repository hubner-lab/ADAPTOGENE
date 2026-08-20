# FINAL figures — full simulation set

**Basis: 32 replicates x 112 gardens x 4 offset methods x 12 marker panels, all computed.**
30 non-control replicates enter every aggregate (2 Est-Clines controls excluded by design).
Gardens: full reciprocal transplant, 100 landscape + 12 climate-novelty.

Superseded folders (partial data, do not cite): `../figures/`, `../figures_prelim/`.

## MAIN FIGURE — GEA marker-selection efficiency

**`MAIN_FIGURE_gea_marker_efficiency.png`** (identical to `abs_D_composition_counts.png`;
duplicated under the MAIN_ name so it is unambiguous which panel carries the headline).

What it shows: for each obtainable selection rule, the ABSOLUTE number of markers delivered,
split into causal / linked / background, faceted by genetic architecture and ordered by usable
markers (causal + linked). Reference panels (true QTNs, size-matched random) are deliberately
excluded — they are not candidates a user can choose, and including them compressed the axis.

Why this is the headline rather than the accuracy plot: it shows the MECHANISM. Percentage
composition made `EMMAX only` look excellent (20% causal, the purest panel) while hiding that
this is 20% of ~12 markers, falling to 3 usable markers on highly-polygenic traits. Counts make
the selection rules comparable on the quantity that actually drives offset accuracy.

Read alongside `abs_E_size_vs_architecture.png`, which pairs the same size trend with accuracy:
`≥2 agree` grows 140→150→158 markers as polygenicity rises and holds 0.711–0.723 accuracy;
`EMMAX only` shrinks 22→10→5 and drops to 0.631; `all 3 agree` shrinks 80→36→18 and falls to
0.513 — below a size-matched random panel, despite carrying zero background noise.



Metric: prediction accuracy = -Kendall's tau vs common-garden fitness, HIGHER IS BETTER.
Per-replicate medians first, then across replicates. Panel contrasts are paired within replicate.

Precision/recall are recomputed per replicate from its own truth table on a chr:pos key --
NOT from snp_sets_summary.tsv, which is overwritten per --sets run and covers different
replicate counts per panel.

| file | shows |
|---|---|
| explain_A_accuracy | which marker set predicts best |
| explain_B_recall | more true loci found => worse prediction |
| explain_C_composition | causal / linked / background make-up of each set (SHARES — see caveat above) |
| **MAIN_FIGURE_gea_marker_efficiency** | **headline: absolute marker counts per selection rule x architecture** |
| abs_D_composition_counts | same figure, original name |
| abs_E_size_vs_architecture | panel size and accuracy vs architecture |
| abs_F_usable_vs_accuracy | usable markers vs accuracy |
| abs_table_by_architecture.tsv | the absolute-count table behind D/E/F |
| dist_accuracy | per-replicate spread of accuracy, one dot per replicate |
| dist_recall, dist_precision, dist_markers | same, for the other quantities |
| fig1-fig3 | pipeline throughput (old per-garden driver vs scenario pipeline) |
| fig4, fig5 | tau by panel x method; tau by architecture stratum |
| panel_summary.tsv, panel_distribution_table.tsv, panel_paired_tests.tsv | numbers behind the panels |
