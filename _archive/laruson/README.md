# Láruson et al. 2022 — archived 2026-08-01

Retired as a GEA benchmark. **Do not revive for GEA testing.** See the "Retired: Láruson et al.
2022" section in `CLAUDE.md` for the full reasoning, and `benchmarks/laruson_eval/FINDINGS.md`
here for the complete measured results.

## Why it was retired

1. **No population structure.** SLiM parameters `m = 0.2`, `n = 100` → Nm = 20, expected
   Fst ≈ 0.012. Measured PC1 = 0.38% of variance. PC1/PC2 are a rotation of the two
   environmental axes (R²_env ≈ 0.84), PC3 onward is noise. Structure correction can therefore
   only remove signal, and the dataset cannot exercise that part of the pipeline at all. It was
   designed to test Gradient Forest genetic-offset prediction, not GEA.
2. **Truth set incompatible with our MAF filter.** The 102 causal positions are already
   MAF-filtered by the authors at MAF ≥ 0.01 (min causal MAF = exactly 0.010000), while the
   shipped VCF is unfiltered. At `Filter.maf: 0.05` only 28/102 are testable.

## Contents

| path | what |
|---|---|
| `results/` | `LARUSON_results` (10k), `LARUSON1K_results` (1k subsample), `LARUSON1Kmaf001_results` (aborted maf 0.01 arm) |
| `logs/` | matching `*_logs` directories |
| `data/` | Dryad archive zip, `laruson/` (raw replicate), `laruson_converted/` (10k pipeline input), `laruson_1k/` (subsample input) |
| `benchmarks/laruson_eval/` | **FINDINGS.md** + all scored TSVs, Pareto tables, precision-recall plots |
| `benchmarks/laruson_sweep/` | per-rung RDA p-value tables (condition_pcs 0..10) |
| `benchmarks/convert_laruson.R` | Dryad → pipeline-input converter (dataset-specific) |
| `benchmarks/subsample_laruson.R` | balanced N-per-site subsampler (dataset-specific) |
| `docs/` | `laruson_benchmark.md`, `laruson_dataset_notes.md`, `laruson_automation_roadmap.md` |
| `configs/` | `config_LARUSON.yaml`, `config_LARUSON1K.yaml`, `config_LARUSON1Kmaf001.yaml` |

## Deliberately NOT archived — reuse these for the next benchmark

These stayed in `benchmarks/` because they are dataset-agnostic: they take any wide p-value
table (`SNPID chr pos <traits...>`) plus a truth table with `chr, pos, category` columns
(`causal` / `linked_neutral` / `background_neutral`).

- `lib_detection.R` — scoring core: TP/FP accounting, window matching, dual recall denominators,
  genomic-control recalibration
- `eval_detection.R` — score one p-value table (`--mode=threshold` or `--mode=rank`)
- `sweep_thresholds.R` — full adjust×threshold grid plus union / intersection / ≥2-of-3 /
  rank-sum combine rules, Pareto frontier + PR plot
- `sweep_rda.sh` — RDA `condition_pcs` ladder retaining per-SNP p-values
- `score_ladders.sh` — score every pregea LFMM-K / EMMAX-n_pcs rung against truth

They reuse the pipeline's own `compute_pval_threshold()`, so threshold calls reproduce
`find_sig_snps.R` exactly.

## Transferable lessons (independent of this dataset)

- `EMMAX.params.n_pcs` defaults to `sNMF.k_best`. Where structure and environment are collinear
  this silently destroys the method, and **λ_gc stays ≈ 1.0 at every rung**, so no QQ plot or
  calibration-based rule can detect it. Only scoring against known causal loci exposes it.
- Cap BLAS threads (`OMP_NUM_THREADS=4`) when running `scripts/rda.R` on a shared machine.
  Unbounded, one rung burned 3.5 CPU-hours in 13 min without finishing; capped, 70 s.
- `axes=<fixed int>` does **not** skip `anova.cca` — `rda.R:318-329` runs full/by-axis/by-margin
  unconditionally, the latter two unparallelized. Permutations drive RDA wall-clock regardless.
- A third pyslim VCF-writer quirk beyond the two in `docs/laruson_dataset_notes.md` §3: three
  rows carry a GT allele index of 2 while declaring one ALT, aborting `bcftools` AC/AN
  recalculation. Use `bcftools view -I`.
