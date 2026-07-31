# Láruson et al. 2022 — Gate G1 confirmed dataset facts

Confirmed empirically from the full Dryad archive (`doi_10_5061_dryad_x95x69pkk__v20220525.zip`,
1.14 GB, user-provided) — `data/laruson/README.txt`, `seeds_source_R.txt`, and one full replicate
(`Case1_2889863491989_*`). This supersedes the placeholder/guessed values in
`docs/laruson_automation_roadmap.md` and `docs/laruson_benchmark.md`.

## 1. Archive layout (corrects `laruson_benchmark.md` §2/§3c)

**The "4 Cases" are NOT Neutral / Monogenic / Polygenic / Additional as originally guessed.**
Per `seeds_source_R.txt`, all four are multilocus (polygenic) scenarios, distinguished by
**landscape shape**, not architecture:

| Case | SLiM script | Landscape |
|------|-------------|-----------|
| 1 | `ML_WF_Cline.slim` | Simple monotonic cline — closest to the "1D gradient" mental model |
| 2 | `ML_WF_MountainRange.slim` | Non-monotonic ("mountain") environmental gradient |
| 3 | `ML_WF_Case3.slim` | Undocumented beyond script name (paper methods needed for detail) |
| 4 | `ML_WF_Case4.slim` | Undocumented beyond script name |

`seeds_Neutral.txt` and `seeds_SingleLocus.txt` are seed lists only — **the true Neutral and
Single-Locus scenarios have no corresponding data files in this archive** (no `.vcf.gz`, no
`_ind.txt`, etc.). Only the 4 multilocus Cases (40 replicates total, 10 seeds each) ship data.

**Chosen replicate for the adapter pilot: `Case1_2889863491989`** (Case 1 / Cline — simplest,
most interpretable landscape, matches the benchmark doc's original "1D gradient" framing even
though the underlying model has 2 traits / 2 environmental axes, see below).

No `src/` analysis scripts are included in this archive (README describes them, but the actual
`.R`/`.py`/`.slim` files are not shipped) — the row-order join described in §3 could not be
cross-checked against the authors' own code and is stated as an explicit, validated-in-practice
assumption.

## 2. Landscape structure

**100 demes (subpopulations) on a 10×10 grid**, confirmed from `Case1_..._ind.txt`: unique
`(subpop, opt0, opt1)` triples show `opt0` cycling through 10 evenly-spaced values in
`{-1, -0.778, ..., 1}` (step 2/9) as `subpop` increases 1→10, then `opt1` stepping to the next
of the same 10 values for subpops 11-20, 21-30, etc. — i.e. `row = (subpop-1) %% 10`,
`col = (subpop-1) %/% 10`, `opt0 = -1 + 2*row/9`, `opt1 = -1 + 2*col/9`. Confirmed for this
replicate; the converter derives row/col from the actual unique pairs per replicate rather than
hardcoding the formula, since Case 2-4 landscapes may differ.

**Two independent environmental optima axes** (`opt0`, `opt1`) and **two traits** (`phen0`,
`phen1`) — this is a 2-env/2-trait model, not literally a single 1D gradient, despite Case 1's
"Cline" name (the cline applies per-axis).

## 3. Per-replicate file contents (confirmed against `Case1_2889863491989_*`)

### `.vcf.gz` — unfiltered, non-subsampled
- 10,000 samples, single contig (`##contig=<ID=1,length=500000>`), 33,541 unfiltered variants (this replicate)
- Sample column names are **not** a simple `30000000 + indID` offset (checked and rejected — the
  pattern held for the first ~5 names then diverged)

### `_ind.txt` — 10,000 rows + header, columns: `indID indSubpopIndex subpop phen0 phen1 opt0 opt1 fitness`
- `indID`: 0-based row counter (not a VCF join key)
- `indSubpopIndex`: **misleadingly named** — empirically only ~100 distinct (rounded) values across
  10,000 rows; it is a PySlim-internal subpop-level index, effectively redundant with `subpop`,
  **not a per-individual identifier**. Do not use it to join to VCF samples.
- `subpop`: 1-100, **sorted** — rows are grouped in contiguous blocks of exactly 100 per subpop
  (rows 1-100 = subpop 1, 101-200 = subpop 2, ..., confirmed for this replicate)
- **VCF join key: row order.** `ind.txt` row *i* (0-based, after header) = VCF sample column *i*
  (0-based, after the 9 fixed VCF columns). Validated: both have exactly 10,000 entries, and the
  README states the same PySlim script emits the paired VCF + subset table together. **This is an
  assumption, not independently re-derived from the authors' code** (no `src/` scripts shipped) —
  the converter hard-validates `n_vcf_samples == n_ind_rows` and stops loudly if they ever diverge.

### `_causal_mutations_pos_filtered.txt` — headerless, single column, 102 rows (this replicate)
- Values are **exact 1-based VCF `POS`** on the single contig — verified 102/102 exact matches
  against this replicate's VCF `POS` column. No chromosome column needed (always the VCF's one
  contig). **No effect-size or per-trait-assignment column** — contra the roadmap's assumed
  4-column schema (`chr, pos, effect_size, category`), this file gives position only. `causal_loci.tsv`
  will carry `chr, pos, category` (no `effect_size`); `category` is derived by the converter
  (causal / linked_neutral via a TP-window / background_neutral), not read from source data.

### `_ML_WF_CG_sum_Gen.txt` — 33 rows (per ~100-generation snapshot) + header
- Columns: `generation sympatry allopatry local_adaptation mean_pheno0 mean_pheno1 corr_pheno0_opt0 corr_pheno1_opt1`
- **This is a population-average scalar time series, NOT a per-deme×garden reciprocal-transplant
  fitness matrix.** `sympatry` = mean home-deme fitness, `allopatry` = mean fitness when
  transplanted elsewhere (aggregated across all other demes, not itemized per garden),
  `local_adaptation = sympatry - allopatry`. Early generations show `sympatry == allopatry`
  (pre-adaptation burn-in); by the final generation (3300, this replicate) they diverge
  meaningfully (0.803 vs 0.535, local_adaptation 0.268).
- **Consequence for the (deferred) eval harness:** there is no ready-made home×garden fitness
  table in this file set. A real Axis-1 (offset-vs-fitness) design will need to be rebuilt from
  `_Freq_ML_WF.txt`'s per-generation, per-deme `P#_fit` columns (100 demes × 33 generations) as a
  substitute for "sympatric fitness per deme", or some other reconstruction — not attempted here.

### `_Freq_ML_WF.txt` — 33 rows × 613 columns
- 13 simulation-parameter columns (`m, n, u, r, mean_Eff, var_Eff, var_Opt, Env_rate, Burnin1,
  Burnin2, Env_shift, Total_Gen, Generation`) + 6 blocks of 100 per-deme columns each
  (`P{1..100}_fit`, `P{1..100}_Freq`, `P{1..100}_Phen1`, `P{1..100}_Phen2`, `P{1..100}_Env1`,
  `P{1..100}_Env2`). Not used by the adapter converter (per-deme env is already available, more
  simply, from `ind.txt`'s `opt0`/`opt1`); kept as a future source for the eval harness.

## 4. Converter design decisions this implies

- **Sites** = `subpop` (1-100). **Samples** = VCF sample column names, taken verbatim.
- **Predictors**: 2 (`bio_1` = `opt0`, `bio_2` = `opt1`).
- **Coordinates**: `latitude = row * 2.0`, `longitude = col * 3.0` (row/col derived per-replicate
  from the actual unique `opt0`/`opt1` grid, not hardcoded) — valid in-range lat/lon, well-spaced
  (≥2° apart) so `Climate.custom.grid_resolution` can safely be small (e.g. 0.5°) without any
  two-sites-one-cell collision.
- **`environments_future.tsv`**: no real "future"/garden concept exists in this static landscape —
  built as `environments_present` + a fixed synthetic shift, explicitly labeled a placeholder to
  smoke-test the maladaptation machinery only. Real reciprocal-transplant evaluation is deferred
  to the eval-harness phase (which needs a different construction per §3 above).
- **`fitness_timeseries.tsv`** (renamed from the roadmap's assumed `fitness.tsv` — different
  shape): raw pass-through of `_ML_WF_CG_sum_Gen.txt`, archived for future eval-harness use, not
  consumed by this adapter-only pipeline run.
- **`TP_WINDOW_KB`**: no LD-decay estimate available without running the pipeline; defaulted to
  5 kb (1% of the 500 kb contig) as a starting value, clearly a placeholder pending real LD-decay
  analysis once `mode=processing`/`mode=structure` have run on the converted data.
