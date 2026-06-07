# Bug: WZA Manhattan rules fail — `NA/NaN/Inf in foreign function call (arg 6)` (scattermore)

**Reported:** 2026-06-07 · **Project that surfaced it:** WBDC_final (GWAS) · **Severity:** Medium
(the `mode=gwas` run aborts, but only 2 leaf plots are affected — every GWAS table and the
per-SNP Manhattans are complete and usable)

> **STATUS: FIXED (2026-06-07, user-authorized).** Both suggested fixes applied. WBDC_final
> `mode=gwas -R assoc_wza` re-ran 35/35 steps (exit 0); both plots regenerated; the
> regenerated WZA tables contain no 0/NA/negative. The previously-underflowed window now
> reads its **true** p = 8.48e-18 (−log10 = 17.07): `1 - pnorm()` was catastrophic
> cancellation rounding a representable value to 0, and `pnorm(lower.tail = FALSE)` recovers
> it — the `pmax` floor never had to trigger. Changed files: `scripts/compute_wza.R`,
> `scripts/R/utils/manhattan_utils.R` (`add_scatter_layer`), `scripts/plot_manhattan_combined.R`.
> The downstream guard also fixes a **latent** crash: WZA tables can legitimately contain NA
> windows, which scattermore equally cannot render (this dataset just had none).

## Summary
Two **WZA-regime** Manhattan jobs failed during `mode=gwas` for WBDC_final; the failure
aborted the run (`WorkflowError: At least one job did not complete successfully`). From the
master log `WBDC_final_logs/.run.log` (lines 1559–1693):

- `assoc_wza_manhattan` — jobid 118 — `source=GWAS, method=BLINK, trait=Tiller_Number`
- `assoc_wza_manhattan_combined` — jobid 131 — GWAS, EMMAX+BLINK pooled

Everything else completed: all other per-trait WZA Manhattans (EMMAX ×11, BLINK ×10) and the
entire non-WZA GWAS run finished (119/129 steps before the abort). These two outputs are leaf
plots — nothing depends on them.

Both rule logs
(`WBDC_final_logs/GWAS/wza_manhattan_BLINK_Tiller_Number_bonf_0.05.log` and
`WBDC_final_logs/GWAS/wza_manhattan_combined.log`) end with the **same** error:

```
INFO: Generating background Manhattan for Shiny overlay
Error in `scattermore::geom_scattermore()`:
! Problem while converting geom to grob.
ℹ Error occurred in the 1st layer.            # (5th layer for the combined plot)
Caused by error in `scatter_points_rgbwt()`:
! NA/NaN/Inf in foreign function call (arg 6)
Execution halted
```

Note the per-trait log shows the simple Manhattan **and** the QQ plot saving successfully —
it only dies at the **background PNG** step.

## Root cause (two levels)

### 1. Upstream: a WZA window p-value of exactly `0`
`scripts/compute_wza.R:235` (and the no-correction fallback `:240`):

```r
wza_p[valid] <- 1 - pnorm(win_zw$Z_W[valid], corr$pred_mean, corr$pred_sd)
...
wza_p        <- 1 - pnorm(win_zw$Z_W, mean(win_zw$Z_W, ...), global_sd)
```

For a strongly significant window the standardized `Z_W` is large enough that
`pnorm(Z_W, ...)` returns exactly `1.0` in IEEE double precision, so `1 - pnorm(...)`
underflows to **exactly `0`** (catastrophic cancellation). No flooring is applied to the
result.

**Verified in the data** — the `Tiller_Number` column of
`WBDC_final_results/GWAS/tables/methods/BLINK/BLINK_wza_K5.tsv` contains **one window with
p = 0**, and it is the *only* `0` / `NA` / negative value across all 11 trait columns
(scan of 4144 rows). The values are raw p-values, confirmed by the threshold math in the
log: `0.05 / 4144 windows = 1.2066e-5`, exactly the logged Bonferroni threshold.

### 2. Downstream: `-log10(0) = +Inf` reaches scattermore unguarded
`prepare_manhattan_data()` at `scripts/R/utils/manhattan_utils.R:49`:

```r
log10p = -log10(.data[[pval_col]])   # -log10(0) -> +Inf, no finite guard
```

- **Simple Manhattan + QQ survive.** The `+Inf` window is significant
  (`Inf >= threshold_log10` → `TRUE`, `plot_manhattan.R:99`), so it is routed to the
  `geom_point(data = filter(is_significant))` layer (`plot_manhattan.R:127`). ggplot2's
  `geom_point` silently drops non-finite rows (with a warning), so those plots render.
- **Background PNG fails.** The background passes **all** rows (threshold-independent — the
  complete static cloud) into scattermore:
  `add_scatter_layer(data = plot_df, ...)` (`plot_manhattan.R:229`) →
  `scattermore::geom_scattermore(...)` (`manhattan_utils.R:114`). scattermore's C routine
  `scatter_points_rgbwt()` does **not** tolerate non-finite coordinates and aborts with
  `NA/NaN/Inf in foreign function call (arg 6)` → non-zero exit → rule fails and Snakemake
  removes the partial outputs.
- **Combined background fails identically** at `plot_manhattan_combined.R:195`
  (`geom_scattermore(data = df_t_bg, ...)`); the combined plot pools all traits/methods and
  therefore includes the same Tiller_Number `+Inf` point ("5th layer" in its backtrace).

No finite guard exists anywhere in the plotting/util chain (confirmed by grep across
`plot_manhattan.R`, `plot_manhattan_combined.R`, `manhattan_utils.R`, `io_pvalues.R`).

## Why only `Tiller_Number`/BLINK and the combined plot failed
The crash needs a window p-value that underflows to `0`. Only the BLINK `Tiller_Number`
column had one. So:
- the per-trait **`assoc_wza_manhattan`** failed *only* for `Tiller_Number`/BLINK;
- the **`assoc_wza_manhattan_combined`** failed because it pools every trait/method and
  thus includes that one `+Inf` point.

The **per-SNP (non-WZA) regime is unaffected** because raw EMMAX/BLINK per-SNP p-values in
this dataset are never exactly `0` — it is the WZA Stouffer-Z → p conversion that underflows
for the strongest window. The non-WZA `assoc_manhattan_combined` (jobid 130) finished fine.

## Affected
- Rules `assoc_wza_manhattan` (WZA per-trait) and `assoc_wza_manhattan_combined`
  (`workflow/rules/_assoc_downstream.smk`), for **any** source / method / trait whenever a
  WZA window p underflows to `0`. This is general (GEA as well as GWAS), not specific to this
  dataset.
- **Not affected:** per-SNP-regime Manhattans (`assoc_manhattan_plot`,
  `assoc_manhattan_combined`), and all GWAS p-value / q-value / WZA / region / gene tables —
  the bad value is only a *display* problem.

## Reproduce
1. Run `mode=gwas` (or `gea`) with the WZA regime enabled on data where at least one window
   is strongly significant (large `Z_W`).
2. `compute_wza` emits a window with p = 0 for some trait/method.
3. `assoc_wza_manhattan` for that trait/method (and the WZA combined plot) die at the
   "Generating background Manhattan for Shiny overlay" step with the scattermore error above.

## Suggested fixes (dev — NOT applied; this session is USE-only)

**Recommend doing both** — (1) fixes the value at its source for every consumer; (2) makes the
plotting rules fail-soft so a single bad value can never abort the whole run again (the same
lesson as the prior `region_id` bug report).

1. **Upstream correctness fix (preferred root fix) — `scripts/compute_wza.R:235,240`:**
   Replace `1 - pnorm(Z_W, m, s)` with the upper-tail form to avoid catastrophic
   cancellation, and floor the result so it can never be exactly `0`:
   ```r
   wza_p <- pnorm(win_zw$Z_W, m, s, lower.tail = FALSE)   # instead of 1 - pnorm(...)
   wza_p <- pmax(wza_p, .Machine$double.xmin)             # ~2.2e-308 -> -log10(p) ~ 307.6
   ```
   This gives a tiny **positive** p (finite `-log10(p)`) for the strongest window and fixes
   every downstream consumer at once — plots, significant-SNP selection, and the Shiny
   overlay. It is also the statistically correct way to evaluate an extreme upper-tail
   probability.

2. **Downstream robustness fix (defense-in-depth, single source of truth) —
   `scripts/R/utils/manhattan_utils.R`:** guard non-finite `log10p` before it reaches
   scattermore. Cleanest single location is `add_scatter_layer()` (line 109), which every
   background caller routes through:
   ```r
   data <- data[is.finite(data$log10p) & is.finite(data$pos_cum), ]
   ```
   To honor the project rule that *every* significant SNP keeps a faint ghost dot in the
   background, **cap rather than drop** — e.g. clamp `log10p` to the finite plot ceiling
   (`y_max`) instead of removing the row. Optionally also cap inside
   `prepare_manhattan_data()` (line 49) so `-log10(p)` is always finite for all callers.

## Workaround (until patched)
None needed for analysis — the WZA Manhattan plots are optional leaf outputs. All GWAS
p-value/q-value tables, WZA window tables, regions, genes, enrichment inputs, and the
per-SNP Manhattans are complete and correct. Only the two WZA Manhattan PNG/SVG/coords
outputs are missing.

## Re-run after fix
- Applying **fix (2) only:** the next `mode=gwas` run regenerates just the two removed plots
  (their inputs are unchanged), e.g.
  `--forcerun assoc_wza_manhattan assoc_wza_manhattan_combined`.
- Applying **fix (1)** (changes `compute_wza.R` output): regenerate the WZA tables first so
  the `0` is replaced by a finite tiny p — `--forcerun assoc_wza` — and the dependent WZA
  Manhattans rerun automatically.
