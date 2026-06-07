# Bug: Per-method Manhattan overlay crashes the app — `object 'region_id' not found`

**Reported:** 2026-06-06 · **Project that surfaced it:** WBDC_final (GEA) · **Severity:** High (uncaught error terminates the Shiny session/app)

## Summary
Opening the **Per-Method Details** panel on the GEA (or GWAS) tab throws an uncaught
error inside `build_manhattan_plotly()` and the app closes:

```
Warning: Error in eval: object 'region_id' not found
 128: paste0
 125: [.data.table
 123: build_manhattan_plotly
 121: func (renderPlotly)
 101: output$gea-method_manhattan-overlay
   3: runApp
```

## Root cause
`build_manhattan_plotly()` unconditionally requires a `region_id` column on `sig_snps`,
but the per-method overlay supplies sig SNPs that don't have one.

- **Consumer** — `scripts/adaptogene.app/R/fct_manhattan.R:393–399` (and `:421`):
  ```r
  sig_snps[, hover_text := paste0(
      "SNP: ", SNPID, ... "Region: ", region_id, ...   # <- region_id hard-referenced
  )]
  ...
  customdata = ~region_id
  ```
- **Producer (per-method)** — `compute_method_sigsnps_cached()` at
  `scripts/adaptogene.app/R/fct_data_loading.R:667–689` returns only
  `SNPID, chr, pos, pvalue, method, trait`. `add_cum_pos()`
  (`fct_manhattan.R:468–473`) then adds `cum_pos` + `log10p`, but **never `region_id`**.
- **Wiring** — `mod_gea.R:495–509` (`per_method_sigsnps_override`) passes that table
  straight into `build_manhattan_plotly` via
  `mod_manhattan_overlay_server("method_manhattan", …)`. Identical wiring in
  `mod_gwas.R:574–588`.
- **Why the Combined Manhattan is fine** — it sources sig SNPs from
  `compute_interactive_sigsnps()` (`fct_combine.R`), whose output schema includes
  `region_id` (empty skeleton at `fct_combine.R:181–191`). Only the **per-method**
  path lacks the column.

## Affected
GEA **and** GWAS → "Per-Method Details" → per-method Manhattan
(`output$<mod>-method_manhattan-overlay`). Combined Manhattan, QQ, and region detail
are unaffected.

## Why it surfaced only after restarting the container
Before the restart, a stale **app-level `bindCache`** held `k_best=6` (the config
`config_WBDC_final.yaml` was finalized with `k_best: 5.0` ~30 h after the app process
started; until then `resolve_k_best()` fell back to `max(find_k_values)` = 6). With
`k_best=6` the per-method coords path (`…_K6_coords.json`) was missing, so the overlay
short-circuited to the "plot not available" placeholder and never reached
`build_manhattan_plotly`. After restarting (`k_best=5`, coords present), the per-method
overlay actually renders → hits the missing-column crash.

## Reproduce
1. Open a project with GEA results (e.g. WBDC_final), ensure `k_best` matches the
   produced outputs (K5) so the per-method plot actually renders.
2. GEA tab → expand **Per-Method Details**.
3. App errors with `object 'region_id' not found` and the session dies.

## Suggested fixes (dev)
1. **Make `region_id` optional in `build_manhattan_plotly`** (preferred — single source
   of truth, protects every caller): guard the hover/`customdata` with
   `if ("region_id" %in% names(sig_snps))`, else omit the "Region:" line and set
   `customdata` to `SNPID`/`NA`. *Alternative:* have `per_method_sigsnps_override`
   add `dt[, region_id := NA_character_]` to match the combined schema — but per-method
   SNPs have no real region assignment, so the guard is cleaner.
2. **Fail soft on render errors** — wrap the `renderPlotly` body (or the
   `build_manhattan_plotly` call) in `tryCatch` so one bad output degrades to a
   placeholder instead of terminating the whole Shiny process. A per-output render
   error killing the app is its own robustness bug.

## Related robustness issues observed in the same session
- **Stale `k_best` via app-level `bindCache`** on `project_data`
  (`app_server.R:66–71`): an external config/output change isn't picked up without a
  process restart or an in-app run bumping `project_data_trigger`; browser reloads don't
  help (cache is shared across sessions). Consider a session-scoped cache or keying on
  the config file mtime.
- **`resolve_k_best()` silent `max(K)` fallback** (`utils_helpers.R:43–48`): when
  `config.sNMF.k_best` is momentarily absent it silently picks the largest available K
  instead of failing/flagging — this is how `K6` got cached.

## Workaround (until patched)
Don't expand **Per-Method Details** — Combined Manhattan, region selection, genes,
enrichment, and QQ all work. Per-method inspection is the only thing that triggers the
crash.
