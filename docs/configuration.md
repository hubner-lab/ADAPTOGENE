# Configuration Reference

All pipeline parameters are defined in `config.yaml`. This document covers every parameter group, GO enrichment details, and the output directory structure.

---

## Parameter Reference

### Project setup

```yaml
project_name: MyProject
Input:
  dir: data/
  vcf: your_genotypes.vcf
  metadata: your_metadata.tsv
  gff: your_annotation.gff3
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `project_name` | Name for output directories (`{name}_results/`, `{name}_logs/`) | — |
| `Input.dir` | Directory containing input files | `"data/"` |
| `Input.vcf` | VCF filename (`.vcf` or `.vcf.gz`) | — |
| `Input.metadata` | TSV with columns: site, sample, latitude, longitude [, traits...] | — |
| `Input.gff` | GFF3 annotation file | — |

CPU cores are auto-detected (`max(1, cpu_count - 2)`). Override with optional `cpu: N`.

---

### Filter — VCF quality filtering

```yaml
Filter:
  maf: 0.05
  snp_miss: 0.1
  sample_miss: 0.5
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Filter.maf` | Minor allele frequency threshold | `0.05` |
| `Filter.snp_miss` | Max SNP missingness (fraction) | `0.1` |
| `Filter.sample_miss` | Max sample missingness (fraction) | `0.5` |

---

### LD — Linkage disequilibrium pruning

```yaml
LD:
  window: 100
  step: 20
  r2: 0.2
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `LD.window` | Sliding window size (kb) | `100` |
| `LD.step` | Step size (kb) | `20` |
| `LD.r2` | R² threshold for pruning | `0.2` |

---

### sNMF — Population structure

```yaml
sNMF:
  k_start: 3
  k_end: 7
  k_best: 5
  repeats: 10
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sNMF.k_start` | Minimum K to test | `3` |
| `sNMF.k_end` | Maximum K to test | `7` |
| `sNMF.k_best` | Selected K (set after reviewing cross-entropy from `mode=prestructure`) | — |
| `sNMF.repeats` | Number of sNMF repetitions per K | `10` |

Ploidy is hardcoded to 2 (diploid). The pipeline does not support haploid organisms.

---

### Map — Geographic extent and climate resolution

```yaml
Map:
  climate_extent: auto
  gap: 0.5
  resolution: 0.5
  zoom_extent: "NULL"
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Map.climate_extent` | Geographic extent: `"auto"` (from sample coords) or `[min_lon, max_lon, min_lat, max_lat]` | `"auto"` |
| `Map.gap` | Buffer around samples (degrees) | `0.5` |
| `Map.resolution` | WorldClim resolution (arc-minutes: 0.5, 2.5, 5, 10) | `0.5` |
| `Map.zoom_extent` | Zoom coordinates: `"xmin,xmax,ymin,ymax"` or `"NULL"` | `"NULL"` |

**Regional zoom**: When `Map.zoom_extent` is set (e.g., `"34.8,35.8,32.0,33.3"`), the pipeline generates both full-extent and zoomed plots for all piemaps and genetic offset maps. Zoomed plots are organized in coordinate-named subdirectories under `zoom/`.

---

### Climate — Bioclimatic variables

```yaml
Climate:
  enabled: true
  predictors: bio_1,bio_12,bio_15
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Climate.enabled` | Enable climate download and climate-dependent analyses. Set `false` to skip GEA, maladaptation, and gea_x_gwas modes. Structure runs without climate plots/piemaps (only imputation + pop stats). GWAS runs without piemaps. | `true` |
| `Climate.predictors` | Comma-separated WorldClim bioclimatic variables for GEA/GWAS analyses (bio_1 through bio_19). All 19 are shown in the Structure tab piemaps; this list controls which are used as predictors in association. Ignored when `enabled: false`. | — |

---

### Population — Population statistics (optional)

```yaml
Population:
  calc_stats: false
  window_size: 10000
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Population.calc_stats` | Calculate Tajima's D, Pi diversity, IBD, Mantel test, AMOVA | `false` |
| `Population.window_size` | Window size (bp) for Tajima's D and Pi diversity | `10000` |

---

### LDdecay — LD decay analysis

```yaml
LDdecay:
  group_by: site
  min_samples: 5
  max_distance: 500
  scope: both
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `LDdecay.group_by` | Group samples by: `"site"` or `"cluster"` (sNMF K) | `"site"` |
| `LDdecay.min_samples` | Minimum samples per group to include | `5` |
| `LDdecay.max_distance` | Maximum pairwise SNP distance (kb) to compute | `500` |
| `LDdecay.scope` | Compute: `"genome_wide"`, `"per_chr"`, `"per_pop"`, or `"both"` (genome + pop) | `"both"` |

LD decay is used to determine `GEA.region_distance`: set `region_distance: auto` to use the r²<0.2 distance from LD decay results.

---

### Piemap — Pie map display

```yaml
Piemap:
  alpha: 0.8
  show_labels: true
  label_size: 3
  pie_scale: 0.5
  use_points: false
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Piemap.alpha` | Transparency of pie charts | `0.8` |
| `Piemap.show_labels` | Show site labels | `true` |
| `Piemap.label_size` | Label font size | `3` |
| `Piemap.pie_scale` | Base pie chart size multiplier | `0.5` |
| `Piemap.use_points` | Use points instead of pies (for very many sites) | `false` |

---

### GEA — Climate association analysis

```yaml
GEA:
  configs:
    - method: EMMAX
      adjust: bonf
      threshold: '0.05'
    - method: LFMM
      adjust: qvalue
      threshold: '0.05'
    - method: BLINK
      adjust: bonf
      threshold: '0.05'
  combine_gap: 200000
  region_distance: 2000000
  sig_snp_distance: 10000
  top_regions: 10
  promoter_length: 3000
  scattermore_threshold: 30000
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GEA.configs` | List of `{method, adjust, threshold}` entries (see below) | — |
| `GEA.combine_gap` | Gap (bp) for merging per-method SNP sets before region clustering | `200000` |
| `GEA.region_distance` | Distance (bp) for SNP clustering into regions (`"auto"` = from LD decay) | `2000000` |
| `GEA.sig_snp_distance` | Distance (bp) for merging overlapping significant SNPs | `10000` |
| `GEA.top_regions` | Maximum number of top regions to report per trait | `10` |
| `GEA.promoter_length` | Promoter length (bp) upstream of gene start for SNP counting | `3000` |
| `GEA.scattermore_threshold` | SNP count above which scattermore is used for Manhattan background rendering | `30000` |

**configs format**: Each entry specifies:
- `method`: `EMMAX`, `LFMM`, or any GAPIT3 model (`GLM`, `MLM`, `CMLM`, `ECMLM`, `SUPER`, `MLMM`, `FarmCPU`, `BLINK`)
- `adjust`: p-value correction — `bonf` (Bonferroni), `qvalue` (FDR), or a numeric threshold (top-N selection)
- `threshold`: significance threshold after correction

**Combine strategy** (set in Shiny, not config): `All` (union of all methods) or `MethodOverlap` (intersection).

---

### GFF — Feature parsing

```yaml
GFF:
  feature: mRNA
  gene_name: Name
  biotype: biotype
  go_field: ontology
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GFF.feature` | GFF feature type to extract genes (`gene`, `mRNA`, etc.) | `"gene"` |
| `GFF.gene_name` | Attribute for gene display name | `"Name"` |
| `GFF.biotype` | Attribute for biotype filtering | `"biotype"` |
| `GFF.go_field` | GFF attribute containing GO terms (or `"NULL"` to skip enrichment) | `"NULL"` |

---

### Enrichment — GO enrichment visualization

```yaml
Enrichment:
  top_terms: 20
  plot_width: 12
  plot_height: 10
  cnet_label: gene_id
  top_plot_regions: 5
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Enrichment.top_terms` | Max GO terms to show in plots | `20` |
| `Enrichment.plot_width` | Plot width (inches) | `12` |
| `Enrichment.plot_height` | Plot height (inches) | `10` |
| `Enrichment.cnet_label` | Gene label in cnetplot: `"gene_id"`, `"Name"`, or `"description"` | `"gene_id"` |
| `Enrichment.top_plot_regions` | Top regions per trait to generate plots for (`0` = all) | `5` |

GO enrichment is computed on-demand in the Shiny dashboard when the user clicks a region.

---

### Future — Climate projections

```yaml
Future:
  ssp: ssp370
  year: 2061-2080
  models:
    - ACCESS-CM2
    - CNRM-CM6-1
    - MIROC6
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Future.ssp` | Shared Socioeconomic Pathway scenario (`ssp126`, `ssp245`, `ssp370`, `ssp585`) | `"ssp370"` |
| `Future.year` | Future time period (`2021-2040`, `2041-2060`, `2061-2080`, `2081-2100`) | `"2061-2080"` |
| `Future.models` | List of CMIP6 models to ensemble-average | — |

---

### Maladaptation — Gradient Forest

```yaml
Maladaptation:
  methods:
    gradient_forest:
      ntree: 500
      cor_threshold: 0.5
      spatial_correction: both
      run_label: my_run
      random_model: false
      combine_gap: 200000
```

All settings are nested under `Maladaptation.methods.gradient_forest`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `ntree` | Number of trees in Gradient Forest | `500` |
| `cor_threshold` | Correlation threshold for variable selection | `0.5` |
| `spatial_correction` | Include PCNM spatial variables: `"with"`, `"without"`, or `"both"` | `"both"` |
| `run_label` | Label for output directory names | — |
| `random_model` | Also run a neutral (random SNP) GF model for comparison | `false` |
| `combine_gap` | Gap (bp) for merging SNPs from multiple GEA methods before GF | `200000` |

**Combine strategy** for adaptive SNPs (set in Shiny): `All` or `MethodOverlap` — controls which GEA sig SNPs are passed to Gradient Forest.

---

### GWAS — Phenotype association

```yaml
GWAS:
  configs:
    - method: EMMAX
      adjust: bonf
      threshold: '0.05'
  missing_strategy: DROP
  traits: height,flowering_time
  region_distance: 2000000
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GWAS.configs` | Same format as `GEA.configs` | — |
| `GWAS.missing_strategy` | How to handle missing values: `DROP`, `MEAN`, or `MEDIAN` | `"DROP"` |
| `GWAS.traits` | Comma-separated trait column names to include (omit to use all cols 5+) | all |
| `GWAS.region_distance` | Distance (bp) for SNP clustering (defaults to `GEA.region_distance`) | inherited |
| `GWAS.combine_gap` | Gap for merging SNP sets (defaults to `GEA.combine_gap`) | inherited |
| `GWAS.promoter_length` | Promoter length (defaults to `GEA.promoter_length`) | inherited |

**Missing value strategies**:
- **DROP** — Per-trait sample subsetting. Samples missing a trait value are excluded from that trait's analysis. Separate VCF, TPED, and kinship are computed per trait.
- **MEAN** — Missing values replaced with trait mean. All samples used together.
- **MEDIAN** — Missing values replaced with trait median. All samples used together.

---

### GEAxGWAS — GEA + GWAS overlap

```yaml
GEAxGWAS:
  pairwise:
    window_size: 500000
    min_snps: 1
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GEAxGWAS.pairwise.window_size` | Window (bp) around SNPs for window-based overlap scoring | `500000` |
| `GEAxGWAS.pairwise.min_snps` | Minimum overlapping SNPs to report a pair | `1` |

`region_distance` for overlap region computation is set interactively in the Shiny dashboard (defaults to `max(GEA.region_distance, GWAS.region_distance)`).

---

## GO Enrichment Analysis

When a GFF file with GO annotations is provided (`GFF.go_field` parameter), the pipeline can perform functional enrichment analysis on genes associated with significant regions.

### How it Works

Enrichment is computed **on-demand** in the Shiny dashboard when a region is selected:
1. **Gene finding** — Genes overlapping the region boundaries are identified from the GFF.
2. **Enrichment testing** — Per-region GO enrichment via clusterProfiler (Fisher's exact test).
3. **Visualization** — Three plot types generated: dotplot, emapplot, cnetplot.
4. **Persistence** — Results saved to `{MODULE}/plots/enrichment/` and `{MODULE}/tables/enrichment/` for instant reuse.

### GFF Requirements

Your GFF file must include GO terms in the attributes column. Specify the attribute via `GFF.go_field`:

```gff
chr1  AUGUSTUS  gene  1000  5000  .  +  .  ID=gene001;Name=GENE1;ontology=GO:0009414,GO:0042538
```

```yaml
GFF:
  go_field: ontology    # attribute containing GO terms
```

Set `GFF.go_field: "NULL"` to skip enrichment.

### Visualization Types

**Dotplot** — Enriched GO terms ranked by significance, with gene ratio and adjusted p-value. Works with ≥1 enriched term.

**Emapplot** — GO term similarity network. Nodes = GO terms sized by gene count. Requires ≥2 enriched terms.

**Cnetplot** — Gene-concept bipartite network. Gene labels configurable via `Enrichment.cnet_label` (`gene_id`, `Name`, or `description`). Requires ≥2 enriched terms.

---

## Output Structure

```
{project_name}_results/
├── Processing/
│   └── tables/                            # metadata.tsv, sample_missing_stats.tsv
│
├── PreStructure/
│   ├── plots/                             # pca.png, tracy_widom.png, cross_entropy_K{s}-{e}.png
│   │   └── K{k}/                          # structure_K{k}.png, pca_structure_K{k}.png
│   └── tables/
│       └── K{k}/                          # clusters_K{k}.tsv (Q-matrices)
│
├── climate/
│   ├── plots/                             # density_plot_present, correlation_heatmap, density_plot_future_*
│   ├── tables/present/                    # climate_present_{all,site,site_scaled}.tsv
│   ├── tables/future/                     # climate_future_year{Y}_ssp{S}_{site,all}.tsv
│   └── rasters/{present,future}/          # WorldClim .tif rasters
│
├── Structure/
│   ├── plots/piemap/                      # piemap_{bio}.png/svg + zoom/
│   ├── plots/piemap/{tajima_d,pi_diversity}/   # trait-scaled piemaps (optional)
│   ├── plots/pop_stats/                   # mantel_test, amova (optional)
│   ├── plots/ld_decay/                    # ld_decay_genome_wide, per_chr, per_pop
│   └── tables/pop_stats/                  # tajima_d_by_pop, pi_diversity_by_pop, ibd_*, amova
│
├── GEA/
│   ├── GAPIT_native_output/{model}/       # raw GAPIT output files
│   ├── plots/manhattan/
│   │   ├── {method}/                      # manhattan_{trait}_K{k}_{adjust}.png/svg
│   │   │                                  # qq_{trait}_K{k}_{adjust}.png/svg
│   │   └── combined/                      # manhattan_combined_K{k}.png/svg
│   ├── plots/enrichment/{trait}/          # region_{id}_dotplot/emapplot/cnetplot (Shiny on-demand)
│   └── tables/
│       ├── methods/{method}/              # {method}_pvalues_K{k}.tsv, _sig_snps_{adjust}.tsv
│       ├── selected_snps.tsv, regions_per_trait.tsv, regions_combined.tsv
│       ├── genes_per_region.tsv, genes_per_region_collapsed.tsv, genes_combined.tsv
│       └── enrichment/{trait}/            # GO enrichment TSVs (Shiny on-demand)
│
├── GWAS/                                  # same structure as GEA/
│   └── plots/piemap/                      # phenomap_{trait}.png/svg + zoom/
│
├── GEAxGWAS/
│   ├── plots/miami_combined_K{k}.{png,svg,_background.png,_coords.json}
│   └── tables/pairwise_{overlap_table,collapsed_snps}.tsv
│
├── Maladaptation/
│   ├── plots/gradient_forest/{run_label}_{spatial_tag}/
│   │   └── cumulative_importance, overall_importance, genetic_offset_piemap[_{tajima_d,pi_diversity}]
│   └── tables/gradient_forest/{run_label}_{spatial_tag}/
│       └── genetic_offset_{map,site}.tsv
│
├── pipeline_summary.tsv                   # all modes append here
│
├── _work/maf{}_miss{}_smiss{}/ld{}_win{}_step{}/    # parameterized intermediates
│   ├── filtered.vcf, pruned.vcf
│   ├── {name}.geno, .lfmm, .vcfsnp
│   └── {name}.snmfProject/
│
└── _intermediate/                         # internal: samples/, annotation/, {method}/, qs_cache/
    └── region_params.json                 # per-region Shiny parameter persistence

{project_name}_logs/{module}/              # per-module log directories
```

**Path parameterization**:
- `_work/maf0.05_miss0.1_smiss0.5/` — VCF filtered with these thresholds
- `_work/.../ld0.2_win100_step20/` — LD-pruned with R²=0.2, 100kb window, 20kb step
- `Maladaptation/.../my_run_spatial/` — GF run_label + spatial_correction tag
