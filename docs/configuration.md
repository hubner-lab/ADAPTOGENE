# Configuration Reference

All pipeline parameters are defined in `config.yaml`. This document covers every parameter group, GO enrichment details, and the output directory structure.

---

## Parameter Reference

### INPUT — Project setup

```yaml
INPUT_DIR: "data/"
INPUT_VCF: 'your_genotypes.vcf'
INPUT_METADATA: 'your_metadata.tsv'
INPUT_GFF: 'your_annotation.gff3'
PROJECT_NAME: 'MyProject'
CPU: 10
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `INPUT_DIR` | Directory containing input files | `"data/"` |
| `INPUT_VCF` | VCF filename (`.vcf` or `.vcf.gz`) | — |
| `INPUT_METADATA` | TSV with columns: site, sample, latitude, longitude [, traits...] | — |
| `INPUT_GFF` | GFF3 annotation file | — |
| `PROJECT_NAME` | Name for output directories (`{name}_results/`, `{name}_logs/`) | — |
| `CPU` | Max CPU cores for parallelizable steps | `10` |

### FILTER — VCF quality filtering

```yaml
FILTER_MAF: 0.05
FILTER_SNP_MISS: 0.1
FILTER_SAMPLE_MISS: 0.5
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `FILTER_MAF` | Minor allele frequency threshold | `0.05` |
| `FILTER_SNP_MISS` | Max SNP missingness (fraction) | `0.1` |
| `FILTER_SAMPLE_MISS` | Max sample missingness (fraction) | `0.5` |

### LD — Linkage disequilibrium pruning

```yaml
LD_WINDOW: 100
LD_STEP: 20
LD_R2: 0.2
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `LD_WINDOW` | Sliding window size (kb) | `100` |
| `LD_STEP` | Step size (kb) | `20` |
| `LD_R2` | R² threshold for pruning | `0.2` |

### SNMF — Population structure

```yaml
SNMF_K_START: 3
SNMF_K_END: 7
SNMF_K_BEST: 5
SNMF_PLOIDY: 2
SNMF_REPEATS: 10
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `SNMF_K_START` | Minimum K to test | `3` |
| `SNMF_K_END` | Maximum K to test | `7` |
| `SNMF_K_BEST` | Selected K (set after reviewing cross-entropy from `mode=structure`) | — |
| `SNMF_PLOIDY` | Ploidy level | `2` |
| `SNMF_REPEATS` | Number of sNMF repetitions per K | `10` |

### MAP — Geographic extent and climate resolution

```yaml
MAP_CROP_REGION: "auto"
MAP_GAP: 0.5
MAP_RESOLUTION: 0.5
MAP_REGIONMAP_EXTENT: "NULL"
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `MAP_CROP_REGION` | Geographic extent: `"auto"` (from sample coords) or `[min_lon, max_lon, min_lat, max_lat]` | `"auto"` |
| `MAP_GAP` | Buffer around samples (degrees) | `0.5` |
| `MAP_RESOLUTION` | WorldClim resolution (arc-minutes: 0.5, 2.5, 5, 10) | `0.5` |
| `MAP_REGIONMAP_EXTENT` | Zoom coordinates: `"xmin,xmax,ymin,ymax"` or `"NULL"` | `"NULL"` |

**Regional zoom**: When `MAP_REGIONMAP_EXTENT` is set to coordinates (e.g., `"34.8,35.8,32.0,33.3"`), the pipeline generates **both** full-extent and zoomed plots for all piemaps and genetic offset maps. Zoomed plots are organized in coordinate-named subdirectories (e.g., `regionmap/34p8_35p8_32p0_33p3/`).

### CLIMATE — Bioclimatic variables

```yaml
climate:
  enabled: true
  predictors: "bio_1,bio_12,bio_15"
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `climate.enabled` | Enable climate data download and climate-dependent analyses. Set `false` to skip GEA, maladaptation, and overlapping modes. Structure K runs without climate plots/piemaps (only imputation + pop stats). Phenotype association runs without piemaps. | `true` |
| `climate.predictors` | Comma-separated WorldClim bioclimatic variables (bio_1 through bio_19). Ignored when `enabled: false`. | — |

### POP — Population statistics (optional)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `POP_CALC_STATS` | Calculate population diversity statistics | `FALSE` |
| `POP_WINDOW_SIZE` | Window size for Tajima's D and Pi | `10000` |
| `POP_CUSTOM_TRAITS` | Custom traits for scaled piemaps | `"NULL"` |

### PIEMAP — Pie map display

| Parameter | Description | Default |
|-----------|-------------|---------|
| `PIEMAP_ALPHA` | Transparency of pie charts | `0.8` |
| `PIEMAP_LABELS` | Show site labels | `TRUE` |
| `PIEMAP_PIE_SCALE` | Base pie chart size | `0.5` |

### ASSOC — Association analysis

```yaml
ASSOC_CONFIGS:
  - method: "EMMAX"
    traits: "climate"
    correction: "bonferroni"
  - method: "LFMM"
    traits: "climate"
    correction: "qvalue"

ASSOC_COMBINE_METHOD: "Sum"
ASSOC_REGION_DISTANCE: 2000000
ASSOC_TOP_REGIONS: 10
ASSOC_GO_FIELD: "ontology"
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `ASSOC_CONFIGS` | List of `{method, traits, correction}` entries (see below) | — |
| `ASSOC_COMBINE_METHOD` | How to merge methods: `"Sum"`, `"Overlap"`, `"EMMAX"`, `"LFMM"` | `"Sum"` |
| `ASSOC_REGION_DISTANCE` | Distance (bp) for SNP clustering into regions | `2000000` |
| `ASSOC_TOP_REGIONS` | Number of top regions to report | `10` |
| `ASSOC_GO_FIELD` | GFF attribute containing GO terms (or `"NULL"` to skip enrichment) | `"NULL"` |

**ASSOC_CONFIGS format**: Each entry specifies a method (`EMMAX` or `LFMM`), trait type, and p-value correction method (`bonferroni`, `qvalue`, or a numeric threshold for top-N selection).

### GFF — Feature parsing

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GFF_FEATURE` | GFF feature type to extract genes | `"gene"` |
| `GFF_GENE_NAME` | Attribute for gene name | `"Name"` |
| `GFF_BIOTYPE` | Attribute for biotype filtering | `"biotype"` |

### ENRICHMENT — GO enrichment visualization

```yaml
ENRICHMENT_TOP_TERMS: 20
ENRICHMENT_PLOT_WIDTH: 12
ENRICHMENT_PLOT_HEIGHT: 10
ENRICHMENT_CNET_LABEL: "gene_id"
ENRICHMENT_TOP_PLOT_REGIONS: 5
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `ENRICHMENT_TOP_TERMS` | Max GO terms to show in plots | `20` |
| `ENRICHMENT_PLOT_WIDTH` | Plot width (inches) | `12` |
| `ENRICHMENT_PLOT_HEIGHT` | Plot height (inches) | `10` |
| `ENRICHMENT_CNET_LABEL` | Gene label in cnetplot: `"gene_id"`, `"Name"`, or `"description"` | `"gene_id"` |
| `ENRICHMENT_TOP_PLOT_REGIONS` | Top regions per trait to generate plots for (`0` = all) | `5` |

### REGIONPLOT — Regional visualization

| Parameter | Description | Default |
|-----------|-------------|---------|
| `REGIONPLOT_*` | Parameters for regional Manhattan plots with gene annotations (see `config.yaml`) | — |

### FUTURE — Climate projections

```yaml
FUTURE_SSP: "ssp370"
FUTURE_YEAR: "2061-2080"
FUTURE_MODELS: ["ACCESS-CM2", "CNRM-CM6-1", "MIROC6"]
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `FUTURE_SSP` | Shared Socioeconomic Pathway scenario | `"ssp370"` |
| `FUTURE_YEAR` | Future time period | `"2061-2080"` |
| `FUTURE_MODELS` | List of CMIP6 models to average | — |

### GF — Gradient Forest

| Parameter | Description | Default |
|-----------|-------------|---------|
| `GF_NTREE` | Number of trees in Gradient Forest | `500` |
| `GF_COR_THRESHOLD` | Correlation threshold for variable selection | `0.5` |
| `GF_PCNM` | Include PCNM spatial variables: `"with"`, `"without"`, or `"both"` | `"both"` |
| `GF_SUFFIX` | Suffix for GF output filenames | — |

### PHENO_ASSOC — Phenotype association (GWAS)

```yaml
PHENO_ASSOC_CONFIGS:
  - method: "EMMAX"
    traits: "phenotype"
    correction: "bonferroni"

PHENO_MISSING_STRATEGY: "DROP"
PHENO_REGION_DISTANCE: 2000000
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `PHENO_ASSOC_CONFIGS` | Association config for phenotype GWAS (same format as `ASSOC_CONFIGS`) | — |
| `PHENO_MISSING_STRATEGY` | How to handle missing phenotype values: `"DROP"`, `"MEAN"`, or `"MEDIAN"` | `"DROP"` |
| `PHENO_REGION_DISTANCE` | Distance (bp) for phenotype SNP clustering | `2000000` |

**Missing value strategies**:
- **DROP** — Per-trait sample subsetting. Samples with missing values for a trait are excluded from that trait's analysis. Separate VCF subset, TPED, and kinship are computed per trait.
- **MEAN** — Missing values imputed with trait mean. All samples included in a single analysis.
- **MEDIAN** — Missing values imputed with trait median. All samples included in a single analysis.

### OVERLAP — GEA + GWAS overlap

| Parameter | Description | Default |
|-----------|-------------|---------|
| `OVERLAP_REGION_DISTANCE` | Distance (bp) for combined region clustering | max(`ASSOC_REGION_DISTANCE`, `PHENO_REGION_DISTANCE`) |
| `OVERLAP_TOP_REGIONS` | Number of top overlap regions | — |
| `OVERLAP_PROMOTER_LENGTH` | Promoter length for SNP counting | — |

---

## GO Enrichment Analysis

When a GFF file with GO annotations is provided (`ASSOC_GO_FIELD` parameter), the pipeline automatically performs functional enrichment analysis on genes associated with significant regions.

### How it Works

1. **Region definition** — Significant SNPs are clustered into genomic regions via single-linkage clustering; boundaries are extended by `REGION_DISTANCE` during creation.
2. **Per-trait regions** — SNPs are clustered separately for each climate variable (e.g., bio_1, bio_12).
3. **Gene finding** — Genes overlapping region boundaries are identified.
4. **Enrichment testing** — Per-region GO enrichment via clusterProfiler.
5. **Visualization** — Three plot types generated for each region (see below).

### GFF Requirements

Your GFF file must include GO terms in the attributes column. Specify which attribute contains GO terms via `ASSOC_GO_FIELD`:

```gff
chr1  AUGUSTUS  gene  1000  5000  .  +  .  ID=gene001;Name=GENE1;ontology=GO:0009414,GO:0042538
```

```yaml
ASSOC_GO_FIELD: "ontology"    # Field name containing GO terms
```

Set to `"NULL"` to skip enrichment analysis.

### Visualization Types

**Dotplot** — Standard GO enrichment visualization. Shows enriched GO terms ranked by significance, with gene ratio, color-coded adjusted p-value, and gene count size. Always created for any region with enrichment.

**Emapplot** — GO term similarity network. Nodes = GO terms sized by gene count, edges = term similarity. Reveals clusters of related biological processes. Requires ≥2 enriched terms.

**Cnetplot** — Gene-concept network. Bipartite network of genes and their GO terms. Identifies hub genes appearing in multiple terms. Gene labels configurable via `ENRICHMENT_CNET_LABEL` (`gene_id`, `Name`, or `description`). Requires ≥2 enriched terms.

All plots saved in both PNG (300 dpi) and SVG (vector) formats.

### Enrichment Output Organization

```
{PROJECT}_results/
├── intermediate/enrichment/{trait}/
│   └── Region_{region_id}_enrichment.qs      # enrichResult objects
├── tables/association/enrichment/{trait}/
│   ├── Region_{region_id}_enrichment.tsv     # Per-region tables
│   └── Enrichment_{trait}_summary.tsv        # Per-trait summaries
└── plots/association/enrichment/{trait}/
    ├── Region_{region_id}_dotplot.{png,svg}
    ├── Region_{region_id}_emapplot.{png,svg}
    └── Region_{region_id}_cnetplot.{png,svg}
```

### Interpretation

- **Per-trait results** — Understand biological processes specific to each climate variable.
- **Regional specificity** — Different regions may show different functional enrichments.
- **Hub genes** (from cnetplot) — Prioritize genes appearing in multiple GO terms.
- **Term clusters** (from emapplot) — Identify overarching biological themes.

---

## Output Structure

```
{PROJECT_NAME}_results/
├── work/                                    # Intermediate files
│   └── maf{MAF}_miss{MISS}_smiss{SAMPLE_MISS}/
│       ├── {VCF_BASE}.vcf                  # Filtered VCF
│       └── ld{r2}_win{WIN}_step{STEP}/
│           ├── {VCF_BASE}.vcf              # LD-pruned VCF
│           ├── {VCF_BASE}.geno             # LEA geno format
│           ├── {VCF_BASE}.lfmm             # LEA lfmm format
│           ├── {VCF_BASE}.snmfProject/     # sNMF results
│           └── {VCF_BASE}_K{K_BEST}imp.lfmm  # Imputed genotypes
│
├── plots/                                   # All visualizations
│   ├── pca/                                # PCA plots
│   ├── structure/                          # Structure plots, cross-entropy
│   ├── climate_{RESOLUTION}/               # Climate maps
│   ├── piemap/                             # Geographic ancestry maps
│   │   └── regionmap/{coords}/             # Zoomed piemaps
│   ├── EMMAX/                              # EMMAX Manhattan plots
│   ├── LFMM/                               # LFMM Manhattan plots
│   ├── association/                        # Combined Manhattan plots
│   │   ├── Manhattan_all_traits_combined_K{K}.[png,svg]
│   │   ├── Manhattan_all_traits_combined_K{K}_regions.[png,svg]
│   │   └── enrichment/{trait}/             # GO enrichment plots
│   ├── association_phenotypes/             # Phenotype Manhattan + piemaps
│   ├── regionplot/                         # Regional plots with genes
│   └── gradientForest/                     # GF importance, offset maps
│       └── regionmap/{coords}/             # Zoomed GF plots
│
├── tables/                                  # TSV outputs
│   ├── structure/                          # Q-matrices, cluster assignments
│   ├── climate_{RESOLUTION}/               # Climate values per sample
│   ├── association/
│   │   ├── {METHOD}/                       # Per-method significant SNPs
│   │   ├── Selected_SNPs.tsv               # Combined significant SNPs
│   │   ├── Regions_per_trait.tsv           # Regions per climate variable
│   │   ├── Regions_climate_combined.tsv    # Combined climate regions
│   │   ├── Genes_per_region.tsv            # Genes in per-trait regions
│   │   ├── Genes_per_region_climate_combined.tsv
│   │   └── enrichment/{trait}/             # GO enrichment results
│   ├── association_phenotypes/             # Phenotype GWAS results
│   ├── overlapping/                        # GEA+GWAS overlap results
│   ├── Pipeline_summary.tsv                # Running summary across all modes
│   └── gradientForest/                     # GF predictions, offset scores
│
└── intermediate/                            # Shared files
    ├── samples.list
    ├── Climate_present_RasterStack.grd     # Current climate rasters
    ├── Climate_future_year{YEAR}_ssp{SSP}_RasterStack.grd
    └── enrichment/{trait}/                 # enrichResult objects (.qs)

{PROJECT_NAME}_logs/                         # Snakemake rule logs
    └── {rule_name}.log
```

**Path organization**: Files are organized by the parameters that created them:
- `maf0.05_miss0.1_smiss0.5/` — VCF filtered with these thresholds
- `ld0.2_win100_step20/` — LD-pruned with R²=0.2, 100kb window, 20kb step
