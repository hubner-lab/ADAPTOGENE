# ADAPTOGENE
**A standardized population genomics pipeline for association mapping and adaptive potential assessment**

---

## Overview

**ADAPTOGENE** is a Dockerized, Snakemake-based bioinformatics pipeline for **population genomics**, designed to standardize and automate complex analyses within a single, reproducible framework. The pipeline integrates three tightly connected analytical layers:

1. **VCF preprocessing and population structure analysis**
   Comprehensive processing of raw genotype data, including:
   - Variant filtering and LD pruning
   - Principal Component Analysis (PCA)
   - Model-based ancestry inference (sNMF) across a range of K values

   This step provides a robust description of population structure, supports informed selection of the optimal number of genetic clusters, and generates the foundation for downstream association and landscape genomics analyses.

2. **GWAS / GEA analyses**
   Identification of loci and candidate genes associated with:
   - Climate variables (bioclimatic predictors)
   - Phenotypic traits

   ADAPTOGENE simplifies comparison across traits and methods and explicitly detects **overlapping association signals**, enabling systematic identification of shared genetic drivers.

3. **Assessment of population adaptive potential**
   Landscape genomics–based evaluation of maladaptation risk under present and future climate scenarios, supporting **population prioritization for conservation and management**.

Overall, ADAPTOGENE standardizes complex population genomic workflows, automates data handling (including climate data acquisition), and ensures full reproducibility through containerization.

---

## Key Features

- Unified framework for **population structure, association, and maladaptation**
- Automated **climate data download and preprocessing** (WorldClim bioclimatic variables)
- Explicit detection of **overlapping association signals** across traits and methods
- **GO enrichment analysis** with comprehensive visualizations:
  - Dotplot - Standard enrichment visualization
  - Emapplot - GO term similarity network showing clustering
  - Cnetplot - Gene-concept network with configurable gene labels (ID, Name, or description)
- **Per-region and per-trait** functional enrichment reporting
- **Cross-trait evidence** - Per-trait regions report overlapping signals from other traits
- **Dual region outputs** - Per-trait regions for specific analysis + combined climate regions with per-trait and per-method SNP breakdowns
- **Exon/promoter SNP counting** - Genes annotated with counts of significant SNPs in exons and promoters
- **Pipeline summary table** - Running summary accumulated across all modes
- **Combined Manhattan plots** - Visualize all traits and methods simultaneously
- **Interactive Shiny app** - Explore all results with interactive Manhattan plots, clickable regions, gene/enrichment tables, and Gradient Forest maps
- Reproducible execution via **Docker** (all dependencies included)
- Designed for **landscape and conservation genomics**
- Minimal manual intervention once configured

---

## Design Principles

### Compute Once, Reuse Downstream

Each pipeline stage produces a **finished entity** that downstream stages consume as-is. Intermediate results should never be recomputed or re-derived in later scripts.

**Regions are the key example:**

| Stage | What happens | Script |
|-------|-------------|--------|
| **Region creation** | SNPs clustered via single-linkage, boundaries extended by `REGION_DISTANCE` on each side | `create_regions.R` |
| **Gene finding** | Genes overlapping region boundaries (used as-is, no re-extension) | `find_genes_around_regions.R` |
| **Cross-trait evidence** | Other-trait SNPs falling within region boundaries (used as-is, no re-extension) | `create_regions.R`, `combine_overlap_snps.R` |
| **Enrichment** | GO analysis on genes already found per region | `run_enrichment.R` |
| **Visualization** | Region rectangles drawn at region start/end coordinates | `plot_manhattan_combined.R`, `app.R` |

`REGION_DISTANCE` is a parameter for `create_regions.R` only. No downstream script receives or uses it. The region table's `start` and `end` columns are the single source of truth for region boundaries.

### Structured Data Flow

Each script should have a clear, single responsibility:
- **Input**: Read finished outputs from upstream stages
- **Process**: Perform one well-defined operation
- **Output**: Write finished outputs for downstream stages

Scripts should not:
- Re-derive values that upstream scripts already computed (e.g., re-extending regions)
- Accept parameters they don't use (pass-through arguments)
- Contain dead code guarded by conditions that are always false (e.g., `if distance > 0` when always called with `0`)

### Normalize Early

Data normalization (chromosome names, sample IDs, coordinate systems) happens at the earliest possible step. All downstream scripts assume normalized inputs. See [Chromosome Naming](#important-chromosome-naming) for the primary example.

---

## Requirements

- Docker
- Unix-like operating system (Linux/macOS recommended)
- Sufficient disk space for climate data (~2GB for WorldClim)

All software dependencies are included inside the Docker image.

---

## Installation

```bash
git clone https://github.com/hubner-lab/ADAPTOGENE.git
cd ADAPTOGENE
docker build -t adaptogene .
```

---

## Input Data

All input files must be placed in the directory specified by `INPUT_DIR` in `config.yaml`.

### Required Inputs

- **VCF** — Genotype data (supports `.vcf` and `.vcf.gz`)
- **Sample metadata** — TSV with columns: `site`, `sample`, `latitude`, `longitude`
- **GFF3** — Genome annotation (required for association mode)

### Optional Inputs

- Phenotypic trait tables
- Custom population-level traits

### Important: Chromosome Naming

The pipeline normalizes chromosome names early by stripping "chr" prefix (chr1→1, chr2H→2H). This ensures:
- Compatibility with LEA package (used for sNMF and imputation)
- Consistency across all outputs (association tables, regions, gene lists)
- Matching between VCF and GFF annotations

If your input files use "chr" prefix, they will be automatically normalized during processing. All output files will use normalized names.

---

## Configuration

All pipeline parameters are defined in `config.yaml`. Only this file needs modification between runs.

### Key Parameter Groups

**INPUT** - Files and project settings
```yaml
INPUT_DIR: "data/"
INPUT_VCF: 'your_genotypes.vcf'
INPUT_METADATA: 'your_metadata.tsv'
INPUT_GFF: 'your_annotation.gff3'
PROJECT_NAME: 'MyProject'
CPU: 10
```

**FILTER** - VCF quality filtering
```yaml
FILTER_MAF: 0.05              # Minor allele frequency threshold
FILTER_SNP_MISS: 0.1          # Max SNP missingness
FILTER_SAMPLE_MISS: 0.5       # Max sample missingness
```

**LD** - Linkage disequilibrium pruning
```yaml
LD_WINDOW: 100                # Window size in kb
LD_STEP: 20                   # Step size for sliding window
LD_R2: 0.2                    # R² threshold
```

**SNMF** - Population structure
```yaml
SNMF_K_START: 3               # Minimum K to test
SNMF_K_END: 7                 # Maximum K to test
SNMF_K_BEST: 5                # Selected K (set after reviewing cross-entropy)
SNMF_PLOIDY: 2
SNMF_REPEATS: 10
```

**MAP** - Geographic extent for climate data
```yaml
MAP_CROP_REGION: "auto"       # or [min_lon, max_lon, min_lat, max_lat]
MAP_GAP: 0.5                  # Buffer around samples (degrees)
MAP_RESOLUTION: 0.5           # Climate resolution (0.5, 2.5, 5, 10 arc-min)
MAP_REGIONMAP_EXTENT: "NULL"  # Optional zoom: "xmin,xmax,ymin,ymax" or "NULL"
```

**Note on Regional Zoom**: When `MAP_REGIONMAP_EXTENT` is set to coordinates (e.g., `"34.8,35.8,32.0,33.3"`), the pipeline automatically generates **both** full extent and zoomed plots for all piemaps and genetic offset maps. Zoomed plots are organized in coordinate-named subdirectories (e.g., `regionmap/34p8_35p8_32p0_33p3/PieMap_bio_1.png`). This keeps filenames short and groups all zoomed plots by extent.

**CLIMATE** - Bioclimatic variables
```yaml
CLIMATE_PREDICTORS: "bio_1,bio_12,bio_15"  # Select from 19 WorldClim variables
```

**ASSOC** - Association analysis
```yaml
ASSOC_CONFIGS:
  - method: "EMMAX"
    traits: "climate"
    correction: "bonferroni"
  - method: "LFMM"
    traits: "climate"
    correction: "qvalue"

ASSOC_COMBINE_METHOD: "Sum"   # How to merge methods: "Sum", "Overlap", "EMMAX", "LFMM"
ASSOC_REGION_DISTANCE: 2000000  # Distance for SNP clustering and gene finding (bp)
ASSOC_TOP_REGIONS: 10         # Number of top regions to report
ASSOC_GO_FIELD: "ontology"    # GFF field containing GO terms (or "NULL")
```

**ENRICHMENT** - GO enrichment visualization
```yaml
ENRICHMENT_TOP_TERMS: 20      # Max GO terms to show in plots
ENRICHMENT_PLOT_WIDTH: 12     # Plot width in inches
ENRICHMENT_PLOT_HEIGHT: 10    # Plot height in inches
ENRICHMENT_CNET_LABEL: "gene_id"  # Gene label in cnetplot: "gene_id", "Name", "description"
ENRICHMENT_TOP_PLOT_REGIONS: 5    # Top regions per trait to generate plots for (0 = all)
```

**FUTURE** - Climate projections for maladaptation
```yaml
FUTURE_SSP: "ssp370"          # Climate scenario
FUTURE_YEAR: "2061-2080"      # Time period
FUTURE_MODELS: ["ACCESS-CM2", "CNRM-CM6-1", "MIROC6"]  # CMIP6 models
```

---

## Running ADAPTOGENE

### Method 1: One-line execution (recommended for production)

```bash
docker run --user $(id -u):$(id -g) --rm --memory=20g -v $PWD:/pipeline adaptogene:latest \
  snakemake -c4 -s Snakefile --config mode=<MODE> --configfile config.yaml --scheduler greedy
```

Replace `<MODE>` with your desired pipeline mode (see below).

### Method 2: Interactive session (recommended for development/debugging)

```bash
# Enter Docker container
docker run -w /pipeline --user $(id -u):$(id -g) -it -v $PWD:/pipeline adaptogene bash

# Inside container, run Snakemake
snakemake -s Snakefile --configfile config.yaml --config mode=<MODE> --cores 24 --scheduler greedy
```

### Snakemake Options

- `-c <N>` or `--cores <N>` - Number of CPU cores
- `-n` - Dry run (show what would be executed)
- `-R <rule>` - Force re-run specific rule and downstream dependencies
- `--forcerun <rule>` - Force re-run only specific rule
- `-F` - Force re-run all rules
- `--scheduler greedy` - Recommended for optimal job scheduling

---

## Pipeline Modes

### 1. Processing (`mode=processing`)

**Purpose**: Initial VCF filtering and LD pruning

**Input**
- Raw VCF file
- Sample metadata

**Operations**
- Sample missingness filtering (removes samples with >FILTER_SAMPLE_MISS missing)
- SNP quality filtering (MAF, missingness thresholds)
- LD pruning (sliding window approach)
- Format conversion to LEA formats (GENO, LFMM)
- Chromosome name normalization (strips "chr" prefix)

**Output**
- Filtered and LD-pruned VCF
- LEA format files (`.geno`, `.lfmm`)
- Sample list

---

### 2. Population Structure (`mode=structure`)

**Purpose**: Population structure inference and PCA

**Input**
- Filtered, LD-pruned VCF
- Sample metadata

**Operations**
- Principal Component Analysis (PCA)
- sNMF ancestry inference across K range (SNMF_K_START to SNMF_K_END)
- Cross-entropy calculation for K selection

**Output**
- PCA plots (colored by site/population)
- Cross-entropy plots for K selection
- Ancestry proportion matrices (Q-matrices) for each K
- Structure plots showing ancestry proportions

**User Action Required**
1. Inspect cross-entropy plots
2. Select optimal K based on cross-entropy minimum
3. Set `SNMF_K_BEST` in `config.yaml`

---

### 3. Spatial Structure Analysis (`mode=structure_K`)

**Purpose**: Visualize ancestry in geographic space

**Requirements**
- `SNMF_K_BEST` must be set in config

**Input**
- Selected K ancestry proportions
- Geographic coordinates
- Climate data (automatically downloaded)

**Operations**
- Imputation of missing genotypes using sNMF Q-matrix
- Climate data download (WorldClim current conditions)
- Climate value extraction at sample locations
- Optional: Population diversity statistics (π, Tajima's D) if `POP_CALC_STATS=TRUE`
- Geographic visualization of ancestry proportions (piemaps)

**Output**
- Geographic ancestry maps (piemaps) for each climate variable
- Climate value maps
- Population diversity metrics (if calculated)
- Cluster assignments

---

### 4. Association Analysis (`mode=association`)

**Purpose**: Identify SNPs and genes associated with traits/climate

**Requirements**
- `SNMF_K_BEST` must be set
- `INPUT_GFF` must be provided
- Climate data downloaded (from structure_K mode)

**Input**
- Full (non-LD-pruned) VCF for EMMAX
- LD-pruned VCF for LFMM training
- Climate variables or phenotypic traits
- Genome annotation (GFF3)

**Methods**
- **EMMAX**: Mixed model accounting for population structure via kinship matrix
- **LFMM**: Latent factor mixed model for genotype-environment associations

**Operations**
1. Imputation of full dataset using sNMF Q-matrix
2. Association testing (EMMAX and/or LFMM)
3. Multiple testing correction (Bonferroni, q-value, or top-N)
4. Combining significant SNPs across methods (Sum/Overlap/single method)
5. **Region clustering**:
   - Per-trait regions (SNPs clustered separately for each climate variable), with cross-trait evidence
   - Combined climate regions (all significant SNPs clustered together), with per-trait and per-method SNP breakdowns
6. **Gene annotation**: Find genes overlapping regions (regions already have extended boundaries from clustering), with exon and promoter SNP counts
7. **GO enrichment analysis** (if ASSOC_GO_FIELD is set):
   - Per-region enrichment testing
   - Organized by trait
   - Configurable number of top regions to plot per trait (`ENRICHMENT_TOP_PLOT_REGIONS`)
   - Three visualization types per region, with configurable gene labels in cnetplot

**Output**

**Tables** (`tables/association/`):
- `Selected_SNPs.tsv` - All significant SNPs (combined across methods)
- `Regions_per_trait.tsv` - Genomic regions per climate trait, with cross-trait evidence (`other_traits`, `other_snp_count`)
- `Regions_climate_combined.tsv` - Combined climate-associated regions, with per-trait (e.g., `bio_1_snps`) and per-method (e.g., `EMMAX_snps`) SNP counts
- `Genes_per_region.tsv` - Genes within per-trait regions, with `exon_snp_count` and `promoter_snp_count`
- `Genes_per_region_climate_combined.tsv` - Genes within combined regions
- `enrichment/{trait}/Region_{region_id}_enrichment.tsv` - Per-region GO results
- `enrichment/{trait}/Enrichment_{trait}_summary.tsv` - Per-trait summaries

**Pipeline Summary** (`tables/Pipeline_summary.tsv`):
- Running summary accumulated across all pipeline modes
- Long format: `step`, `metric`, `value`
- Each mode appends its metrics (processing: sample/SNP counts; structure: K range; association: region/gene counts; maladaptation: offset statistics)

**Plots** (`plots/`):
- `{METHOD}/Manhattan_{trait}_K{K}.png` - Per-trait Manhattan plots
- `{METHOD}/Manhattan_{trait}_K{K}_regions.png` - With regions highlighted
- `association/Manhattan_all_traits_combined_K{K}.png` - All traits and methods
- `association/Manhattan_all_traits_combined_K{K}_regions.png` - With combined regions
- `association/enrichment/{trait}/Region_{region_id}_{dotplot,emapplot,cnetplot}.png` - GO visualizations

**Manhattan Plot Features**:
- Per-trait plots show associations for single climate variable
- Combined plots show ALL traits and methods simultaneously
- Color-coded by trait-method combination
- Multi-association SNPs highlighted (significant in multiple trait-method combos)
- Region rectangles for visual identification of clustered SNPs

---

### 5. Phenotype Association (`mode=association_phenotypes`)

**Purpose**: GWAS on phenotypic traits stored in metadata columns 5+ (after site, sample, latitude, longitude)

**Requirements**
- `SNMF_K_BEST` must be set
- Phenotype columns present in metadata (numeric, columns 5+)

**Missing Value Strategies**
- **DROP**: Per-trait sample subsetting (separate PCA/kinship per trait)
- **MEAN/MEDIAN**: Impute missing values, single sample set for all traits

**Operations**
1. Extract and validate phenotype traits from metadata
2. Handle missing values per configured strategy
3. Per-trait VCF subsetting and PCA (DROP mode)
4. EMMAX association testing per trait
5. Combine per-trait results
6. Region clustering, gene annotation, GO enrichment (reuses existing scripts)
7. Per-trait Manhattan plots and geographic piemaps

**Output**
- `tables/association_phenotypes/` — p-values, sig SNPs, regions, genes, enrichment
- `plots/association_phenotypes/` — Manhattan plots, piemaps, enrichment visualizations
- `Pipeline_summary.tsv` — per-trait missing value info, dropped samples

**Config**: Independent `PHENO_*` parameters (`PHENO_ASSOC_CONFIGS`, `PHENO_MISSING_STRATEGY`, `PHENO_REGION_DISTANCE`, etc.)

---

### 6. Overlapping Regions (`mode=overlapping`)

**Purpose**: Combine GEA (climate association) and GWAS (phenotype association) results to identify shared genomic signals

**Requirements**
- Completed `mode=association` (Selected_SNPs.tsv, Regions_climate_combined.tsv)
- Completed `mode=association_phenotypes` (Selected_SNPs.tsv, Regions_phenotype_combined.tsv)

**Operations**
1. Merge all significant SNPs from both GEA and GWAS into a unified set
2. Compute overlap between GEA and GWAS regions (regions already have extended boundaries from clustering)
3. Create NEW combined regions from ALL climate + phenotype SNPs (same single-linkage clustering)
4. Find genes around new combined regions
5. Run GO enrichment on new combined regions
6. Output overlap statistics (overlap pairs, percentages, shared SNPs, contributing factors)

**Output**
- `tables/overlapping/Selected_SNPs_all.tsv` — Union of all significant SNPs
- `tables/overlapping/Regions_per_trait_all.tsv` — Per-trait regions with `source` column (climate/phenotype)
- `tables/overlapping/Regions_all_combined.tsv` — All SNPs clustered together with `source`, `climate_snps`, `phenotype_snps` columns
- `tables/overlapping/Overlap_summary.tsv` — Pairwise overlap between GEA and GWAS extended regions
- `tables/overlapping/Genes_per_region.tsv` — Genes within new combined regions
- `tables/overlapping/enrichment/` — Per-trait GO enrichment results

**Config**: Optional `OVERLAP_REGION_DISTANCE`, `OVERLAP_TOP_REGIONS`, `OVERLAP_PROMOTER_LENGTH` (defaults to max of ASSOC/PHENO values)

---

### 7. Regional Visualization (`mode=regionplot`)

**Purpose**: Detailed visualization of specific genomic regions

**Input**
- Association results
- GFF annotation
- Top regions from association analysis

**Output**
- Regional Manhattan plots with gene annotations
- Publication-ready figures for candidate regions

---

### 8. Adaptive Potential and Maladaptation (`mode=maladaptation`)

**Purpose**: Model genotype-environment relationships and predict maladaptation risk

**Requirements**
- Completed association analysis (`Selected_SNPs.tsv`)
- Climate data (present and future)

**Input**
- Adaptive SNP set from association mode
- Present climate layers (from structure_K mode)
- Future climate layers (automatically downloaded)
- Population coordinates

**Methods**
- **Gradient Forest**: Non-parametric ensemble model
- Climate scenarios from CMIP6 (multiple models averaged)
- Optional PCNM spatial variables (`GF_PCNM: 'with'` or `'without'`, producing `_PCNM` or `_noPCNM` suffix in output files)

**Operations**
1. Download future climate projections for selected SSP/year/models
2. Average across climate models (ensemble projection)
3. Train Gradient Forest on adaptive SNPs and current climate
4. Predict genomic offset under future climate
5. Calculate maladaptation scores per population

**Output**
- Adaptive potential maps
- Genetic offset maps (difference between current and future predicted genotypes)
- Population-level maladaptation scores
- Variable importance plots
- Cumulative importance plots
- Population prioritization tables for conservation
- Zoomed regional maps in coordinate subdirectories (when `MAP_REGIONMAP_EXTENT` is set)

---

## GO Enrichment Analysis

When a GFF file with GO annotations is provided (`ASSOC_GO_FIELD` parameter), the pipeline automatically performs functional enrichment analysis on genes associated with significant regions.

### How it Works

1. **Region Definition**: Significant SNPs are clustered into genomic regions using `ASSOC_REGION_DISTANCE` (regions are extended by this distance during creation)
2. **Per-Trait Regions**: SNPs are clustered separately for each climate variable (e.g., bio_1, bio_12)
3. **Gene Finding**: Genes overlapping regions are identified (regions already have extended boundaries from clustering)
4. **Enrichment Testing**: Per-region GO enrichment using clusterProfiler
5. **Visualization**: Three plot types generated for each region

### GFF Requirements

Your GFF file must include GO terms in the attributes column. Specify which attribute contains GO terms via `ASSOC_GO_FIELD`:

```gff
chr1  AUGUSTUS  gene  1000  5000  .  +  .  ID=gene001;Name=GENE1;ontology=GO:0009414,GO:0042538
```

In `config.yaml`:
```yaml
ASSOC_GO_FIELD: "ontology"    # Field name containing GO terms
```

Set to `"NULL"` to skip enrichment analysis.

### Output Organization

```
results/
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

### Visualization Types

**Dotplot** - Standard GO enrichment visualization
- Shows enriched GO terms ranked by significance
- Gene ratio (fraction of genes in term)
- Color-coded by adjusted p-value
- Size represents gene count
- Always created for any region with enrichment

**Emapplot** - GO term similarity network
- Network showing GO term similarities and clustering
- Node size = number of genes in term
- Edge thickness = term similarity
- Identifies clusters of related biological processes
- Useful for understanding functional themes
- Created when ≥2 terms are enriched

**Cnetplot** - Gene-concept network
- Bipartite network showing genes and their GO terms
- Identifies hub genes (genes in multiple terms)
- Shows which specific genes drive enrichment
- Gene labels configurable via `ENRICHMENT_CNET_LABEL`: gene ID, Name, or description
- Useful for candidate gene prioritization
- Created when ≥2 terms are enriched

All plots saved in both PNG (300 dpi) and SVG (vector) formats for publication.

### Configuration Parameters

```yaml
ENRICHMENT_TOP_TERMS: 20        # Maximum GO terms to show in plots
ENRICHMENT_PLOT_WIDTH: 12       # Plot width in inches
ENRICHMENT_PLOT_HEIGHT: 10      # Plot height in inches
ENRICHMENT_CNET_LABEL: "gene_id"  # Gene label in cnetplot: "gene_id", "Name", "description"
ENRICHMENT_TOP_PLOT_REGIONS: 5    # Top regions per trait to generate plots for (0 = all)
```

### Interpretation

- **Per-trait results**: Understand biological processes specific to each climate variable
- **Regional specificity**: Different regions may show different functional enrichments
- **Hub genes** (from cnetplot): Prioritize genes appearing in multiple GO terms
- **Term clusters** (from emapplot): Identify overarching biological themes

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
│   │   └── regionmap/{coords}/             # Zoomed piemaps in coordinate subdirs
│   ├── EMMAX/                              # EMMAX Manhattan plots
│   ├── LFMM/                               # LFMM Manhattan plots
│   ├── association/                        # Combined Manhattan plots
│   │   ├── Manhattan_all_traits_combined_K{K}.[png,svg]
│   │   ├── Manhattan_all_traits_combined_K{K}_regions.[png,svg]
│   │   └── enrichment/{trait}/             # GO enrichment plots
│   ├── regionplot/                         # Regional plots with genes
│   └── gradientForest/                     # GF importance, offset maps
│       └── regionmap/{coords}/             # Zoomed plots in coordinate subdirs
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
│   │   ├── Genes_per_region_climate_combined.tsv  # Genes in combined regions
│   │   └── enrichment/{trait}/             # GO enrichment results
│   │       ├── Region_{region_id}_enrichment.tsv
│   │       └── Enrichment_{trait}_summary.tsv
│   ├── Pipeline_summary.tsv                # Running summary across all modes
│   └── gradientForest/                     # GF predictions, offset scores
│
└── intermediate/                            # Shared files
    ├── samples.list
    ├── Climate_present_RasterStack.grd     # Current climate rasters
    ├── Climate_future_year{YEAR}_ssp{SSP}_RasterStack.grd  # Future climate
    └── enrichment/{trait}/                 # enrichResult objects (.qs)

{PROJECT_NAME}_logs/                         # Snakemake rule logs
    └── {rule_name}.log                     # One log per rule
```

**Path organization**: Files are organized by the parameters that created them:
- `maf0.05_miss0.1_smiss0.5/` - VCF filtered with these thresholds
- `ld0.2_win100_step20/` - LD-pruned with R²=0.2, 100kb window, 20kb step

---

## Workflow Summary

**Typical workflow for climate adaptation study**:

1. **Processing**: Filter VCF, LD prune
   ```bash
   mode=processing
   ```

2. **Structure**: Run sNMF, inspect cross-entropy, select K_BEST
   ```bash
   mode=structure
   # Review cross-entropy plot, set SNMF_K_BEST in config
   ```

3. **Spatial Structure**: Visualize ancestry, download climate
   ```bash
   mode=structure_K
   ```

4. **Association**: Identify adaptive SNPs and genes
   ```bash
   mode=association
   # Produces regions, genes, GO enrichment
   ```

5. **Regional Plots** (optional): Detailed region visualization
   ```bash
   mode=regionplot
   ```

6. **Phenotype Association** (optional): GWAS on phenotypic traits
   ```bash
   mode=association_phenotypes
   # Requires phenotype columns (5+) in metadata file
   ```

7. **Maladaptation**: Predict future vulnerability
   ```bash
   mode=maladaptation
   ```

---

## Interactive Results Viewer (Shiny App)

ADAPTOGENE includes an interactive Shiny dashboard for exploring pipeline results. The app auto-discovers all projects in the pipeline directory and organizes outputs into seven analysis tabs.

### Running the App

```bash
docker run --user $(id -u):$(id -g) --rm -p 3838:3838 -v $PWD:/pipeline adaptogene:latest \
  Rscript -e "shiny::runApp('/pipeline/scripts/app.R', host='0.0.0.0', port=3838)"
```

Then open http://localhost:3838 in your browser.

### Features

**Structure Tab** — PCA, Tracy-Widom test, cross-entropy plot, and K-dependent visualizations (ancestry piemaps, population differentiation, SNMF structure barplots) with a K slider.

**Structure K Tab** — Climate correlation heatmap, bio-variable density plots and scaled piemaps with dropdown selectors for bio variable and pie-scaling metric.

**Association Tab** — Interactive Manhattan plot combining all traits and methods:
- Color encodes **trait** (e.g., bio_1 = red, bio_12 = blue), shape encodes **method** (circles = EMMAX, triangles = LFMM)
- Significant SNPs are interactive with hover tooltips; non-significant SNPs are static for performance
- Region rectangles show regions (boundaries include `ASSOC_REGION_DISTANCE` extension from clustering)
- Click a significant SNP to select its region and filter the gene and GO enrichment tables below
- Summary boxes show unique sig SNPs, sig hits on plot (SNP × trait × method), regions, genes, GO terms, and methods
- Controls: toggle individual trait×method layers, switch region type (combined/per-trait), adjust top N regions

**Phenotype Association Tab** — Interactive GWAS results for phenotypic traits:
- Same Manhattan architecture as Association tab (color = trait, shape = method)
- Per-trait phenotype distribution histograms and geographic piemaps
- Gene and GO enrichment tables linked to region selection
- Handles missing data via DROP/MEAN/MEDIAN strategies

**Overlapping Regions Tab** — Miami Manhattan plot comparing GEA and GWAS signals:
- Requires `mode=overlapping` pipeline output (pre-computed overlap stats, new combined regions)
- Miami plot: GEA significant SNPs (triangles) above x-axis, GWAS significant SNPs (circles) below
- Green rectangles highlight new combined regions containing both climate and phenotype sources
- Summary valueBoxes: GEA regions, GWAS regions, overlap pairs, new combined, genes, GO terms
- Overlap pairs table with overlap length, percentage, SNP counts, and contributing traits
- Region selector with genes and GO enrichment tables (same interactive pattern as other tabs)

**Maladaptation Tab** — Gradient Forest model results:
- Model selector to switch between PCNM variants
- Overall and Cumulative Importance diagnostic plots
- Genetic Offset piemaps (base, with Tajima's D, with Pi Diversity)
- Auto-discovered regional zoom maps
- Interactive site-level genetic offset table with summary statistics

---

## Intended Applications

- **Landscape genomics** - Understanding spatial genetic variation
- **Climate adaptation studies** - Identifying adaptive loci and vulnerable populations
- **Conservation genomics** - Population prioritization for management
- **Multi-trait association analysis** - Detecting shared genetic architecture

---

## Citation

If you use ADAPTOGENE in published work, please cite:

- This repository: `https://github.com/hubner-lab/ADAPTOGENE`
- Underlying methods:
  - **sNMF**: Frichot et al. 2014, Molecular Biology and Evolution
  - **EMMAX**: Kang et al. 2010, Nature Genetics
  - **LFMM**: Frichot et al. 2013, Molecular Biology and Evolution
  - **Gradient Forest**: Ellis et al. 2012, Ecology
  - **WorldClim**: Fick & Hijmans 2017, International Journal of Climatology
  - **clusterProfiler**: Yu et al. 2012, OMICS

---

## Support and Development

- **Issues**: Report bugs or request features at https://github.com/hubner-lab/ADAPTOGENE/issues
- **Developer Guide**: See `CLAUDE.md` for development instructions
- **Contact**: For questions about methods or interpretation, consult the original publications

---

## License

This pipeline is open-source software. See LICENSE file for details.
