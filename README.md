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
  - Cnetplot - Gene-concept network identifying hub genes
- **Per-region and per-trait** functional enrichment reporting
- **Dual region outputs** - Per-trait regions for specific analysis + combined climate regions for overview
- **Combined Manhattan plots** - Visualize all traits and methods simultaneously
- Reproducible execution via **Docker** (all dependencies included)
- Designed for **landscape and conservation genomics**
- Minimal manual intervention once configured

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
```

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
   - Per-trait regions (SNPs clustered separately for each climate variable)
   - Combined climate regions (all significant SNPs clustered together)
6. **Gene annotation**: Find genes within extended regions (±ASSOC_REGION_DISTANCE)
7. **GO enrichment analysis** (if ASSOC_GO_FIELD is set):
   - Per-region enrichment testing
   - Organized by trait
   - Three visualization types per region

**Output**

**Tables** (`tables/association/`):
- `Selected_SNPs.tsv` - All significant SNPs (combined across methods)
- `Regions_per_trait.tsv` - Genomic regions per climate trait
- `Regions_climate_combined.tsv` - Combined climate-associated regions
- `Genes_per_region.tsv` - Genes within per-trait regions
- `Genes_per_region_climate_combined.tsv` - Genes within combined regions
- `enrichment/{trait}/Region_{region_id}_enrichment.tsv` - Per-region GO results
- `enrichment/{trait}/Enrichment_{trait}_summary.tsv` - Per-trait summaries

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

### 5. Regional Visualization (`mode=regionplot`)

**Purpose**: Detailed visualization of specific genomic regions

**Input**
- Association results
- GFF annotation
- Top regions from association analysis

**Output**
- Regional Manhattan plots with gene annotations
- Publication-ready figures for candidate regions

---

### 6. Adaptive Potential and Maladaptation (`mode=maladaptation`)

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

---

## GO Enrichment Analysis

When a GFF file with GO annotations is provided (`ASSOC_GO_FIELD` parameter), the pipeline automatically performs functional enrichment analysis on genes associated with significant regions.

### How it Works

1. **Region Definition**: Significant SNPs are clustered into genomic regions using `ASSOC_REGION_DISTANCE`
2. **Per-Trait Regions**: SNPs are clustered separately for each climate variable (e.g., bio_1, bio_12)
3. **Gene Finding**: Genes within extended regions (region ± ASSOC_REGION_DISTANCE) are identified
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
- Useful for candidate gene prioritization
- Created when ≥2 terms are enriched

All plots saved in both PNG (300 dpi) and SVG (vector) formats for publication.

### Configuration Parameters

```yaml
ENRICHMENT_TOP_TERMS: 20        # Maximum GO terms to show in plots
ENRICHMENT_PLOT_WIDTH: 12       # Plot width in inches
ENRICHMENT_PLOT_HEIGHT: 10      # Plot height in inches
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
│   ├── EMMAX/                              # EMMAX Manhattan plots
│   ├── LFMM/                               # LFMM Manhattan plots
│   ├── association/                        # Combined Manhattan plots
│   │   ├── Manhattan_all_traits_combined_K{K}.[png,svg]
│   │   ├── Manhattan_all_traits_combined_K{K}_regions.[png,svg]
│   │   └── enrichment/{trait}/             # GO enrichment plots
│   ├── regionplot/                         # Regional plots with genes
│   └── gradientForest/                     # GF importance, offset maps
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

6. **Maladaptation**: Predict future vulnerability
   ```bash
   mode=maladaptation
   ```

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
