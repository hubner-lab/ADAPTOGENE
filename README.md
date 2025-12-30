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

## Key features

- Unified framework for **population structure, association, and maladaptation**
- Automated **climate data download and preprocessing**
- Explicit detection of **overlapping association signals**
- Reproducible execution via **Docker**
- Designed for **landscape and conservation genomics**
- Minimal manual intervention once configured

---

## Requirements

- Docker  
- Unix-like operating system (Linux/macOS recommended)

All software dependencies are included inside the Docker image.

---

## Installation

```bash
git clone https://github.com/hubner-lab/ADAPTOGENE.git
cd ADAPTOGENE
```

---

## Input data

All input files must be placed in the directory specified by `INDIR` in `config.yaml`.
Pipeline should be run in the ADAPTOGENE directory, so it will find all underlying scripts, in the same directory you should create INDIR directory with all input files which you specify in config.yaml file.

### Required inputs

- **VCF** — genotype data  
- **Sample metadata** — sample IDs, populations, geographic coordinates  
- **GFF3** — genome annotation  

### Optional inputs

- Phenotypic trait tables  
- Custom population-level traits  

---

## Configuration

All pipeline parameters are defined in `config.yaml`, including:

- Variant filtering thresholds  
- Population structure parameters  
- Association methods and correction strategies  
- Climate and maladaptation settings  

Only the configuration file needs to be modified between runs.

---

## Running ADAPTOGENE

### 1. Build Docker image

```bash
docker build -t adaptogene .
```

### 2. Enter Docker container

```bash
docker run -w /pipeline --user $(id -u):$(id -g) -it -v $PWD:/pipeline adaptogene bash
```

### 3. Run Snakemake inside docker 

Example command:

```bash
snakemake -s Snakefile --configfile config.yaml --config mode=structure --cores 24  --scheduler greedy
```


---

## Pipeline modes

### 1. Population structure (`mode=structure`)

**Input**
- Filtered VCF  
- Sample metadata  

**Purpose**
- LD pruning  
- Principal Component Analysis (PCA)  
- Ancestry inference using sNMF across a range of K values  

**Output**
- PCA plots  
- Cross-entropy plots for K selection  
- Ancestry proportion matrices  

**User action**
- Inspect cross-entropy plots  
- Select the optimal K  
- Set `K_BEST` in `config.yaml`  

---

### 2. Spatial structure analysis (`mode=structure_K`)

**Input**
- Selected `K_BEST`  
- Geographic coordinates  
- Optional custom traits  

**Purpose**
- Visualize ancestry proportions in geographic space  
- Assess spatial distribution of genetic clusters  
- Compute population diversity statistics (π, Tajima’s D)  

**Output**
- Geographic ancestry maps  
- Population-level diversity metrics  

---

### 3. Association analysis (`mode=association`)

**Input**
- VCF  
- Environmental predictors and/or phenotypic traits  
- GFF annotation  

**Methods**
- EMMAX  
- LFMM (optional)  
- Multiple testing correction (Bonferroni, q-value, top-N)  

**Purpose**
- Identify SNPs associated with traits or environmental variables  
- Aggregate signals across methods  
- Detect overlapping loci across traits  

**Output**
- Significant SNP tables  
- Candidate gene lists  
- Overlap summaries across traits and methods  

---

### 4. Regional visualization (`mode=regionplot`) *(optional)*

**Purpose**
- Visualize SNP–gene relationships in selected genomic regions  
- Highlight candidate genes and functional annotations  

**Output**
- Publication-ready regional plots  

---

### 5. Adaptive potential and maladaptation (`mode=maladaptation`)

**Input**
- Adaptive SNP set from association analysis  
- Present and future climate layers (automatically downloaded)  
- Population coordinates  

**Methods**
- Gradient Forest  
- Climate scenarios (SSP)  

**Purpose**
- Model genotype–environment relationships  
- Estimate genomic offset under future climate conditions  
- Rank populations by maladaptation risk  

**Output**
- Adaptive potential maps  
- Population-level maladaptation scores  
- Population prioritization tables for conservation planning  

---

## Intended applications

- Landscape genomics  
- Climate adaptation studies  
- Conservation genomics  
- Multi-trait association analysis  

---


## Citation

If you use ADAPTOGENE in published work, please cite:
- This repository  
- sNMF, EMMAX, LFMM, Gradient Forest, and other underlying methods  

---

## Contact

For questions, bug reports, or contributions, please use the GitHub issue tracker.
