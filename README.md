# ADAPTOGENE

A Dockerized pipeline for population genomics: structure, association, and maladaptation.

---

## Overview

**ADAPTOGENE** automates population genomic analyses from filtered VCF to candidate genes and conservation-relevant predictions. It integrates three analytical layers: (1) population structure inference via PCA and sNMF, (2) genotype–environment (GEA) and genotype–phenotype (GWAS) association mapping with region clustering, gene annotation, and GO enrichment, and (3) maladaptation assessment via Gradient Forest genetic offset modeling under future climate scenarios. All dependencies ship in a single Docker image; climate data is downloaded automatically.

---

## Key Features

- Unified framework for **population structure, association, and maladaptation**
- Automated **climate data download** (WorldClim bioclimatic variables, CMIP6 projections)
- **EMMAX + LFMM** association methods with configurable combination strategies
- **Region clustering → gene annotation → GO enrichment** pipeline with per-region visualizations
- Detection of **overlapping GEA/GWAS signals** across traits and methods
- **Interactive Shiny dashboard** for exploring all results
- **Pipeline summary table** accumulated across all modes
- Fully reproducible via **Docker + Snakemake**

---

## Requirements and Installation

Requires **Docker** and a Unix-like OS (Linux/macOS). All software dependencies are inside the Docker image.

```bash
git clone https://github.com/hubner-lab/ADAPTOGENE.git
cd ADAPTOGENE
docker build -t adaptogene .
```

---

## Input Data

| File | Format | Required | Description |
|------|--------|----------|-------------|
| Genotypes | `.vcf` / `.vcf.gz` | Yes | Biallelic SNP data |
| Metadata | TSV | Yes | Columns: `site`, `sample`, `latitude`, `longitude` [, trait columns...] |
| Annotation | GFF3 | For association | Gene models; include GO terms in attributes for enrichment |

Chromosome names are normalized automatically during processing (`chr1` → `1`). All outputs use stripped names.

---

## Quick Start

```bash
# Run processing mode
docker run --user $(id -u):$(id -g) --rm --memory=20g -v $PWD:/pipeline adaptogene:latest \
  snakemake -c4 -s Snakefile --config mode=processing --configfile config.yaml --scheduler greedy
```

Replace `mode=processing` with subsequent modes from the table below. Use `-n` for dry runs, `-R <rule>` to force a rule and its downstream, `-p` to print shell commands.

For interactive development:

```bash
docker run -w /pipeline --user $(id -u):$(id -g) -it -v $PWD:/pipeline adaptogene bash
snakemake -s Snakefile --configfile config.yaml --config mode=<MODE> --cores 4 --scheduler greedy
```

---

## Pipeline Modes

Run modes sequentially. Each mode is invoked via `--config mode=<MODE>`.

| # | Mode | Purpose | Key Config | Main Outputs |
|---|------|---------|------------|--------------|
| 1 | `processing` | Filter VCF, LD prune, normalize | `FILTER_*`, `LD_*` | Filtered VCF, LEA formats |
| 2 | `structure` | PCA, sNMF ancestry (K range) | `SNMF_K_START/END` | PCA plots, cross-entropy, Q-matrices |
| 3 | `structure_K` | Impute, climate download, piemaps | `SNMF_K_BEST`, `MAP_*` | Piemaps, climate tables, pop stats |
| 4 | `association` | GEA: EMMAX/LFMM → regions → genes → GO | `ASSOC_*`, `ENRICHMENT_*` | Manhattan plots, regions, genes, enrichment |
| 5 | `association_phenotypes` | GWAS on metadata traits (cols 5+) | `PHENO_*` | Manhattan, piemaps, regions, genes |
| 6 | `overlapping` | GEA ∩ GWAS overlap analysis | `OVERLAP_*` | Miami plot, overlap stats, combined regions |
| 7 | `regionplot` | Regional Manhattan + gene annotations | `REGIONPLOT_*` | Publication-ready region figures |
| 8 | `maladaptation` | Gradient Forest → genetic offset | `GF_*`, `FUTURE_*` | Offset maps, importance plots, site scores |
| 9 | `haplotype_scan` | crosshap epsilon scanning | `HAPLOTYPE_SCAN_*` | HapObjects, clustree plots |
| 10 | `haplotype` | Haplotype viz at selected epsilon | `HAPLOTYPE_EPSILON_SELECTED` | crosshap_viz, boxplots, piemaps |

> **After mode 2**: inspect the cross-entropy plot and set `SNMF_K_BEST` in `config.yaml` before continuing.
> **After mode 9**: inspect clustree plots and set `HAPLOTYPE_EPSILON_SELECTED` before running mode 10.

### Pipeline Flow

```mermaid
flowchart TB
    classDef input fill:#FFFDE7,stroke:#F9A825,color:#333
    classDef preproc fill:#E3F2FD,stroke:#1565C0,color:#333
    classDef struct fill:#E8F5E9,stroke:#2E7D32,color:#333
    classDef assoc fill:#FFF3E0,stroke:#E65100,color:#333
    classDef malad fill:#FFEBEE,stroke:#C62828,color:#333

    VCF["VCF"]:::input
    META["Metadata"]:::input
    GFF["GFF3"]:::input

    proc["1 · processing\nFilter, LD prune, normalize"]:::preproc
    struct["2 · structure\nPCA, sNMF (K range)"]:::struct
    structK["3 · structure_K\nImpute, climate, piemaps"]:::struct
    assoc["4 · association\nGEA → regions → genes → GO"]:::assoc
    pheno["5 · association_phenotypes\nGWAS → regions → genes → GO"]:::assoc
    overlap["6 · overlapping\nGEA ∩ GWAS → combined regions"]:::assoc
    regplot["7 · regionplot\nRegional Manhattan + genes"]:::assoc
    malad["8 · maladaptation\nGradient Forest → offset"]:::malad
    hapscan["9 · haplotype_scan\ncrosshap epsilon scan"]:::assoc
    hapviz["10 · haplotype\ncrosshap viz + boxplots"]:::assoc

    VCF & META & GFF --> proc --> struct
    struct -. "Set K_BEST" .-> structK
    structK --> assoc & pheno
    assoc --> regplot & malad
    assoc & pheno --> overlap
    assoc & pheno & overlap -.-> hapscan
    hapscan -. "Set epsilon" .-> hapviz
```

### Detailed Rule DAGs

<details>
<summary><b>Processing Mode</b> — VCF filtering, LD pruning, format conversion</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5

    VCF[/"Input VCF"/]:::data
    META[/"Metadata"/]:::data
    GFF[/"GFF3"/]:::data

    extract_samples["extract_samples\nSample IDs from metadata"]
    calc_missing["calculate_sample_missing\nPer-sample missingness filter"]
    filter_vcf["filter_vcf\nMAF + SNP miss + chr normalization"]
    ld_prune["ld_prune\nLD pruning (sliding window)"]
    sample_order["extract_vcf_sample_order\nSample order from VCF header"]
    align_meta["align_metadata\nMatch metadata to VCF order"]
    norm_gff["normalize_gff\nStrip chr prefix from GFF"]
    vcf2lfmm["vcf_to_lfmm\nVCF → .geno + .lfmm + .vcfsnp"]
    summary["write_summary\nLog processing stats"]

    VCF --> calc_missing
    META --> extract_samples
    extract_samples --> calc_missing --> filter_vcf
    filter_vcf --> ld_prune --> vcf2lfmm
    filter_vcf --> sample_order --> align_meta
    META --> align_meta
    GFF --> norm_gff
    filter_vcf & ld_prune & extract_samples & calc_missing --> summary
```
</details>

<details>
<summary><b>Structure Mode</b> — Population structure inference</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    GENO[/".geno (LD-pruned)"/]:::data
    LFMM[/".lfmm (LD-pruned)"/]:::data
    META[/"Aligned Metadata"/]:::data

    snmf["snmf\nsNMF ancestry K_START..K_END"]
    pca["pca_plot\nPCA + Tracy-Widom test"]
    cross_ent["cross_entropy_plot\nK selection guide"]
    clusters["extract_clusters ×K\nQ-matrix per K"]
    struct_bar["structure_barplot ×K\nAncestry proportions"]
    pca_struct["pca_structure_plot ×K\nPCA with pie charts"]
    pop_diff["pop_diff_test ×K\nPopulation differentiation"]
    summary["write_summary"]

    GENO --> snmf
    LFMM --> pca
    META --> pca & clusters
    snmf --> cross_ent & clusters & pop_diff
    clusters --> struct_bar & pca_struct
    pca -- "projections\neigenvalues" --> pca_struct
    cross_ent --> summary
```
</details>

<details>
<summary><b>Structure K Mode</b> — Imputation, climate, population statistics</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef optional fill:#F3E5F5,stroke:#9C27B0,color:#333,stroke-dasharray:3 3

    SNMF[/"sNMF project"/]:::data
    LFMM[/".lfmm"/]:::data
    VCF_LD[/"LD-pruned VCF"/]:::data
    VCF_FILT[/"Filtered VCF"/]:::data
    META[/"Aligned Metadata"/]:::data
    CLUST[/"Clusters (K_BEST)"/]:::data

    impute_ld["impute_ld\nImpute LD-pruned (sNMF Q)"]
    lfmm2vcf["lfmm2vcf_ld\nImputed LFMM → VCF"]
    dl_climate["download_climate_present\nWorldClim bioclim rasters"]
    density["density_plot\nClimate density curves"]
    corr_heat["correlation_heatmap\nClimate × trait correlations"]
    piemap_s["piemap_simple ×bio\nPieMaps (uniform pie size)"]

    tajima["tajima_d\nTajima D per population"]:::optional
    pi_div["pi_diversity\nPi diversity per population"]:::optional
    ibd_rule["ibd\nIsolation by distance"]:::optional
    mantel["mantel_test\nIBD/IBE Mantel test"]:::optional
    amova_rule["amova\nAMOVA analysis"]:::optional
    piemap_t["piemap_plot ×bio\nPieMaps (trait-scaled)"]:::optional
    summary["write_summary"]

    SNMF --> impute_ld
    LFMM --> impute_ld
    impute_ld --> lfmm2vcf
    VCF_LD --> lfmm2vcf
    META --> dl_climate
    dl_climate --> density & corr_heat & piemap_s & piemap_t & mantel & summary
    META --> corr_heat & tajima & pi_div & ibd_rule
    CLUST --> piemap_s & piemap_t & ibd_rule & mantel
    VCF_FILT --> tajima & pi_div
    tajima & pi_div --> piemap_t
    lfmm2vcf --> amova_rule
    META --> amova_rule
```

Purple nodes = optional, require `POP_CALC_STATS: TRUE`.
</details>

<details>
<summary><b>Association Mode (GEA)</b> — Climate association analysis</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef enrich fill:#E8F5E9,stroke:#2E7D32,color:#333,stroke-dasharray:3 3

    VCF_FILT[/"Filtered VCF"/]:::data
    CLIMATE[/"Climate (scaled)"/]:::data
    PCA[/"PCA projections"/]:::data
    META[/"Aligned Metadata"/]:::data
    LFMM_IMP[/"Imputed LFMM (LD)"/]:::data
    GFF[/"Normalized GFF"/]:::data

    vcf2lfmm_full["vcf_to_lfmm_full\nFull VCF → LEA formats"]
    snmf_full["snmf_full\nsNMF on full dataset"]
    impute_full["impute_full\nImpute full dataset"]
    lfmm2vcf_full["lfmm2vcf_full\nImputed → VCF"]
    emmax["emmax_analysis\nEMMAX: kinship + PCA covariates"]
    lfmm_rule["lfmm_analysis\nLFMM: latent factor model"]
    sig_snps["find_sig_snps ×method\nPer-method significant SNPs"]
    combine["combine_selected_snps\nMerge methods (Sum/Overlap)"]
    regions["create_regions\nCluster SNPs → per-trait + combined"]
    genes["find_genes_around_regions\nGenes in per-trait regions"]
    genes_c["find_genes_combined_regions\nGenes in combined regions"]
    enrichment["run_enrichment\nGO enrichment per region"]:::enrich
    enrich_plots["plot_enrichment\nDotplot / emapplot / cnetplot"]:::enrich
    manh_simple["manhattan_plot ×method×trait\nPer-trait Manhattan"]
    manh_regions["manhattan_plot_regions ×method×trait\nManhattan + regions highlighted"]
    manh_combined["manhattan_combined\nAll traits + methods combined"]
    summary["write_summary"]

    VCF_FILT --> vcf2lfmm_full --> snmf_full --> impute_full
    vcf2lfmm_full -- ".lfmm" --> impute_full
    VCF_FILT --> lfmm2vcf_full
    impute_full --> lfmm2vcf_full
    VCF_FILT --> emmax
    CLIMATE --> emmax & lfmm_rule
    PCA --> emmax
    META --> emmax
    LFMM_IMP --> lfmm_rule
    impute_full --> lfmm_rule
    vcf2lfmm_full -- ".vcfsnp" --> lfmm_rule & genes & genes_c
    emmax & lfmm_rule --> sig_snps --> combine --> regions
    regions --> genes & genes_c & manh_regions & manh_combined
    GFF --> genes & genes_c & enrichment
    genes --> enrichment --> enrich_plots
    emmax & lfmm_rule --> manh_simple
    combine & regions & genes --> summary
```
</details>

<details>
<summary><b>Phenotype Association Mode</b> — GWAS on phenotypic traits</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef enrich fill:#E8F5E9,stroke:#2E7D32,color:#333,stroke-dasharray:3 3
    classDef drop fill:#FFEBEE,stroke:#C62828,color:#333,stroke-dasharray:3 3

    META[/"Aligned Metadata"/]:::data
    VCF_FILT[/"Filtered VCF"/]:::data
    PCA[/"PCA projections"/]:::data
    GFF[/"Normalized GFF"/]:::data
    VCFSNP[/".vcfsnp (full)"/]:::data
    CLIMATE[/"Climate raster"/]:::data
    CLUST[/"Clusters (K_BEST)"/]:::data

    prep["prepare_phenotypes\nExtract traits, handle missing (MEAN/MEDIAN/DROP)"]

    subgraph pathA["Path A: MEAN/MEDIAN"]
        tped_a["tped_pheno\nVCF → TPED/TFAM"]
        kin_a["kinship_pheno\nCompute kinship"]
        emmax_a["emmax_pheno\nEMMAX all traits at once"]
    end

    subgraph pathB["Path B: DROP"]
        subset["subset_vcf_pheno ×trait\nPer-trait VCF subset"]:::drop
        tped_b["tped_pheno_trait ×trait\nPer-trait TPED"]:::drop
        kin_b["kinship_pheno_trait ×trait\nPer-trait kinship"]:::drop
        emmax_b["emmax_pheno_trait ×trait\nEMMAX per trait"]:::drop
        combine_pv["combine_pheno_pvalues\nMerge per-trait p-values"]:::drop
    end

    sig["find_sig_snps_pheno\nSignificant phenotype SNPs"]
    combine_sel["combine_selected_snps_pheno\nMerge methods"]
    reg["create_regions_pheno\nCluster phenotype SNPs"]
    genes_p["find_genes_pheno\nGenes in per-trait regions"]
    genes_pc["find_genes_combined_pheno\nGenes in combined regions"]
    enrichment_p["run_enrichment_pheno"]:::enrich
    enrich_plots_p["plot_enrichment_pheno"]:::enrich
    manh_p["manhattan_pheno ×trait"]
    manh_pr["manhattan_pheno_regions ×trait"]
    manh_pc["manhattan_combined_pheno\nAll traits combined"]
    piemap_p["piemap_pheno ×trait\nTrait PieMaps"]
    summary["write_summary"]

    META --> prep
    VCF_FILT --> tped_a --> kin_a --> emmax_a
    PCA --> emmax_a & emmax_b
    prep --> emmax_a

    VCF_FILT & prep --> subset --> tped_b --> kin_b --> emmax_b
    prep --> emmax_b
    emmax_b --> combine_pv

    emmax_a & combine_pv --> sig --> combine_sel --> reg
    reg --> genes_p & genes_pc & manh_pr & manh_pc
    GFF --> genes_p & genes_pc & enrichment_p
    VCFSNP --> genes_p & genes_pc
    genes_p --> enrichment_p --> enrich_plots_p
    emmax_a & combine_pv --> manh_p
    CLIMATE --> piemap_p
    CLUST --> piemap_p
    prep --> piemap_p
    prep & combine_sel & reg & genes_p --> summary
```

Red nodes = DROP path only (per-trait subsetting). Either Path A or Path B runs, not both.
</details>

<details>
<summary><b>Overlapping Mode</b> — GEA + GWAS combined</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef enrich fill:#E8F5E9,stroke:#2E7D32,color:#333,stroke-dasharray:3 3

    GEA_SNPS[/"GEA Selected_SNPs"/]:::data
    GWAS_SNPS[/"GWAS Selected_SNPs"/]:::data
    GEA_REG[/"GEA Regions_combined"/]:::data
    GWAS_REG[/"GWAS Regions_combined"/]:::data
    GFF[/"Normalized GFF"/]:::data
    VCFSNP[/".vcfsnp (full)"/]:::data

    combine["combine_overlap_snps\nMerge GEA+GWAS SNPs, compute overlaps\n→ new combined regions"]
    genes_o["find_genes_overlap\nGenes in per-trait regions"]
    genes_oc["find_genes_combined_overlap\nGenes in combined regions"]
    enrich_o["run_enrichment_overlap"]:::enrich
    enrich_plots_o["plot_enrichment_overlap"]:::enrich
    summary["write_summary"]

    GEA_SNPS & GWAS_SNPS --> combine
    GEA_REG & GWAS_REG --> combine
    combine --> genes_o & genes_oc
    GFF --> genes_o & genes_oc & enrich_o
    VCFSNP --> genes_o & genes_oc
    genes_o --> enrich_o --> enrich_plots_o
    combine & genes_o & enrich_o --> summary
```
</details>

<details>
<summary><b>Regionplot Mode</b> — Regional Manhattan with gene annotations</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5

    GFF[/"Normalized GFF"/]:::data
    REGIONS[/"Regions_per_trait"/]:::data
    PVALS[/"Association p-values\n(all methods)"/]:::data

    gff2topr["gff2topr\nGFF → topr annotation format"]
    regplot["regionplot\nRegional Manhattan + gene overlays"]

    GFF --> gff2topr --> regplot
    REGIONS --> regplot
    PVALS --> regplot
```
</details>

<details>
<summary><b>Maladaptation Mode</b> — Gradient Forest and genetic offset</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef optional fill:#F3E5F5,stroke:#9C27B0,color:#333,stroke-dasharray:3 3

    META[/"Aligned Metadata"/]:::data
    CLIMATE[/"Climate present (raster + tables)"/]:::data
    LFMM_FULL[/".lfmm (full)"/]:::data
    VCFSNP[/".vcfsnp (full)"/]:::data
    SIG_SNPS[/"Selected_SNPs"/]:::data
    CLUST[/"Clusters (K_BEST)"/]:::data
    TAJIMA[/"Tajima D"/]:::data
    PI[/"Pi Diversity"/]:::data

    dl_future["download_climate_future_model ×model\nCMIP6 future climate per model"]
    merge["merge_climate_future\nEnsemble average across models"]
    gf_adapt["gradient_forest_adaptive\nAdaptive GF model (sig SNPs)"]
    gf_random["gradient_forest_random\nNeutral GF model (random SNPs)"]:::optional
    offset["gradient_forest_offset\nGenetic offset present → future"]
    cumimp["plot_gf_cumimp\nCumulative importance curves"]
    importance["plot_gf_importance\nOverall R²-weighted importance"]
    piemap_go["plot_gf_offset_piemap\nOffset PieMap (uniform)"]
    piemap_taj["plot_gf_offset_piemap_tajima\nOffset PieMap + Tajima D"]:::optional
    piemap_pi["plot_gf_offset_piemap_diversity\nOffset PieMap + Pi Diversity"]:::optional
    density_f["density_plot_future\nFuture climate density"]
    summary["write_summary"]

    META --> dl_future --> merge
    CLIMATE --> merge & gf_adapt & gf_random & offset
    LFMM_FULL & VCFSNP --> gf_adapt & gf_random
    SIG_SNPS --> gf_adapt & gf_random
    META --> gf_adapt & gf_random
    gf_adapt --> offset & cumimp & importance
    gf_random --> cumimp & importance
    merge --> offset & density_f
    offset --> piemap_go & piemap_taj & piemap_pi
    CLUST --> piemap_go & piemap_taj & piemap_pi
    TAJIMA --> piemap_taj
    PI --> piemap_pi
    gf_adapt & offset --> summary
```

Purple nodes = optional (random model requires `GF_RANDOM_MODEL: TRUE`; Tajima/Pi piemaps require `POP_CALC_STATS: TRUE`).
</details>

---

## Configuration

All parameters are defined in `config.yaml`. See the **[Configuration Reference](docs/configuration.md)** for full parameter documentation, YAML examples, GO enrichment details, and output directory structure.

---

## Interactive Results Viewer

ADAPTOGENE includes a Shiny dashboard that auto-discovers all projects and organizes outputs into tabs: Structure, Structure K, Association, Phenotype Association, Overlapping Regions, Maladaptation. Interactive Manhattan plots support click-to-select-region with linked gene and GO enrichment tables.

```bash
docker run --user $(id -u):$(id -g) --rm -p 3838:3838 -v $PWD:/pipeline adaptogene:latest \
  Rscript -e "shiny::runApp('/pipeline/scripts/app.R', host='0.0.0.0', port=3838)"
```

Open http://localhost:3838 in your browser.

---

## Citation

If you use ADAPTOGENE in published work, please cite:

- This repository: `https://github.com/hubner-lab/ADAPTOGENE`
- Underlying methods:
  - **sNMF**: Frichot et al. 2014, *Molecular Biology and Evolution*
  - **EMMAX**: Kang et al. 2010, *Nature Genetics*
  - **LFMM**: Frichot et al. 2013, *Molecular Biology and Evolution*
  - **Gradient Forest**: Ellis et al. 2012, *Ecology*
  - **WorldClim**: Fick & Hijmans 2017, *International Journal of Climatology*
  - **clusterProfiler**: Yu et al. 2012, *OMICS*

---

## License

This pipeline is open-source software. See LICENSE file for details.

