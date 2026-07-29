# ADAPTOGENE

A Dockerized pipeline for population genomics: structure, association, and maladaptation.

---

## Overview

**ADAPTOGENE** automates population genomic analyses from filtered VCF to candidate genes and conservation-relevant predictions. It integrates three analytical layers: (1) population structure inference via PCA and sNMF, (2) genotype–environment (GEA) and genotype–phenotype (GWAS) association mapping with region clustering, gene annotation, and GO enrichment, and (3) maladaptation assessment via Gradient Forest genetic offset modeling under future climate scenarios. All dependencies ship in a single Docker image; climate data is downloaded automatically.

---

## Key Features

- Unified framework for **population structure, association, and maladaptation**
- Automated **climate data download** (WorldClim bioclimatic variables, CMIP6 projections)
- **EMMAX + LFMM + GAPIT3** (8 models) association methods with configurable combination strategies
- **Region clustering → gene annotation → GO enrichment** with on-demand interactive visualization
- Detection of **overlapping GEA/GWAS signals** across traits and methods
- **Interactive Shiny dashboard** for exploring all results (region-centric, on-demand enrichment/regionplot/haplotype)
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
| 1 | `processing` | Filter VCF, LD prune, normalize | `Filter.*`, `LD.*` | Filtered VCF, LEA formats |
| 2 | `prestructure` | PCA, sNMF ancestry (K range) | `sNMF.k_start`, `sNMF.k_end` | PCA plots, cross-entropy, Q-matrices |
| 3 | `structure` | Impute, climate download, piemaps | `sNMF.k_best`, `Map.*` | Piemaps, climate tables, pop stats, LD decay |
| 4 | `climate` | Predictor characterization: correlation/density, invariant-predictor detection, dbMEM + variance partitioning | `Climate.*`, `Climate.Varpart.*`, `Climate.dbMEM.*` | Correlation heatmap, density plots, dbMEM eigenvectors, varpart fractions |
| 5 | `pregea` | Hyperparameter exploration (LD-pruned): LFMM-K / EMMAX-#PC / RDA-Condition()-PC ladders (one shared PC range) | `PreGEA.*` | Ladder grids, RDA per-model artifacts, one recommendation per (method, param) |
| 6 | `gea` | GEA: EMMAX/LFMM/GAPIT → regions → genes | `GEA.*`, `GFF.*`, `Enrichment.*` | Manhattan plots, regions, genes |
| 7 | `gwas` | GWAS on metadata traits (cols 5+) | `GWAS.*` | Manhattan, piemaps, regions, genes |
| 8 | `gea_x_gwas` | GEA + GWAS overlap analysis | `GEAxGWAS.*` | Miami plot, pairwise overlap tables |
| 9 | `maladaptation` | Gradient Forest → genetic offset | `Maladaptation.*`, `Future.*` | Offset maps, importance plots, site scores |

> **After mode 2**: inspect the cross-entropy plot and set `sNMF.k_best` in `config.yaml` before continuing.

> **Mode 4 (`climate`)** characterizes the environmental predictors (correlation/density curves, invariant-predictor detection) and computes the shared spatial artifacts — dbMEM eigenvectors and climate/structure/geography variance partitioning — reused by both `pregea` and the spatial Gradient Forest.

> **Mode 5 (`pregea`) is optional but recommended** before mode 6: it sweeps LFMM-K, EMMAX-#PC, and RDA Condition()-PC (one shared PC range) on the cheap LD-pruned dataset and writes one recommended value per (method, hyperparameter), surfaced with pre-fill/Apply badges in the GEA tab's method editor — read it before committing to the expensive full-SNP `gea` run.

> **Run `mode=climate` before requesting a spatial (or `both`) Gradient Forest** in `mode=maladaptation`: the GF spatial variant reads `climate/tables/{spatial,varpart}/` (dbMEM vectors + forward-selection mask). This ordering is no longer enforced by a config-time guard — a missing `climate` run surfaces only as a Snakemake missing-input error.

> **Regionplot, GO enrichment, haplotype scanning, and haplotype visualization** are computed on-demand in the Shiny dashboard rather than as pipeline modes.

### Fixed constants (PreGEA / Climate)

The preGEA/climate split slimmed `PreGEA.*` from ~30 fields to ~10. The following values are now fixed in code (`workflow/rules/common.smk`) rather than exposed as config keys:

| Constant | Value | Why |
|----------|-------|-----|
| PreGEA random seed | `42` | Reproducible ladders; no longer a `PreGEA.seed` key |
| EMMAX/GAPIT kinship | `BN` | EMMAX and GAPIT share the one BN matrix from `emmax-kin`; the old `EMMAX.kinship` IBS alternative was never exercised |
| PreGEA LFMM `genomic_control` | `FALSE` | Keeps λ_GC informative for K selection (GC forces λ≈1). The production run has its own switch, `GEA.configs[LFMM].params.genomic_control`, default `TRUE` |
| LFMM/EMMAX hit-count FDR | `0.1` | Hit-count vs K/#PC is a trend read, not a threshold decision. Bonferroni is still computed into the ladder TSV but no longer plotted |
| EMMAX minimum #PCs | `0` | `n_pcs=0` is the kinship-only reference panel at the bottom of the #PC ladder |
| RDA minimum predictors | `2` | RDA cannot fit with fewer than 2 predictors — a hard requirement, not a preference |
| Recommender λ tolerance | `0.15` | Width of the recommender's λ calibration band; the user reads the actual λ plot to decide |
| Recommender deflation floor | `0.90` | Lower λ bound the recommender flags as over-correction |
| `ordiR2step` `Pin` / `R2permutations` | `0.01` / `999` | Blanchet double-stopping-rule standard for forward selection (Climate varpart + RDA setup) |

Removed switches and consolidations (no longer config keys anywhere): LFMM/EMMAX/RDA have no per-block `enabled` switches — they always run together in one `mode=pregea` pass; `Varpart` likewise has no switch (it moved into `mode=climate` and runs unconditionally). `LFMM.k_offset_low`/`k_offset_high`/`k_min` consolidated into a single `PreGEA.k_offset` (sweep = `sNMF.k_best ± k_offset`); `RDA.condition_pcs_min`/`condition_pcs_max` consolidated into the shared `PreGEA.n_pcs_max` (EMMAX #PCs and RDA Condition()-PCs now sweep the identical `0..n_pcs_max` range). `LFMM.fdr`/`bonf_alpha`, `EMMAX.fdr`/`bonf_alpha`, `EMMAX.lambda_tol`/`deflation_floor`, and `EMMAX.n_pcs_min` are dropped per the fixed values above.

### Pipeline Flow

```mermaid
flowchart TB
    classDef input fill:#FFFDE7,stroke:#F9A825,color:#333
    classDef preproc fill:#E3F2FD,stroke:#1565C0,color:#333
    classDef struct fill:#E8F5E9,stroke:#2E7D32,color:#333
    classDef climate fill:#E0F7FA,stroke:#00838F,color:#333
    classDef pregea fill:#EDE7F6,stroke:#5E35B1,color:#333,stroke-dasharray:3 3
    classDef assoc fill:#FFF3E0,stroke:#E65100,color:#333
    classDef malad fill:#FFEBEE,stroke:#C62828,color:#333

    VCF["VCF"]:::input
    META["Metadata"]:::input
    GFF["GFF3"]:::input

    proc["1 · processing\nFilter, LD prune, normalize"]:::preproc
    prestruct["2 · prestructure\nPCA, sNMF (K range)"]:::struct
    struct["3 · structure\nImpute, climate, piemaps"]:::struct
    climate["4 · climate\nPredictor characterization + spatial varpart"]:::climate
    pregea["5 · pregea (optional)\nLFMM-K/EMMAX-#PC/RDA ladders"]:::pregea
    gea["6 · gea\nGEA → regions → genes"]:::assoc
    gwas["7 · gwas\nGWAS → regions → genes"]:::assoc
    geaxgwas["8 · gea_x_gwas\nGEA ∩ GWAS → Miami + overlap"]:::assoc
    malad["9 · maladaptation\nGradient Forest → offset"]:::malad

    VCF & META & GFF --> proc --> prestruct
    prestruct -. "Set sNMF.k_best" .-> struct
    struct -.-> climate
    struct -.-> pregea
    climate -. "dbMEM / varpart" .-> malad
    pregea -. "Recommended hyperparameters" .-> gea
    struct --> gea & gwas
    gea --> malad
    gea & gwas --> geaxgwas
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
<summary><b>PreStructure Mode</b> — Population structure inference</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    GENO[/".geno (LD-pruned)"/]:::data
    LFMM[/".lfmm (LD-pruned)"/]:::data
    META[/"Aligned Metadata"/]:::data

    snmf["snmf\nsNMF ancestry k_start..k_end"]
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

**Key outputs**:
- `PreStructure/plots/pca.png/svg`, `tracy_widom.png/svg`
- `PreStructure/plots/cross_entropy_K{start}-{end}.png/svg`
- `PreStructure/plots/K{k}/structure_K{k}.png/svg`
- `PreStructure/tables/K{k}/clusters_K{k}.tsv`
</details>

<details>
<summary><b>Structure Mode</b> — Imputation, climate, population statistics</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef optional fill:#F3E5F5,stroke:#9C27B0,color:#333,stroke-dasharray:3 3

    SNMF[/"sNMF project"/]:::data
    LFMM[/".lfmm"/]:::data
    VCF_LD[/"LD-pruned VCF"/]:::data
    VCF_FILT[/"Filtered VCF"/]:::data
    META[/"Aligned Metadata"/]:::data
    CLUST[/"Clusters (k_best)"/]:::data

    impute_ld["impute_ld\nImpute LD-pruned (sNMF Q)"]
    lfmm2vcf["lfmm2vcf_ld\nImputed LFMM → VCF"]
    dl_climate["download_climate_present\nWorldClim bioclim rasters"]
    piemap_s["piemap_simple ×bio\nPieMaps (uniform pie size)"]

    tajima["tajima_d\nTajima D per population"]:::optional
    pi_div["pi_diversity\nPi diversity per population"]:::optional
    ibd_rule["ibd\nIsolation by distance"]:::optional
    mantel["mantel_test\nIBD/IBE Mantel test"]:::optional
    amova_rule["amova\nAMOVA analysis"]:::optional
    piemap_t["piemap_plot ×bio\nPieMaps (trait-scaled)"]:::optional
    ld_decay["ld_decay_*\nLD decay analysis"]:::optional
    summary["write_summary"]

    SNMF --> impute_ld
    LFMM --> impute_ld
    impute_ld --> lfmm2vcf
    VCF_LD --> lfmm2vcf
    META --> dl_climate
    dl_climate --> piemap_s & piemap_t & mantel & summary
    META --> tajima & pi_div & ibd_rule
    CLUST --> piemap_s & piemap_t & ibd_rule & mantel
    VCF_FILT --> tajima & pi_div & ld_decay
    tajima & pi_div --> piemap_t
    lfmm2vcf --> amova_rule
    META --> amova_rule & ld_decay
```

Purple nodes = optional (`Population.calc_stats: true`).

**Key outputs**:
- `climate/rasters/present/`, `climate/tables/present/climate_present_{all,site,site_scaled}.tsv` (predictor characterization tables are produced by `mode=climate`)
- `Structure/plots/piemap/piemap_{bio}.png/svg` + `zoom/`
- `Structure/tables/pop_stats/tajima_d_by_pop.tsv`, `pi_diversity_by_pop.tsv`
- `Structure/plots/ld_decay/ld_decay_genome_wide.png/svg`
</details>

<details>
<summary><b>Climate Mode</b> — Predictor characterization and shared spatial artifacts</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef optional fill:#F3E5F5,stroke:#9C27B0,color:#333,stroke-dasharray:3 3

    CLIM_SITE[/"Climate site (present)"/]:::data
    CLIM_SCALED[/"Climate site (scaled)"/]:::data
    META[/"Aligned Metadata"/]:::data
    PCA[/"PCA projections + eigenvalues"/]:::data
    CLUST[/"Clusters (k_best)"/]:::data

    check_variance["check_climate_variance\nInvariant-predictor detection"]
    density["density_plot\nPredictor density curves"]
    density_pheno["density_plot_phenotypes\nPhenotype density curves"]:::optional
    corr_heat["correlation_heatmap\nPredictor × trait correlations"]
    dbmem["climate_dbmem\nadespatial::dbmem, MST-truncated geodesic distance"]
    varpart["climate_varpart\nclimate/structure/geography varpart + Px (Lasky et al. 2012)"]
    summary["write_summary"]

    CLIM_SITE --> check_variance & density & corr_heat
    META --> corr_heat & density_pheno & dbmem
    PCA & CLIM_SCALED & CLUST --> varpart
    dbmem --> varpart
    check_variance & corr_heat & varpart --> summary
```

Purple node = optional (`density_plot_phenotypes` needs trait columns in the metadata). The dbMEM + varpart artifacts are the shared spatial products reused by the spatial Gradient Forest (`mode=maladaptation`).

**Key outputs**:
- `climate/plots/correlation_heatmap.png`, `density_plot_present.png`, `density_plot_phenotypes.png`
- `climate/tables/present/climate_invariant_predictors.tsv`
- `climate/tables/spatial/dbmem_vectors.tsv`, `dbmem_diagnostics.tsv` + `climate/plots/spatial/dbmem_screeplot.png/svg`
- `climate/tables/varpart/varpart_fractions.tsv`, `varpart_anova.tsv`, `px_per_variable.tsv`, `dbmem_selected.tsv`
- `climate/plots/varpart/varpart_venn.png`, `varpart_fractions_bar.png`, `px_barplot.png`, `dbmem_selection_path.png`
</details>

<details>
<summary><b>PreGEA Mode</b> (optional) — LD-pruned hyperparameter exploration before the full GEA run</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef optional fill:#F3E5F5,stroke:#9C27B0,color:#333,stroke-dasharray:3 3

    VCF_LD[/"LD-pruned VCF"/]:::data
    LFMM_LD[/".lfmm (LD-pruned, imputed)"/]:::data
    CLIMATE[/"Climate (scaled)"/]:::data
    PCA[/"PCA projections + eigenvalues"/]:::data
    META[/"Aligned Metadata"/]:::data

    subset["subset_vcf_pruned_climate\nLD-pruned VCF -> climate-valid subset"]
    tped["tped_pruned\nVCF -> TPED/TFAM"]
    kinship["kinship_pruned\nBN kinship, fixed across the EMMAX ladder"]
    screeplot["pca_screeplot\nBroken-stick null + LFMM K band"]

    lfmm_rung["lfmm_rung ×K×trait\nLFMM per K (genomic_control=FALSE)"]
    lfmm_ladder["lfmm_ladder\nlambda_GC, hit counts, histogram shape -> plot grid"]

    emmax_rung["emmax_rung ×n_pcs×trait\nEMMAX per #PCs (kinship held fixed)"]
    emmax_ladder["emmax_ladder\nlambda_GC, hit counts -> plot grid"]

    rda_setup["rda_setup\nCollinearity screen + Condition()-PC model ladder + model comparison + ordiR2step path"]
    rda_models[/"rda/models/pc{n}/\nbiplot + axis screeplot + axis anova (per Condition()-PC)"/]:::data

    recommend["pregea_recommend\nOne row per (method, param): rule + evidence"]
    guard["transfer_guard\nFull-set lambda re-check at recommended hyperparameters"]:::optional
    summary["write_summary"]

    VCF_LD --> subset --> tped --> kinship
    META --> subset
    PCA --> screeplot
    kinship & tped & CLIMATE & PCA & META --> emmax_rung --> emmax_ladder
    LFMM_LD & CLIMATE --> lfmm_rung --> lfmm_ladder
    LFMM_LD & CLIMATE & PCA --> rda_setup --> rda_models
    lfmm_ladder & emmax_ladder & rda_setup --> recommend
    recommend & LFMM_LD & CLIMATE & kinship & tped --> guard
    recommend & lfmm_ladder & emmax_ladder & rda_setup --> summary
```

Purple node = opt-in (`PreGEA.TransferGuard.enabled: true`) — the one PreGEA rule that pulls the full-set imputation chain into this otherwise LD-pruned-only mode. `lambda_GC` is deliberately fit with `genomic_control=FALSE` for the LFMM ladder — with it on, `lfmm2.test`'s recalibration forces lambda toward 1 regardless of whether K is correct, so K selection uses p-value-histogram flatness instead.

**Key outputs**:
- `PreGEA/tables/{lfmm,emmax}/{lfmm,emmax}_ladder.tsv` + `PreGEA/plots/{lfmm,emmax}/*_histogram_grid.png`, `*_qq_grid.png`, `*_lambda_vs_*.png`, `*_hits_vs_*.png`
- `PreGEA/tables/rda/rda_condition_ladder.tsv`, `rda_predictor_collinearity.tsv`, `rda_ordir2step_path.tsv` + `PreGEA/plots/rda/rda_model_comparison.png`, `rda_ordir2step_path.png`
- `PreGEA/plots/rda/models/pc{n}/biplot.png/svg`, `axis_screeplot.png/svg` + `PreGEA/tables/rda/models/pc{n}/axis_anova.tsv` — one set per Condition()-PC value, selectable in the Shiny RDA tab
- `PreGEA/tables/pregea_recommendations.tsv` — read by the GEA tab's method editor for pre-fill/Apply badges
- `PreGEA/tables/pregea_transfer_guard.tsv` (opt-in only)
</details>

<details>
<summary><b>GEA Mode</b> — Climate association analysis</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5

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
    tped_gea["tped_gea\nVCF → TPED/TFAM"]
    kinship_gea["kinship_gea\nBN kinship matrix"]
    emmax["emmax_gea\nEMMAX: kinship + PCA covariates"]
    lfmm_rule["lfmm_gea\nLFMM: latent factor model"]
    gapit["gapit_gea ×model\nGAPIT3 (GLM/MLM/BLINK/FarmCPU/...)"]
    sig_snps["assoc_find_sig_snps ×method\nPer-method significant SNPs"]
    combine["assoc_combine_selected_snps\nMerge methods"]
    regions["assoc_create_regions\nCluster SNPs → per-trait + combined"]
    genes["assoc_find_genes_per_region\nGenes in per-trait regions"]
    genes_c["assoc_find_genes_combined\nGenes in combined regions"]
    manh_plot["assoc_manhattan_plot ×method×trait\nPer-trait Manhattan + QQ"]
    manh_combined["assoc_manhattan_combined\nAll traits combined"]
    summary["write_summary"]

    VCF_FILT --> vcf2lfmm_full --> snmf_full --> impute_full
    vcf2lfmm_full -- ".lfmm" --> impute_full
    VCF_FILT --> lfmm2vcf_full & tped_gea
    impute_full --> lfmm2vcf_full
    tped_gea --> kinship_gea
    VCF_FILT --> emmax
    CLIMATE --> emmax & lfmm_rule & gapit
    PCA --> emmax & gapit
    META --> emmax & gapit
    kinship_gea --> emmax & gapit
    LFMM_IMP --> lfmm_rule
    impute_full --> lfmm_rule
    vcf2lfmm_full -- ".vcfsnp" --> lfmm_rule & genes & genes_c
    emmax & lfmm_rule & gapit --> sig_snps --> combine --> regions
    regions --> genes & genes_c & manh_combined
    GFF --> genes & genes_c
    emmax & lfmm_rule & gapit --> manh_plot
    combine & regions & genes --> summary
```

**Key outputs**:
- `GEA/GAPIT_native_output/{model}/`
- `GEA/tables/methods/{method}/` — per-method p-values + significant SNPs
- `GEA/tables/selected_snps.tsv`, `regions_per_trait.tsv`, `regions_combined.tsv`
- `GEA/tables/genes_per_region.tsv`, `genes_combined.tsv`
- `GEA/plots/manhattan/{method}/manhattan_{trait}_K{k}_{adjust}.png/svg`
- `GEA/plots/manhattan/combined/manhattan_combined_K{k}.png/svg`
</details>

<details>
<summary><b>GWAS Mode</b> — Association on phenotypic traits</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef drop fill:#FFEBEE,stroke:#C62828,color:#333,stroke-dasharray:3 3

    META[/"Aligned Metadata"/]:::data
    VCF_FILT[/"Filtered VCF"/]:::data
    PCA[/"PCA projections"/]:::data
    GFF[/"Normalized GFF"/]:::data
    CLIMATE[/"Climate raster"/]:::data
    CLUST[/"Clusters (k_best)"/]:::data

    prep["prepare_phenotypes\nExtract traits, handle missing (MEAN/MEDIAN/DROP)"]

    subgraph pathA["Path A: MEAN/MEDIAN"]
        tped_a["tped_gwas\nVCF → TPED/TFAM"]
        kin_a["kinship_gwas\nCompute kinship"]
        emmax_a["emmax_gwas\nEMMAX all traits"]
        gapit_a["gapit_gwas ×model\nGAPIT3 all traits"]
    end

    subgraph pathB["Path B: DROP"]
        subset["subset_vcf_gwas ×trait\nPer-trait VCF subset"]:::drop
        tped_b["tped_gwas_trait ×trait\nPer-trait TPED"]:::drop
        kin_b["kinship_gwas_trait ×trait\nPer-trait kinship"]:::drop
        emmax_b["emmax_gwas_trait ×trait\nEMMAX per trait"]:::drop
        gapit_b["gapit_gwas_trait ×model×trait\nGAPIT3 per trait"]:::drop
        combine_pv["combine_gapit_gwas_pvalues\nMerge per-trait p-values"]:::drop
    end

    sig["assoc_find_sig_snps ×method"]
    combine_sel["assoc_combine_selected_snps"]
    reg["assoc_create_regions"]
    genes_p["assoc_find_genes_per_region"]
    genes_pc["assoc_find_genes_combined"]
    manh_p["assoc_manhattan_plot ×method×trait"]
    manh_pc["assoc_manhattan_combined"]
    piemap_p["piemap_gwas ×trait\nTrait PieMaps"]
    summary["write_summary"]

    META --> prep
    VCF_FILT --> tped_a --> kin_a --> emmax_a & gapit_a
    PCA --> emmax_a & gapit_a & emmax_b & gapit_b
    prep --> emmax_a & gapit_a
    VCF_FILT & prep --> subset --> tped_b --> kin_b --> emmax_b & gapit_b
    prep --> emmax_b & gapit_b
    emmax_b & gapit_b --> combine_pv
    emmax_a & gapit_a & combine_pv --> sig --> combine_sel --> reg
    reg --> genes_p & genes_pc & manh_pc
    GFF --> genes_p & genes_pc
    emmax_a & gapit_a & combine_pv --> manh_p
    CLIMATE --> piemap_p
    CLUST & prep --> piemap_p
    prep & combine_sel & reg & genes_p --> summary
```

Red nodes = DROP path only.

**Key outputs** (same structure as GEA but under `GWAS/`):
- `GWAS/tables/selected_snps.tsv`, `regions_per_trait.tsv`, `regions_combined.tsv`
- `GWAS/plots/manhattan/{method}/`, `GWAS/plots/manhattan/combined/`
- `GWAS/plots/piemap/phenomap_{trait}.png/svg` + `zoom/`
</details>

<details>
<summary><b>GEAxGWAS Mode</b> — GEA + GWAS overlap analysis</summary>

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5

    GEA_SNPS[/"GEA Selected_SNPs"/]:::data
    GWAS_SNPS[/"GWAS Selected_SNPs"/]:::data
    GEA_REG[/"GEA Regions_combined"/]:::data
    GWAS_REG[/"GWAS Regions_combined"/]:::data

    miami["miami_plot\nGEA above / GWAS below"]
    pairwise["compute_pairwise_overlaps\nAll-vs-all trait SNP overlap"]
    summary["write_summary"]

    GEA_SNPS & GWAS_SNPS & GEA_REG & GWAS_REG --> miami
    GEA_SNPS & GWAS_SNPS --> pairwise
    miami & pairwise --> summary
```

Region overlap detection, gene annotation, and GO enrichment are computed on-demand in the Shiny dashboard.

**Key outputs**:
- `GEAxGWAS/plots/miami_combined_K{k}.png/svg` + `_background.png` + `_coords.json`
- `GEAxGWAS/tables/pairwise_collapsed_snps.tsv`, `pairwise_overlap_table.tsv`
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
    CLUST[/"Clusters (k_best)"/]:::data
    TAJIMA[/"Tajima D"/]:::data
    PI[/"Pi Diversity"/]:::data

    dl_future["download_climate_future_model ×model\nCMIP6 future climate per model"]
    merge["merge_climate_future\nEnsemble average across models"]
    combine_gf["combine_gf_snps\nFilter + combine adaptive SNPs"]
    gf_adapt["gradient_forest_adaptive\nAdaptive GF model (sig SNPs)"]
    gf_random["gradient_forest_random\nNeutral GF model (random SNPs)"]:::optional
    offset["gradient_forest_offset\nGenetic offset present → future"]
    cumimp["plot_gf_cumimp\nCumulative importance curves"]
    importance["plot_gf_importance\nOverall R²-weighted importance"]
    piemap["plot_gf_offset_piemap\n{size_trait}: notrait | tajima_d | pi_diversity"]
    density_f["density_plot_future\nFuture climate density"]
    summary["write_summary"]

    META --> dl_future --> merge
    SIG_SNPS --> combine_gf --> gf_adapt & gf_random
    CLIMATE --> merge & gf_adapt & gf_random & offset
    LFMM_FULL & VCFSNP --> gf_adapt & gf_random
    META --> gf_adapt & gf_random
    gf_adapt --> offset & cumimp & importance
    gf_random --> cumimp & importance
    merge --> offset & density_f
    offset --> piemap
    CLUST & TAJIMA & PI --> piemap
    gf_adapt & offset --> summary
```

Purple nodes = optional (`Maladaptation.methods.gradient_forest.random_model: true`; `Population.calc_stats: true` for Tajima/Pi piemaps).

**Key outputs** (under `Maladaptation/{plots,tables}/gradient_forest/{run_label}_{spatial_tag}/`):
- `Maladaptation/plots/gradient_forest/{suffix}/cumulative_importance.png/svg`
- `Maladaptation/plots/gradient_forest/{suffix}/overall_importance.png/svg`
- `Maladaptation/plots/gradient_forest/{suffix}/genetic_offset_piemap[_{tajima_d,pi_diversity}].png/svg`
- `Maladaptation/tables/gradient_forest/{suffix}/genetic_offset_{map,site}.tsv`
</details>

---

## Configuration

All parameters are defined in `config.yaml`. See the **[Configuration Reference](docs/configuration.md)** for full parameter documentation, YAML examples, and output directory structure.

---

## Interactive Results Viewer

ADAPTOGENE includes a Shiny dashboard that auto-discovers all projects and organizes outputs into tabs: **Home**, **Processing**, **PreStructure**, **Structure**, **GEA**, **GWAS**, **GEAxGWAS**, **Maladaptation**.

The interface is **region-centric**: clicking a significant SNP in any Manhattan plot selects a genomic region and reveals a consolidated detail panel with GO enrichment plots, gene annotations, enrichment tables, regionplot, and haplotype analysis — all computed on-demand, directly below the Manhattan. Results persist to disk across sessions.

```bash
docker run --user $(id -u):$(id -g) --rm -p 3838:3838 -v $PWD:/pipeline adaptogene:latest \
  R -e "adaptogene.app::run_app(options = list(host = '0.0.0.0', port = 3838))"
```

Open http://localhost:3838 in your browser.

---

## Validation and Benchmarking

ADAPTOGENE is validated on two complementary categories of datasets:

**Empirical benchmarks** (known ground-truth loci, real organisms):
- [Benchmark Datasets](benchmarks/README.md) — Arabidopsis 1001G, Populus balsamifera, Rice RDP1, Sorghum SAP. Download and preparation scripts in `benchmarks/`.
- [Benchmark Testing Plan](benchmarks/TESTING_PLAN.md) — step-by-step execution sequence for all modes.

**Simulation benchmarks** (full ground truth: causal loci + fitness):
- [Simulation Benchmark Guide](docs/laruson_benchmark.md) — end-to-end validation on Láruson et al. 2022 (SLiM simulations, Dryad DOI: 10.5061/dryad.x95x69pkk). Tests both offset-vs-fitness correlation and causal-SNP detection rates.
- [Headless Automation Roadmap](docs/laruson_automation_roadmap.md) — specification for running the full benchmark without the Shiny GUI, including the custom-environment source, headless SNP-set promotion, run-all driver, and two-axis evaluation harness.

---

## Citation

If you use ADAPTOGENE in published work, please cite:

- This repository: `https://github.com/hubner-lab/ADAPTOGENE`
- Underlying methods:
  - **sNMF**: Frichot et al. 2014, *Molecular Biology and Evolution*
  - **EMMAX**: Kang et al. 2010, *Nature Genetics*
  - **LFMM**: Frichot et al. 2013, *Molecular Biology and Evolution*
  - **GAPIT3**: Wang & Zhang 2021, *Genomics, Proteomics & Bioinformatics*
  - **Gradient Forest**: Ellis et al. 2012, *Ecology*
  - **WorldClim**: Fick & Hijmans 2017, *International Journal of Climatology*
  - **clusterProfiler**: Yu et al. 2012, *OMICS*

---

## License

This pipeline is open-source software. See LICENSE file for details.
