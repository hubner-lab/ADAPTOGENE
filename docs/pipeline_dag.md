# ADAPTOGENE Pipeline DAGs

## Overview

```mermaid
flowchart TB
    subgraph inputs["Input Files"]
        direction LR
        VCF["VCF Genotypes"]
        META["Sample Metadata"]
        GFF["GFF3 Annotation"]
    end

    subgraph proc["1 · processing"]
        direction LR
        P["Filter + LD prune + normalize"]
    end

    subgraph prestruct["2 · prestructure"]
        direction LR
        S["sNMF + PCA + cross-entropy"]
    end

    subgraph struct["3 · structure"]
        direction LR
        SK["Impute + climate + piemaps + pop stats"]
    end

    subgraph gea_mode["4 · gea"]
        direction LR
        A["EMMAX/LFMM/GAPIT → regions → genes"]
    end

    subgraph gwas_mode["5 · gwas"]
        direction LR
        AP["EMMAX/GAPIT per trait → regions → genes"]
    end

    subgraph geaxgwas_mode["6 · gea_x_gwas"]
        direction LR
        OV["Miami plot + pairwise overlap tables"]
    end

    subgraph malad["7 · maladaptation"]
        direction LR
        M["Gradient Forest → genetic offset"]
    end

    inputs --> proc --> prestruct
    prestruct -. "Set sNMF.k_best" .-> struct
    struct --> gea_mode & gwas_mode
    gea_mode --> malad
    gea_mode & gwas_mode --> geaxgwas_mode
```

Each mode is run separately via `--config mode=<MODE>`.

> **Regionplot, GO enrichment, haplotype scanning, and haplotype visualization** are computed on-demand in the Shiny dashboard rather than as pipeline modes.

---

## Detailed Rule DAGs

Node labels match Snakefile rule names exactly. `×N` indicates wildcard-expanded rules.

### Processing Mode — VCF filtering, LD pruning, format conversion

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

**Key outputs**:
- `_work/.../filtered.vcf`, `_work/.../ld.../pruned.vcf` — filtered + LD-pruned VCFs
- `_work/.../ld.../{name}.geno`, `.lfmm`, `.vcfsnp` — LEA format conversions
- `Processing/tables/metadata.tsv` — aligned metadata (sample order matches VCF)
- `Processing/tables/sample_missing_stats.tsv` — per-sample missingness
- `_intermediate/annotation/normalized.gff3` — chr-stripped GFF

---

### PreStructure Mode — Population structure inference

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
- `PreStructure/plots/pca.png/svg`, `tracy_widom.png/svg` — PCA scatter + eigenvalue test
- `PreStructure/plots/cross_entropy_K{start}-{end}.png/svg` — K selection guide
- `PreStructure/plots/K{k}/structure_K{k}.png/svg` — ancestry barplots (per K)
- `PreStructure/plots/K{k}/pca_structure_K{k}.png/svg` — PCA with ancestry pies
- `PreStructure/plots/K{k}/pop_diff_K{k}.png/svg` — population differentiation
- `PreStructure/tables/K{k}/clusters_K{k}.tsv` — Q-matrices per K

---

### Structure Mode — Imputation, climate, population statistics

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
    check_variance["check_climate_variance\nDetect invariant predictors"]
    density["density_plot\nClimate density curves"]
    corr_heat["correlation_heatmap\nClimate × trait correlations"]
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
    dl_climate --> check_variance & density & corr_heat & piemap_s & piemap_t & mantel & summary
    META --> corr_heat & tajima & pi_div & ibd_rule
    CLUST --> piemap_s & piemap_t & ibd_rule & mantel
    VCF_FILT --> tajima & pi_div & ld_decay
    tajima & pi_div --> piemap_t
    lfmm2vcf --> amova_rule
    META --> amova_rule & ld_decay
```

**Purple nodes** = optional, require `Population.calc_stats: true`; LD decay always runs.

**Key outputs**:
- `climate/rasters/present/` — WorldClim bioclim rasters (.tif)
- `climate/tables/present/climate_present_{all,site,site_scaled}.tsv`
- `climate/plots/density_plot_present.png/svg`, `correlation_heatmap.png/svg`
- `Structure/plots/piemap/piemap_{bio}.png/svg` + `zoom/` variants
- `Structure/plots/piemap/{tajima_d,pi_diversity}/piemap_{bio}_{metric}.png/svg` (optional)
- `Structure/tables/pop_stats/tajima_d_by_pop.tsv`, `pi_diversity_by_pop.tsv`, `ibd_raw.tsv`, `ibd_pairs.tsv`, `amova.tsv`
- `Structure/plots/pop_stats/mantel_test.png/svg`, `amova.png/svg`
- `Structure/plots/ld_decay/ld_decay_genome_wide.png/svg`, `ld_decay_per_chr.png/svg`, `ld_decay_per_pop.png/svg`

---

### GEA Mode — Climate association analysis

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
    combine["assoc_combine_selected_snps\nMerge methods (All/MethodOverlap)"]
    regions["assoc_create_regions\nCluster SNPs → per-trait + combined"]
    genes["assoc_find_genes_per_region\nGenes in per-trait regions"]
    genes_c["assoc_find_genes_combined\nGenes in combined regions"]
    manh_plot["assoc_manhattan_plot ×method×trait\nPer-trait Manhattan + QQ"]
    manh_combined["assoc_manhattan_combined\nAll traits + methods combined"]
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

> **GO enrichment and regionplot** are computed on-demand in the Shiny dashboard (GEA tab → region detail). Results persist to `GEA/plots/enrichment/` and `GEA/tables/enrichment/`.

**Key outputs**:
- `GEA/GAPIT_native_output/{model}/` — raw GAPIT output files
- `GEA/tables/methods/{method}/{method}_pvalues_K{k}.tsv` — per-method p-values
- `GEA/tables/methods/{method}/{method}_sig_snps_{adjust}.tsv` — significant SNPs per method
- `GEA/tables/selected_snps.tsv` — merged across methods
- `GEA/tables/regions_per_trait.tsv`, `regions_combined.tsv` — clustered regions
- `GEA/tables/genes_per_region.tsv`, `genes_per_region_collapsed.tsv`, `genes_combined.tsv`
- `GEA/plots/manhattan/{method}/manhattan_{trait}_K{k}_{adjust}.png/svg`
- `GEA/plots/manhattan/{method}/qq_{trait}_K{k}_{adjust}.png/svg`
- `GEA/plots/manhattan/combined/manhattan_combined_K{k}.png/svg`

---

### GWAS Mode — Association on phenotypic traits

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5
    classDef drop fill:#FFEBEE,stroke:#C62828,color:#333,stroke-dasharray:3 3

    META[/"Aligned Metadata"/]:::data
    VCF_FILT[/"Filtered VCF"/]:::data
    PCA[/"PCA projections"/]:::data
    GFF[/"Normalized GFF"/]:::data
    VCFSNP[/".vcfsnp (full)"/]:::data
    CLIMATE[/"Climate raster"/]:::data
    CLUST[/"Clusters (k_best)"/]:::data

    prep["prepare_phenotypes\nExtract traits, handle missing (MEAN/MEDIAN/DROP)"]

    subgraph pathA["Path A: MEAN/MEDIAN"]
        tped_a["tped_gwas\nVCF → TPED/TFAM"]
        kin_a["kinship_gwas\nCompute kinship"]
        emmax_a["emmax_gwas\nEMMAX all traits at once"]
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

    sig["assoc_find_sig_snps ×method\nSignificant phenotype SNPs"]
    combine_sel["assoc_combine_selected_snps\nMerge methods"]
    reg["assoc_create_regions\nCluster phenotype SNPs"]
    genes_p["assoc_find_genes_per_region\nGenes in per-trait regions"]
    genes_pc["assoc_find_genes_combined\nGenes in combined regions"]
    manh_p["assoc_manhattan_plot ×method×trait\nPer-trait Manhattan + QQ"]
    manh_pc["assoc_manhattan_combined\nAll traits combined"]
    piemap_p["piemap_gwas ×trait\nTrait PieMaps"]
    summary["write_summary"]

    META --> prep
    VCF_FILT --> tped_a --> kin_a --> emmax_a
    PCA --> emmax_a & emmax_b & gapit_a & gapit_b
    prep --> emmax_a & gapit_a

    VCF_FILT & prep --> subset --> tped_b --> kin_b --> emmax_b
    prep --> emmax_b & gapit_b
    emmax_b & gapit_b --> combine_pv

    emmax_a & gapit_a & combine_pv --> sig --> combine_sel --> reg
    reg --> genes_p & genes_pc & manh_pc
    GFF --> genes_p & genes_pc
    VCFSNP --> genes_p & genes_pc
    emmax_a & gapit_a & combine_pv --> manh_p
    CLIMATE --> piemap_p
    CLUST --> piemap_p
    prep --> piemap_p
    prep & combine_sel & reg & genes_p --> summary
```

**Red nodes** = DROP path only (per-trait subsetting). Either Path A or Path B runs, not both.

> **GO enrichment** is computed on-demand in the Shiny dashboard (GWAS tab → region detail). Results persist to `GWAS/plots/enrichment/` and `GWAS/tables/enrichment/`.

**Key outputs** (same structure as GEA but under `GWAS/`):
- `GWAS/tables/methods/{method}/` — per-method p-values and significant SNPs
- `GWAS/tables/selected_snps.tsv`, `regions_per_trait.tsv`, `regions_combined.tsv`
- `GWAS/tables/genes_per_region.tsv`, `genes_combined.tsv`
- `GWAS/plots/manhattan/{method}/manhattan_{trait}_K{k}_{adjust}.png/svg`
- `GWAS/plots/manhattan/combined/manhattan_combined_K{k}.png/svg`
- `GWAS/plots/piemap/phenomap_{trait}.png/svg` + `zoom/`

---

### GEAxGWAS Mode — GEA + GWAS overlap analysis

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5

    GEA_SNPS[/"GEA Selected_SNPs"/]:::data
    GWAS_SNPS[/"GWAS Selected_SNPs"/]:::data
    GEA_REG[/"GEA Regions_combined"/]:::data
    GWAS_REG[/"GWAS Regions_combined"/]:::data

    miami["miami_plot\nGEA above / GWAS below\nscattermore background + coords JSON"]
    pairwise["compute_pairwise_overlaps\nAll-vs-all trait SNP overlap"]
    summary["write_summary"]

    GEA_SNPS & GWAS_SNPS --> miami
    GEA_REG & GWAS_REG --> miami
    GEA_SNPS & GWAS_SNPS --> pairwise
    miami & pairwise --> summary
```

> **Region overlap detection, gene annotation, and GO enrichment** for overlapping regions are computed on-demand in the Shiny dashboard (GEAxGWAS tab). The Shiny app dynamically computes regions from the merged GEA+GWAS SNP set with user-configurable overlap bounds (Union / Intersection / GEA-only / GWAS-only).

**Key outputs**:
- `GEAxGWAS/plots/miami_combined_K{k}.png/svg` — static Miami background
- `GEAxGWAS/plots/miami_combined_K{k}_background.png` — scattermore non-sig layer
- `GEAxGWAS/plots/miami_combined_K{k}_coords.json` — sig SNP coordinates for plotly overlay
- `GEAxGWAS/tables/pairwise_collapsed_snps.tsv` — per-trait SNP sets
- `GEAxGWAS/tables/pairwise_overlap_table.tsv` — pairwise overlap scores

---

### Maladaptation Mode — Gradient Forest and genetic offset

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
    combine_gf["combine_gf_snps\nFilter + combine adaptive SNPs\n(Maladaptation combine strategy)"]
    gf_adapt["gradient_forest_adaptive\nAdaptive GF model (sig SNPs)"]
    gf_random["gradient_forest_random\nNeutral GF model (random SNPs)"]:::optional
    offset["gradient_forest_offset\nGenetic offset present → future"]
    cumimp["plot_gf_cumimp\nCumulative importance curves"]
    importance["plot_gf_importance\nOverall R²-weighted importance"]
    piemap["plot_gf_offset_piemap\nOffset PieMap\n{size_trait} ∈ notrait | tajima_d | pi_diversity"]
    density_f["density_plot_future\nFuture climate density"]
    summary["write_summary"]

    META --> dl_future --> merge
    SIG_SNPS --> combine_gf
    combine_gf --> gf_adapt & gf_random
    CLIMATE --> merge & gf_adapt & gf_random & offset
    LFMM_FULL & VCFSNP --> gf_adapt & gf_random
    META --> gf_adapt & gf_random
    gf_adapt --> offset & cumimp & importance
    gf_random --> cumimp & importance
    merge --> offset & density_f
    offset --> piemap
    CLUST --> piemap
    TAJIMA --> piemap
    PI --> piemap
    gf_adapt & offset --> summary
```

**Purple nodes** = optional (random model requires `Maladaptation.methods.gradient_forest.random_model: true`; `tajima_d`/`pi_diversity` piemaps require `Population.calc_stats: true`)

**Key outputs** (under `Maladaptation/{plots,tables}/gradient_forest/{run_label}_{spatial_tag}/`):
- `climate/tables/future/climate_future_year{Y}_ssp{S}_{site,all}.tsv`
- `climate/rasters/future/` — CMIP6 future rasters (.tif)
- `climate/plots/density_plot_future_ssp{S}_{Y}.png/svg`
- `Maladaptation/plots/gradient_forest/{suffix}/cumulative_importance.png/svg`
- `Maladaptation/plots/gradient_forest/{suffix}/overall_importance.png/svg`
- `Maladaptation/plots/gradient_forest/{suffix}/genetic_offset_piemap[_{tajima_d,pi_diversity}].png/svg`
- `Maladaptation/tables/gradient_forest/{suffix}/genetic_offset_{map,site}.tsv`
