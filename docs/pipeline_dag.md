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

    subgraph proc["1 · Processing"]
        direction LR
        P["Filter + LD prune + normalize"]
    end

    subgraph struct["2 · Structure"]
        direction LR
        S["sNMF + PCA + cross-entropy"]
    end

    subgraph structK["3 · Structure K"]
        direction LR
        SK["Impute + climate + piemaps + pop stats"]
    end

    subgraph assoc["4 · Association — GEA"]
        direction LR
        A["EMMAX/LFMM → regions → genes → GO"]
    end

    subgraph pheno["5 · Phenotype Association"]
        direction LR
        AP["EMMAX per trait → regions → genes → GO"]
    end

    subgraph overlap_mode["6 · Overlapping"]
        direction LR
        OV["GEA∩GWAS → new regions → genes → GO"]
    end

    subgraph regionplot_mode["7 · Regionplot"]
        direction LR
        RP["Regional Manhattan + gene annotations"]
    end

    subgraph malad["8 · Maladaptation"]
        direction LR
        M["Gradient Forest → genetic offset"]
    end

    subgraph hap_scan["9 · Haplotype Scan"]
        direction LR
        HS["crosshap epsilon scan → clustree"]
    end

    subgraph hap_viz["10 · Haplotype"]
        direction LR
        HV["crosshap_viz + boxplots + piemaps"]
    end

    inputs --> proc --> struct
    struct -. "Select K_BEST" .-> structK
    structK --> assoc & pheno
    assoc --> regionplot_mode & malad
    assoc & pheno --> overlap_mode
    assoc & pheno & overlap_mode -.-> hap_scan
    hap_scan -. "Select epsilon" .-> hap_viz
```

Each mode is run separately via `--config mode=<MODE>`.

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
- `processing/tables/metadata.tsv` — aligned metadata (sample order matches VCF)
- `processing/tables/sample_missing_stats.tsv` — per-sample missingness
- `_intermediate/annotation/normalized.gff3` — chr-stripped GFF

---

### Structure Mode — Population structure inference

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

**Key outputs**:
- `structure/plots/pca.png/svg`, `tracy_widom.png/svg` — PCA scatter + eigenvalue test
- `structure/plots/cross_entropy_K{start}-{end}.png/svg` — K selection guide
- `structure/plots/structure_K{k}.png/svg` — ancestry barplots (per K)
- `structure/plots/pca_structure_K{k}.png/svg` — PCA with ancestry pies
- `structure/plots/pop_diff_K{k}.png/svg` — population differentiation
- `structure/tables/clusters_K{k}.tsv` — Q-matrices per K

---

### Structure K Mode — Imputation, climate, population statistics

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

**Purple nodes** = optional, require `pop.calc_stats: TRUE`

**Key outputs**:
- `climate/rasters/present/climate_present_rasterstack.grd` — WorldClim bioclim rasters
- `climate/tables/present/climate_present_{all,site,site_scaled}.tsv`
- `climate/plots/density_plot_present.png/svg`, `correlation_heatmap.png/svg`
- `structure_k/plots/piemap/piemap_{bio}.png/svg/qs` + `zoom/` variants
- `structure_k/plots/piemap/{tajima_d,pi_diversity}/piemap_{bio}_{metric}.png/svg/qs` (optional)
- `structure_k/tables/pop_stats/tajima_d_by_pop.tsv`, `pi_diversity_by_pop.tsv`, `ibd_raw.tsv`, `ibd_pairs.tsv`, `amova.tsv`
- `structure_k/plots/pop_stats/mantel_test.png/svg`, `amova.png/svg`

---

### Association Mode (GEA) — Climate association analysis

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

**Key outputs**:
- `association/tables/{METHOD}/{method}_pvalues_K{k}.tsv` — per-method p-values
- `association/tables/{METHOD}/{method}_pvalues_K{k}_sig_snps_{adjust}.tsv` — significant SNPs
- `association/tables/selected_snps.tsv` — merged across methods
- `association/tables/regions_per_trait.tsv`, `regions_combined.tsv` — clustered regions
- `association/tables/genes_per_region.tsv`, `genes_per_region_collapsed.tsv`, `genes_combined.tsv`
- `association/plots/manhattan/{METHOD}/manhattan_{trait}_K{k}_{adjust}.png/svg` (plain + regions)
- `association/plots/manhattan_combined_K{k}.png/svg`
- `association/plots/enrichment/{trait}/region_{id}_dotplot.png/svg` (+ emapplot, cnetplot)
- `association/tables/enrichment/{trait}/region_{id}_enrichment.tsv`

---

### Phenotype Association Mode — GWAS on phenotypic traits

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

**Red nodes** = DROP path only (per-trait subsetting). Either Path A or Path B runs, not both.

**Key outputs**: Same structure as association/ but under `phenotype_association/`. Additionally:
- `phenotype_association/tables/phenotype_missing_summary.tsv`
- `phenotype_association/plots/piemap/phenomap_{trait}.png/svg/qs` + `zoom/`

---

### Overlapping Mode — GEA + GWAS combined

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

**Key outputs**:
- `overlapping/tables/selected_snps_all.tsv` — merged GEA+GWAS SNPs
- `overlapping/tables/regions_per_trait_all.tsv`, `regions_combined.tsv`
- `overlapping/tables/overlap_summary.tsv` — overlap statistics
- `overlapping/tables/genes_per_region.tsv`, `genes_per_region_collapsed.tsv`, `genes_combined.tsv`
- `overlapping/plots/enrichment/{trait}/region_{id}_dotplot.png/svg`

---

### Regionplot Mode — Regional Manhattan with gene annotations

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

**Key outputs**:
- `regionplot/plots/regionplot_{region}_{trait}.png` + `_main.svg`, `_overview.svg`, `_genes.svg`

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
    piemap["plot_gf_offset_piemap\nOffset PieMap\n{size_trait} ∈ notrait | tajima_d | pi_diversity"]
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
    offset --> piemap
    CLUST --> piemap
    TAJIMA --> piemap
    PI --> piemap
    gf_adapt & offset --> summary
```

**Purple nodes** = optional (random model requires `gradient_forest.random_model: TRUE`; `tajima_d`/`pi_diversity` piemaps require `pop.calc_stats: TRUE`)

**Key outputs** (under `Maladaptation/{plots,tables}/{method}/{run_label}_{spatial_tag}/`):
- `climate/tables/future/climate_future_year{Y}_ssp{S}_{site,all}.tsv`
- `climate/rasters/future/` — CMIP6 future rasters
- `climate/plots/density_plot_future_ssp{S}_{Y}.png/svg`
- `Maladaptation/plots/gradient_forest/{suffix}/cumulative_importance.png/svg`
- `Maladaptation/plots/gradient_forest/{suffix}/overall_importance.png/svg`
- `Maladaptation/plots/gradient_forest/{suffix}/genetic_offset_piemap[_{tajima_d,pi_diversity}].png/svg`
- `Maladaptation/tables/gradient_forest/{suffix}/genetic_offset_{map,site}.tsv`

---

### Haplotype Scan Mode — crosshap epsilon scanning

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5

    REGIONS[/"Combined Regions\n(assoc/pheno/overlap)"/]:::data
    META[/"Aligned Metadata"/]:::data
    VCF_FILT[/"Filtered VCF"/]:::data
    VCF_IMP[/"Imputed VCF (full)"/]:::data
    CLUST[/"Clusters (K_BEST)"/]:::data

    select["haplotype_select_regions\nTop N regions + suggested MGmin"]
    scan["haplotype_scan\nbcftools region extract\nplink LD matrix\ncrosshap run_haplotyping\nclustreeviz"]
    summary["write_summary"]

    REGIONS --> select
    META --> select & scan
    select --> scan
    VCF_FILT --> scan
    VCF_IMP --> scan
    CLUST -.-> scan
    scan --> summary
```

**Key outputs**:
- `haplotype_scan/{tag}/tables/selected_regions.tsv` — top regions for scanning
- `haplotype_scan/{tag}/tables/scan_status.tsv` — per-region scan status
- `haplotype_scan/{tag}/plots/clustree/region_{id}_clustree.png/svg`
- `_intermediate/haplotype/{tag}/region_{id}_hapobject.qs` — crosshap HapObjects

---

### Haplotype Mode — Visualization at selected epsilon

```mermaid
flowchart TB
    classDef data fill:#FFFDE7,stroke:#F9A825,color:#333,stroke-dasharray:5 5

    SCAN[/"HapObjects (from scan)"/]:::data
    REGIONS[/"Selected Regions"/]:::data
    META[/"Aligned Metadata"/]:::data
    RASTER[/"Climate Raster"/]:::data
    CLUST[/"Clusters (K_BEST)"/]:::data

    viz["haplotype_viz\ncrosshap_viz at selected epsilon\nphenotype boxplots (ggdist + ggpubr)\nhaplotype piemaps (plot_piemap.R)\nfrequency tables"]
    summary["write_summary"]

    SCAN --> viz
    REGIONS --> viz
    META --> viz
    RASTER --> viz
    CLUST --> viz
    viz --> summary
```

**Key outputs**:
- `haplotype/{tag}/plots/region_{id}_crosshap_viz.png/svg`
- `haplotype/{tag}/plots/region_{id}_boxplot_{trait}.png/svg`
- `haplotype/{tag}/plots/region_{id}_piemap_{trait}.png/svg/qs` + `zoom/`
- `haplotype/{tag}/tables/region_{id}_assignments.tsv`, `_frequencies.tsv`
