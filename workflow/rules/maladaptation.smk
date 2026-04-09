#=============================================================================
# MODULE 6: MALADAPTATION
#=============================================================================

# Per-model future climate download (runs in parallel via Snakemake)
rule download_climate_future_model:
    """Download CMIP6 future climate data for a single model."""
    input: samples = O['metadata']
    output: f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_{{model}}.tif"
    wildcard_constraints: model = r"[A-Za-z0-9_-]+"
    params:
        crop = CLIMATE_EXTENT,
        gap = GAP,
        resolution = RESOLUTION,
        data_dir = INDIR
    log: f"{LOGDIR}maladaptation/download_climate_future_{{model}}.log"
    shell:
        """
        Rscript /pipeline/scripts/download_climate_future_model.R \
            {input.samples} {params.crop} {params.gap} {params.data_dir} \
            {SSP} {YEAR} {wildcards.model} {params.resolution} \
            {output} > {log} 2>&1
        """

# Merge per-model rasters into averaged future climate
rule merge_climate_future:
    """Average future climate across models and extract site values."""
    input:
        samples = O['metadata'],
        model_rasters = [f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_{model}.tif" for model in MODELS_LIST],
        present_raster = W['climate_raster'],
        present_all = O['climate_all']
    output:
        raster = W['climate_future_raster'],
        all_vals = O['climate_future_all'],
        site_vals = O['climate_future_site']
    params:
        raster_str = lambda wc, input: ','.join(input.model_rasters),
        n_models = len(MODELS_LIST)
    log: f"{LOGDIR}maladaptation/merge_climate_future.log"
    shell:
        """
        Rscript /pipeline/scripts/merge_climate_future.R \
            {input.samples} {params.raster_str} {params.n_models} \
            {input.present_raster} {input.present_all} \
            {output.raster} {output.all_vals} {output.site_vals} > {log} 2>&1
        """

rule density_plot_future:
    """Generate combined density plot for future climate predictors."""
    input: climate = O['climate_future_site']
    output: O['density_future']
    params:
        predictors = PREDICTORS_SELECTED,
        inter_dir = INTER
    log: f"{LOGDIR}maladaptation/density_plot_future.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_density.R \
            {input.climate} {params.predictors} {output} {params.inter_dir} > {log} 2>&1
        """

# Combine sig SNPs for Gradient Forest using GF-specific strategy
rule combine_gf_snps:
    """Combine significant SNPs for Gradient Forest using the gradient_forest.combine_method strategy."""
    input:
        sigsnps = lambda wc: [assoc_sigsnps(method, adjust) for method, adjust in GEA_CONFIGS.items()]
    output: W['gf_selected_snps']
    params:
        sigsnps_str = lambda wc, input: ' '.join(input.sigsnps),
        method = GF_COMBINE_METHOD,
        gap = GF_COMBINE_GAP,
        predictors = PREDICTORS_SELECTED
    log: f"{LOGDIR}maladaptation/combine_gf_snps.log"
    shell:
        """
        Rscript /pipeline/scripts/combine_selected_snps.R \
            "{params.sigsnps_str}" {params.method} {params.gap} \
            {params.predictors} {output} > {log} 2>&1
        """

# Gradient Forest - adaptive model
rule gradient_forest_adaptive:
    """Build adaptive Gradient Forest model using significant SNPs."""
    input:
        lfmm = W['lfmm_full'],
        sigsnps = W['gf_selected_snps'],
        vcfsnp = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = O['metadata'],
        climate = O['climate_site']
    output: W['gf_adaptive']
    params:
        predictors = PREDICTORS_SELECTED,
        ntree = NTREE,
        cor_threshold = COR_THRESHOLD,
        pcnm = SPATIAL_CORRECTION
    log: f"{LOGDIR}maladaptation/gradient_forest_adaptive.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_model.R \
            {input.lfmm} {input.sigsnps} {input.vcfsnp} {input.removed} \
            {input.samples} {input.climate} {params.predictors} \
            {params.ntree} {params.cor_threshold} {params.pcnm} \
            adaptive {output} > {log} 2>&1
        """

# Gradient Forest - random/neutral model (optional)
rule gradient_forest_random:
    """Build neutral Gradient Forest model using random SNPs."""
    input:
        lfmm = W['lfmm_full'],
        sigsnps = W['gf_selected_snps'],
        vcfsnp = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = O['metadata'],
        climate = O['climate_site']
    output: W['gf_random']
    params:
        predictors = PREDICTORS_SELECTED,
        ntree = NTREE,
        cor_threshold = COR_THRESHOLD,
        pcnm = SPATIAL_CORRECTION
    log: f"{LOGDIR}maladaptation/gradient_forest_random.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_model.R \
            {input.lfmm} {input.sigsnps} {input.vcfsnp} {input.removed} \
            {input.samples} {input.climate} {params.predictors} \
            {params.ntree} {params.cor_threshold} {params.pcnm} \
            random {output} > {log} 2>&1
        """

# Genetic offset calculation
rule gradient_forest_offset:
    """Calculate genetic offset between present and future climate."""
    input:
        gf = W['gf_adaptive'],
        future_all = O['climate_future_all'],
        present_all = O['climate_all'],
        present_raster = W['climate_raster'],
        samples = O['metadata']
    output:
        raster = W['gf_offset_raster'],
        map_values = O['gf_offset_map_values'],
        site_values = O['gf_offset_site_values']
    params:
        predictors = PREDICTORS_SELECTED
    log: f"{LOGDIR}maladaptation/gradient_forest_offset.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_offset.R \
            {input.gf} {params.predictors} {input.future_all} {input.present_all} \
            {input.present_raster} {input.samples} \
            {output.raster} {output.map_values} {output.site_values} > {log} 2>&1
        """

# Cumulative importance plot
rule plot_gf_cumimp:
    """Plot cumulative importance curves for adaptive (and optionally neutral) GF model."""
    input:
        gf = W['gf_adaptive'],
        gf_random = W['gf_random'] if GF_RANDOM_MODEL else []
    output: O['gf_cumimp']
    params:
        gf_random_path = W['gf_random'] if GF_RANDOM_MODEL else 'NULL',
        predictors = PREDICTORS_SELECTED,
        inter_dir = INTER
    log: f"{LOGDIR}maladaptation/plot_gf_cumimp.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_gf_cumimp.R \
            {input.gf} {params.gf_random_path} {params.predictors} \
            {output} {params.inter_dir} > {log} 2>&1
        """

# Overall importance plot
rule plot_gf_importance:
    """Plot R2-weighted importance for adaptive (and optionally neutral) GF model."""
    input:
        gf = W['gf_adaptive'],
        gf_random = W['gf_random'] if GF_RANDOM_MODEL else []
    output: O['gf_importance']
    params:
        gf_random_path = W['gf_random'] if GF_RANDOM_MODEL else 'NULL',
        inter_dir = INTER
    log: f"{LOGDIR}maladaptation/plot_gf_importance.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_gf_importance.R \
            {input.gf} {params.gf_random_path} \
            {output} {params.inter_dir} > {log} 2>&1
        """

# Genetic offset PieMap
rule plot_gf_offset_piemap:
    """Plot genetic offset on map with population structure pie charts (uniform pie size)."""
    input:
        offset_raster = W['gf_offset_raster'],
        samples = O['metadata'],
        clusters = clusters_table(K_BEST)
    output: O['gf_offset_piemap']
    params:
        pie_alpha = PIEMAP_ALPHA,
        pop_label = PIEMAP_SHOW_LABELS,
        pop_label_size = PIEMAP_LABEL_SIZE,
        pie_scale = PIEMAP_PIE_SCALE,
        use_points = PIEMAP_USE_POINTS,
        plot_dir = f"{MOD_MALAD}plots/{GF_RUN_LABEL}_{SPATIAL_TAG}/",
        inter_dir = INTER,
        suffix = f"{GF_RUN_LABEL}_{SPATIAL_TAG}",
        regionmap_extent = REGIONMAP_EXTENT
    log: f"{LOGDIR}maladaptation/plot_gf_offset_piemap.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.offset_raster} 1 "Genetic Offset" \
            {input.samples} {input.clusters} \
            NULL NULL \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            genetic_offset_piemap {params.regionmap_extent} {params.pie_scale} Clusters {params.use_points} > {log} 2>&1
        """

rule plot_gf_offset_piemap_tajima:
    """Plot genetic offset with Tajima's D-scaled pie sizes (requires CALC_POP_STATS=TRUE)."""
    input:
        offset_raster = W['gf_offset_raster'],
        samples = O['metadata'],
        clusters = clusters_table(K_BEST),
        tajima = O['tajima']
    output: O['gf_offset_piemap_tajima']
    params:
        pie_alpha = PIEMAP_ALPHA,
        pop_label = PIEMAP_SHOW_LABELS,
        pop_label_size = PIEMAP_LABEL_SIZE,
        pie_scale = PIEMAP_PIE_SCALE,
        use_points = PIEMAP_USE_POINTS,
        plot_dir = f"{MOD_MALAD}plots/{GF_RUN_LABEL}_{SPATIAL_TAG}/",
        inter_dir = INTER,
        suffix = f"{GF_RUN_LABEL}_{SPATIAL_TAG}",
        regionmap_extent = REGIONMAP_EXTENT
    log: f"{LOGDIR}maladaptation/plot_gf_offset_piemap_tajima.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.offset_raster} 1 "Genetic Offset" \
            {input.samples} {input.clusters} \
            {input.tajima} "Tajima's D" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            genetic_offset_piemap_tajima_d {params.regionmap_extent} {params.pie_scale} Clusters {params.use_points} > {log} 2>&1
        """

rule plot_gf_offset_piemap_diversity:
    """Plot genetic offset with Pi Diversity-scaled pie sizes (requires CALC_POP_STATS=TRUE)."""
    input:
        offset_raster = W['gf_offset_raster'],
        samples = O['metadata'],
        clusters = clusters_table(K_BEST),
        diversity = O['pi_div']
    output: O['gf_offset_piemap_diversity']
    params:
        pie_alpha = PIEMAP_ALPHA,
        pop_label = PIEMAP_SHOW_LABELS,
        pop_label_size = PIEMAP_LABEL_SIZE,
        pie_scale = PIEMAP_PIE_SCALE,
        use_points = PIEMAP_USE_POINTS,
        plot_dir = f"{MOD_MALAD}plots/{GF_RUN_LABEL}_{SPATIAL_TAG}/",
        inter_dir = INTER,
        suffix = f"{GF_RUN_LABEL}_{SPATIAL_TAG}",
        regionmap_extent = REGIONMAP_EXTENT
    log: f"{LOGDIR}maladaptation/plot_gf_offset_piemap_diversity.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.offset_raster} 1 "Genetic Offset" \
            {input.samples} {input.clusters} \
            {input.diversity} "Pi Diversity" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            genetic_offset_piemap_pi_diversity {params.regionmap_extent} {params.pie_scale} Clusters {params.use_points} > {log} 2>&1
        """
