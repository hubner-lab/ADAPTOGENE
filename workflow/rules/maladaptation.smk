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
    output: mala_selected_snps('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG)
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
        sigsnps = mala_selected_snps('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG),
        vcfsnp = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = O['metadata'],
        climate = O['climate_site']
    output: mala_model('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, 'adaptive')
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
        sigsnps = mala_selected_snps('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG),
        vcfsnp = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = O['metadata'],
        climate = O['climate_site']
    output: mala_model('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, 'random')
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
        gf = mala_model('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, 'adaptive'),
        future_all = O['climate_future_all'],
        present_all = O['climate_all'],
        present_raster = W['climate_raster'],
        samples = O['metadata']
    output:
        raster = mala_offset_raster('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG),
        map_values = mala_offset_map_values('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG),
        site_values = mala_offset_site_values('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG)
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
        gf = mala_model('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, 'adaptive'),
        gf_random = mala_model('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, 'random') if GF_RANDOM_MODEL else []
    output: mala_cumimp('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG)
    params:
        gf_random_path = mala_model('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, 'random') if GF_RANDOM_MODEL else 'NULL',
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
        gf = mala_model('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, 'adaptive'),
        gf_random = mala_model('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, 'random') if GF_RANDOM_MODEL else []
    output: mala_importance('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG)
    params:
        gf_random_path = mala_model('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, 'random') if GF_RANDOM_MODEL else 'NULL',
        inter_dir = INTER
    log: f"{LOGDIR}maladaptation/plot_gf_importance.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_gf_importance.R \
            {input.gf} {params.gf_random_path} \
            {output} {params.inter_dir} > {log} 2>&1
        """

# Genetic offset PieMaps (collapsed: notrait / tajima_d / pi_diversity)
def _piemap_trait_input(wc):
    """Return trait file path for size-scaled piemaps, or [] for notrait."""
    if wc.size_trait == 'tajima_d':
        return O['tajima']
    elif wc.size_trait == 'pi_diversity':
        return O['pi_div']
    return []

def _piemap_trait_path(wc):
    """Return trait file arg for plot_piemap.R shell command."""
    if wc.size_trait == 'tajima_d':
        return f"{O['tajima']} \"Tajima's D\""
    elif wc.size_trait == 'pi_diversity':
        return f"{O['pi_div']} \"Pi Diversity\""
    return "NULL NULL"

def _piemap_output_prefix(wc):
    """Return OUTPUT_PREFIX arg for plot_piemap.R."""
    return f'genetic_offset_piemap_{wc.size_trait}'

rule plot_gf_offset_piemap:
    """Plot genetic offset piemap, optionally scaled by a population statistic."""
    wildcard_constraints: size_trait = "notrait|tajima_d|pi_diversity"
    input:
        offset_raster = mala_offset_raster('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG),
        samples = O['metadata'],
        clusters = clusters_table(K_BEST),
        trait_file = _piemap_trait_input
    output: mala_offset_piemap('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG, '{size_trait}')
    params:
        pie_alpha = PIEMAP_ALPHA,
        pop_label = PIEMAP_SHOW_LABELS,
        pop_label_size = PIEMAP_LABEL_SIZE,
        pie_scale = PIEMAP_PIE_SCALE,
        use_points = PIEMAP_USE_POINTS,
        plot_dir = mala_plot_dir('gradient_forest', GF_RUN_LABEL, SPATIAL_TAG),
        inter_dir = INTER,
        regionmap_extent = REGIONMAP_EXTENT,
        trait_path = _piemap_trait_path,
        output_prefix = _piemap_output_prefix
    log: f"{LOGDIR}maladaptation/plot_gf_offset_piemap_{{size_trait}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.offset_raster} 1 "Genetic Offset" \
            {input.samples} {input.clusters} \
            {params.trait_path} \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            {params.output_prefix} {params.regionmap_extent} {params.pie_scale} Clusters {params.use_points} > {log} 2>&1
        """
