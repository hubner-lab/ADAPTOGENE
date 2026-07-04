#=============================================================================
# MODULE 6: MALADAPTATION
#=============================================================================

# Wildcard constraints for all Gradient Forest rules.
# run_label = set name (matches Shiny charset [A-Za-z0-9_.-]+, no slashes).
# spatial_tag = exactly 'spatial' or 'nospatial' — disambiguates trailing token
#   even when run_label itself contains underscores (e.g. 'EMMAX_bonf005').
wildcard_constraints:
    run_label   = r"[A-Za-z0-9_.-]+",
    spatial_tag = r"spatial|nospatial",
    method      = r"gradient_forest|geometric_offset"

if CLIMATE_SOURCE == 'worldclim':
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
else:
    rule stage_custom_climate_future:
        """Stage user-supplied future climate table into pipeline-standard outputs."""
        input:
            samples        = O['metadata'],
            present_all    = O['climate_all'],
            present_raster = W['climate_raster']
        output:
            raster    = W['climate_future_raster'],
            all_vals  = O['climate_future_all'],
            site_vals = O['climate_future_site']
        params:
            env_table = CUSTOM_FUTURE_TABLE,
            columns   = CUSTOM_CLIMATE_COLUMNS,
            key       = CUSTOM_CLIMATE_KEY
        log: f"{LOGDIR}maladaptation/stage_custom_climate_future.log"
        shell:
            """
            Rscript /pipeline/scripts/stage_custom_climate.R future \
                {input.samples} {params.env_table} {params.columns} {params.key} \
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

# Gradient Forest - adaptive model
# sigsnps is an ancestor-less source file produced by the Shiny GEA tab.
# {run_label} = saved SNP-set name; {spatial_tag} = spatial|nospatial.
# pcnm param translates spatial_tag -> 'with'/'without' for the R script.
rule gradient_forest_adaptive:
    """Build adaptive Gradient Forest model using a user-curated SNP set."""
    input:
        lfmm    = W['lfmm_full'],
        sigsnps = lambda wc: snp_set_file(wc.run_label),
        vcfsnp  = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = O['metadata'],
        climate = O['climate_site']
    output: mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'adaptive')
    params:
        predictors    = PREDICTORS_SELECTED,
        ntree         = NTREE,
        cor_threshold = COR_THRESHOLD,
        pcnm          = lambda wc: 'with' if wc.spatial_tag == 'spatial' else 'without'
    log: f"{LOGDIR}maladaptation/gradient_forest_adaptive_{{run_label}}_{{spatial_tag}}.log"
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
        lfmm    = W['lfmm_full'],
        sigsnps = lambda wc: snp_set_file(wc.run_label),
        vcfsnp  = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = O['metadata'],
        climate = O['climate_site']
    output: mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'random')
    params:
        predictors    = PREDICTORS_SELECTED,
        ntree         = NTREE,
        cor_threshold = COR_THRESHOLD,
        pcnm          = lambda wc: 'with' if wc.spatial_tag == 'spatial' else 'without'
    log: f"{LOGDIR}maladaptation/gradient_forest_random_{{run_label}}_{{spatial_tag}}.log"
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
        gf             = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'adaptive'),
        future_all     = O['climate_future_all'],
        present_all    = O['climate_all'],
        present_raster = W['climate_raster'],
        samples        = O['metadata']
    output:
        raster      = mala_offset_raster('gradient_forest', '{run_label}', '{spatial_tag}'),
        map_values  = mala_offset_map_values('gradient_forest', '{run_label}', '{spatial_tag}'),
        site_values = mala_offset_site_values('gradient_forest', '{run_label}', '{spatial_tag}')
    params:
        predictors = PREDICTORS_SELECTED
    log: f"{LOGDIR}maladaptation/gradient_forest_offset_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_offset.R \
            {input.gf} {params.predictors} {input.future_all} {input.present_all} \
            {input.present_raster} {input.samples} \
            {output.raster} {output.map_values} {output.site_values} > {log} 2>&1
        """

# Geometric Genetic Offset (Gain et al. 2023, MBE) — single-call rule.
# No separate model artifact: genetic.gap() fits LFMM2 + computes offset in one pass.
# Nospatial-only: genetic.gap() has no spatial-correction argument.
# candidate.loci = INTEGER INDEX vector into the full imputed matrix (never subset the matrix).
rule geometric_offset:
    """Compute geometric genetic offset using LEA::genetic.gap()."""
    input:
        lfmm_full   = W['lfmm_imp_full'],
        vcfsnp      = W['vcfsnp_full'],
        removed     = W['removed_full'],
        sigsnps     = lambda wc: snp_set_file(wc.run_label),
        env_site_pres  = O['climate_site'],
        env_site_fut   = O['climate_future_site'],
        env_all_pres   = O['climate_all'],
        env_all_fut    = O['climate_future_all'],
        pres_raster    = W['climate_raster'],
        samples        = O['metadata']
    output:
        site_values = mala_offset_site_values('geometric_offset', '{run_label}', '{spatial_tag}'),
        map_values  = mala_offset_map_values('geometric_offset', '{run_label}', '{spatial_tag}'),
        raster      = mala_offset_raster('geometric_offset', '{run_label}', '{spatial_tag}'),
        importance  = mala_importance('geometric_offset', '{run_label}', '{spatial_tag}')
    wildcard_constraints:
        spatial_tag = r"nospatial"   # geometric_offset is nospatial-only
    params:
        predictors = PREDICTORS_SELECTED,
        k          = lambda wc: GO_K if GO_K != '' else K_BEST,
        scale      = GO_SCALE
    log: f"{LOGDIR}maladaptation/geometric_offset_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/geometric_offset.R \
            {input.lfmm_full} {input.vcfsnp} {input.removed} {input.sigsnps} \
            {input.env_site_pres} {input.env_site_fut} \
            {input.env_all_pres} {input.env_all_fut} \
            {input.pres_raster} {input.samples} \
            {params.predictors} {params.k} {params.scale} \
            {output.site_values} {output.map_values} \
            {output.raster} {output.importance} > {log} 2>&1
        """

# Cumulative importance plot
rule plot_gf_cumimp:
    """Plot cumulative importance curves for adaptive (and optionally neutral) GF model."""
    input:
        gf        = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'adaptive'),
        gf_random = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'random') if GF_RANDOM_MODEL else []
    output: mala_cumimp('gradient_forest', '{run_label}', '{spatial_tag}')
    params:
        gf_random_path = lambda wc: mala_model('gradient_forest', wc.run_label, wc.spatial_tag, 'random') if GF_RANDOM_MODEL else 'NULL',
        predictors     = PREDICTORS_SELECTED,
        inter_dir      = INTER
    log: f"{LOGDIR}maladaptation/plot_gf_cumimp_{{run_label}}_{{spatial_tag}}.log"
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
        gf        = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'adaptive'),
        gf_random = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'random') if GF_RANDOM_MODEL else []
    output: mala_importance('gradient_forest', '{run_label}', '{spatial_tag}')
    params:
        gf_random_path = lambda wc: mala_model('gradient_forest', wc.run_label, wc.spatial_tag, 'random') if GF_RANDOM_MODEL else 'NULL',
        inter_dir      = INTER
    log: f"{LOGDIR}maladaptation/plot_gf_importance_{{run_label}}_{{spatial_tag}}.log"
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
    """Plot genetic offset piemap for any maladaptation method, optionally scaled by a population statistic."""
    wildcard_constraints: size_trait = "notrait|tajima_d|pi_diversity"
    input:
        offset_raster = lambda wc: mala_offset_raster(wc.method, wc.run_label, wc.spatial_tag),
        samples       = O['metadata'],
        clusters      = clusters_table(K_BEST),
        trait_file    = _piemap_trait_input
    output: mala_offset_piemap('{method}', '{run_label}', '{spatial_tag}', '{size_trait}')
    params:
        pie_alpha       = PIEMAP_ALPHA,
        pop_label       = PIEMAP_SHOW_LABELS,
        pop_label_size  = PIEMAP_LABEL_SIZE,
        pie_scale       = PIEMAP_PIE_SCALE,
        use_points      = PIEMAP_USE_POINTS,
        plot_dir        = lambda wc: mala_plot_dir(wc.method, wc.run_label, wc.spatial_tag),
        inter_dir       = INTER,
        regionmap_extent = REGIONMAP_EXTENT,
        trait_path      = _piemap_trait_path,
        output_prefix   = _piemap_output_prefix
    log: f"{LOGDIR}maladaptation/plot_offset_piemap_{{method}}_{{run_label}}_{{spatial_tag}}_{{size_trait}}.log"
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
