#=============================================================================
# MODULE 3: STRUCTURE_K (requires K_BEST selection)
#=============================================================================

rule impute_ld:
    """Impute missing genotypes using SNMF results for K_BEST."""
    input:  snmf = W['snmf'], lfmm = W['lfmm']
    output: W['lfmm_imp']
    params: k = K_BEST
    log:    f"{LOGDIR}structure_k/impute_ld.log"
    shell:
        """
        Rscript /pipeline/scripts/impute.R {input.snmf} {input.lfmm} {params.k} {output} > {log} 2>&1
        """

rule lfmm2vcf_ld:
    """Convert imputed LFMM back to VCF format."""
    input:  vcf = W['vcf_ld'], lfmm_imp = W['lfmm_imp']
    output: W['vcf_imp']
    params: blocksize = 20000
    log:    f"{LOGDIR}structure_k/lfmm2vcf_ld.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm2vcf.R {input.lfmm_imp} {input.vcf} {params.blocksize} > {log} 2>&1
        """
    # NOTE: lfmm2vcf.R already exists, outputs .vcf with same base as .lfmm input

rule download_climate_present:
    """Download and process present climate data for sampling locations."""
    input:  meta = O['metadata']
    output:
        site = O['climate_site'],
        site_scaled = O['climate_site_scaled'],
        all_vals = O['climate_all'],
        raster = W['climate_raster']
    params:
        crop = CROP_REGION,
        gap = GAP,
        resolution = RESOLUTION,
        data_dir = f"{INDIR}",
        raster_dir = f"{MOD_CLIMATE}rasters/present/",
        tables_dir = f"{MOD_CLIMATE}tables/present/"
    log: f"{LOGDIR}structure_k/download_climate_present.log"
    shell:
        """
        Rscript /pipeline/scripts/download_climate_present.R \
            {input.meta} {params.crop} {params.gap} {params.data_dir} {params.resolution} \
            {params.raster_dir} {params.tables_dir} > {log} 2>&1
        """

rule density_plot:
    """Generate combined density plot for all climate predictors."""
    input:  climate = O['climate_site']
    output: DENSITY_PLOT_COMBINED
    params:
        predictors = PREDICTORS_SELECTED,
        inter_dir = INTER
    log:    f"{LOGDIR}structure_k/density_plot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_density.R \
            {input.climate} {params.predictors} {output} {params.inter_dir} > {log} 2>&1
        """

rule tajima_d:
    """Calculate Tajima's D per population."""
    input:  vcf = W['vcf_filt'], meta = O['metadata']
    output: O['tajima']
    params: winsize = METRICS_WINSIZE, inter_dir = INTER
    log:    f"{LOGDIR}structure_k/tajima_d.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/tajima_d.R \
            {input.vcf} {input.meta} {params.winsize} {threads} \
            {output} {params.inter_dir} > {log} 2>&1
        """

rule pi_diversity:
    """Calculate nucleotide diversity (Pi) per population."""
    input:  vcf = W['vcf_filt'], meta = O['metadata']
    output: O['pi_div']
    params: winsize = METRICS_WINSIZE, inter_dir = INTER
    log:    f"{LOGDIR}structure_k/pi_diversity.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/pi_diversity.R \
            {input.vcf} {input.meta} {params.winsize} {threads} \
            {output} {params.inter_dir} > {log} 2>&1
        """

rule ibd:
    """Calculate isolation by distance between populations."""
    input:  clusters = clusters_table(K_BEST), meta = O['metadata']
    output: raw = O['ibd_raw'], pairs = O['ibd_pairs']
    params: tables_dir = f"{MOD_STRUCTK}tables/pop_stats/"
    log:    f"{LOGDIR}structure_k/ibd.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/ibd.R \
            {input.clusters} {input.meta} {threads} {params.tables_dir} > {log} 2>&1
        """

rule correlation_heatmap:
    """Generate correlation heatmap of climate variables and traits."""
    input:  climate = O['climate_site'], meta = O['metadata']
    output: O['corr_heatmap']
    params: inter_dir = INTER
    log:    f"{LOGDIR}structure_k/correlation_heatmap.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_correlation_heatmap.R \
            {input.climate} {input.meta} {output} {params.inter_dir} > {log} 2>&1
        """

rule mantel_test:
    """Perform Mantel test for IBD/IBE."""
    input:  meta = O['metadata'], climate = O['climate_site'], clusters = clusters_table(K_BEST)
    output: O['mantel']
    params:
        predictors = PREDICTORS_SELECTED,
        plot_dir = f"{MOD_STRUCTK}plots/pop_stats/",
        inter_dir = INTER
    log: f"{LOGDIR}structure_k/mantel_test.log"
    shell:
        """
        Rscript /pipeline/scripts/mantel_test.R \
            {input.meta} {input.clusters} {input.climate} {params.predictors} \
            {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule amova:
    """Perform AMOVA analysis."""
    input:  vcf = W['vcf_imp'], meta = O['metadata']
    output: table = O['amova'], plot = O['amova_plot']
    params: plot_dir = f"{MOD_STRUCTK}plots/pop_stats/", tables_dir = f"{MOD_STRUCTK}tables/pop_stats/", inter_dir = INTER
    log:    f"{LOGDIR}structure_k/amova.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/amova.R \
            {input.vcf} {input.meta} {threads} \
            {params.plot_dir} {params.tables_dir} {params.inter_dir} > {log} 2>&1
        """

rule piemap_plot:
    """Generate PieMap plots for a climate predictor with trait-scaled pie sizes."""
    input:
        meta = O['metadata'],
        clusters = clusters_table(K_BEST),
        raster = W['climate_raster'],
        tajima = O['tajima'],
        diversity = O['pi_div']
    output:
        tajima_plot = piemap_tajima("{bio}"),
        diversity_plot = piemap_diversity("{bio}")
    wildcard_constraints: bio = r"bio_\d+"
    params:
        bio = lambda wc: wc.bio,
        pie_alpha = PIEMAP_ALPHA,
        pop_label = PIEMAP_SHOW_LABELS,
        pop_label_size = PIEMAP_LABEL_SIZE,
        pie_scale = PIEMAP_PIE_SCALE,
        use_points = PIEMAP_USE_POINTS,
        plot_dir = f"{MOD_STRUCTK}plots/piemap/",
        inter_dir = INTER,
        regionmap_extent = REGIONMAP_EXTENT
    log: f"{LOGDIR}structure_k/piemap_{{bio}}.log"
    shell:
        """
        # TajimaD piemap
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            {input.tajima} "Tajima's D" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            piemap_{params.bio}_tajima_d {params.regionmap_extent} {params.pie_scale} Clusters {params.use_points} > {log} 2>&1

        # PiDiversity piemap
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            {input.diversity} "Pi Diversity" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            piemap_{params.bio}_pi_diversity {params.regionmap_extent} {params.pie_scale} Clusters {params.use_points} >> {log} 2>&1
        """

rule piemap_simple:
    """Generate basic PieMap plots with autoscaled uniform pie size (no trait overlay)."""
    input:
        meta = O['metadata'],
        clusters = clusters_table(K_BEST),
        raster = W['climate_raster']
    output: piemap_notrait("{bio}")
    wildcard_constraints: bio = r"bio_\d+"
    params:
        bio = lambda wc: wc.bio,
        pie_alpha = PIEMAP_ALPHA,
        pop_label = PIEMAP_SHOW_LABELS,
        pop_label_size = PIEMAP_LABEL_SIZE,
        pie_scale = PIEMAP_PIE_SCALE,
        use_points = PIEMAP_USE_POINTS,
        plot_dir = f"{MOD_STRUCTK}plots/piemap/",
        inter_dir = INTER,
        regionmap_extent = REGIONMAP_EXTENT
    log: f"{LOGDIR}structure_k/piemap_simple_{{bio}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            NULL NULL \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            piemap_{params.bio} {params.regionmap_extent} {params.pie_scale} Clusters {params.use_points} > {log} 2>&1
        """
