# MODULE: CLIMATE — predictor characterization (mode=climate)
#
# Owns the full predictor story: correlation/density plots, invariant-
# predictor detection, and spatial structure (dbMEM + variance partitioning).
# Split out of the old preGEA module (which mixed predictor characterization
# with LFMM-K/EMMAX-#PC/RDA hyperparameter sweeps) so each module answers one
# question: Climate answers "what do my predictors look like and how much do
# they overlap with geography/structure", PreGEA answers "how many K / #PCs".
#
# Climate DATA ACQUISITION (WorldClim download / custom-table staging) stays
# in structure.smk — climate_site is an input to Structure's own piemaps, and
# the boundary here is "predictor characterization", not "data acquisition".
#
# Two blocks:
#   Block 1 — predictor characterization (moved from structure.smk unchanged)
#   Block 2 — dbMEM (shared spatial artifact) + varpart (moved from pregea.smk,
#             renamed climate_dbmem/climate_varpart)

# =============================================================================
# Block 1 — predictor characterization
# =============================================================================

rule check_climate_variance:
    """Detect invariant (zero-variance or all-NA) bioclimatic predictors at sample sites."""
    input:  site = O['climate_site']
    output: O['climate_invariant']
    log:    f"{LOGDIR}climate/check_climate_variance.log"
    shell:
        """
        Rscript /pipeline/scripts/check_climate_variance.R \
            {input.site} {output} > {log} 2>&1
        """

rule density_plot:
    """Generate combined density plot for all climate predictors (all BIO columns)."""
    input:  climate = O['climate_site']
    output: DENSITY_PLOT_COMBINED
    params:
        predictors = "all",
        inter_dir = INTER
    log:    f"{LOGDIR}climate/density_plot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_density.R \
            {input.climate} {params.predictors} {output} {params.inter_dir} > {log} 2>&1
        """

if META_HAS_PHENO:
    rule density_plot_phenotypes:
        """Generate combined density plot for phenotype traits (metadata columns 5+)."""
        input:  meta = O['metadata']
        output: DENSITY_PLOT_PHENOTYPES
        params:
            predictors = "all",
            inter_dir = INTER
        log:    f"{LOGDIR}climate/density_plot_phenotypes.log"
        shell:
            """
            Rscript /pipeline/scripts/plot_density.R \
                {input.meta} {params.predictors} {output} {params.inter_dir} > {log} 2>&1
            """

rule correlation_heatmap:
    """Generate correlation heatmap of climate variables and traits.
    Uses metadata_climate_valid.tsv (climate-valid samples) -- plot_correlation_heatmap.R
    cbind()s climate values to trait columns positionally, so meta must have the same row
    count/order as climate (bonus fix: was O['metadata'], the full cohort, which silently
    misaligned whenever any sample had missing coordinates -- a pre-existing bug independent
    of the climate-NA warn-and-exclude feature this rule is now consistent with)."""
    input:  climate = O['climate_site'], meta = W['metadata_climate_valid']
    output: O['corr_heatmap']
    params: inter_dir = INTER
    log:    f"{LOGDIR}climate/correlation_heatmap.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_correlation_heatmap.R \
            {input.climate} {input.meta} {output} {params.inter_dir} > {log} 2>&1
        """

# =============================================================================
# Block 2 — dbMEM (shared spatial artifact) + varpart. Split into two rules
# specifically so the shared artifact is cheap and standalone: climate_dbmem
# depends only on coordinates, so Gradient Forest's 'spatial' variant can
# request it in maladaptation mode without dragging any genotype work into
# the DAG.
# =============================================================================

rule climate_dbmem:
    """SHARED spatial-vector artifact — the single dbMEM computation reused
    by varpart AND Gradient Forest (A20/C4). Replaces gradient_forest_model.R's
    old ad hoc raw-Euclidean-dist + "keep half the positive eigenvalues"
    heuristic — see scripts/pregea_dbmem.R header."""
    input:
        metadata = W['metadata_climate_valid'],
        samples  = W['climate_valid_samples'],
    output:
        vectors = O['climate_dbmem_vectors'], diag = O['climate_dbmem_diag'],
        png = O['climate_dbmem_png'], svg = O['climate_dbmem_svg'],
    params: level = CLIMATE_DBMEM_LEVEL, min_sites = CLIMATE_DBMEM_MIN_SITES, inter_dir = INTER
    log:    f"{LOGDIR}climate/dbmem.log"
    shell:
        """
        Rscript /pipeline/scripts/pregea_dbmem.R \
            {input.metadata} {input.samples} {params.level} {params.min_sites} \
            {output.vectors} {output.diag} {output.png} {params.inter_dir} > {log} 2>&1
        """

rule climate_varpart:
    """Shared climate/structure/geography variance-partitioning panel.
    scripts/pregea_varpart.R — see its header for the full port-and-correction
    notes from rda_varpart_lasky.Rmd."""
    input:
        projections = W['pca_projections'],
        eigenvalues = W['pca_eigenvalues'],
        dbmem       = O['climate_dbmem_vectors'],
        clusters    = clusters_table(K_BEST) if CLIMATE_VP_STRUCT == 'qmatrix' else [],
        climate     = O['climate_site_scaled'],
        climate_valid = W['climate_valid_samples'],
        samples_order = W['samples_order'],
        lfmm_pruned = W['lfmm_imp_climate'] if CLIMATE_VP_RESPONSE == 'snps' else [],
    output:
        selection  = O['climate_vp_selection'], selected = O['climate_vp_selected'],
        fractions  = O['climate_vp_fractions'], anova    = O['climate_vp_anova'],
        px         = O['climate_vp_px'],
        p_path = O['climate_vp_path_png'], p_venn = O['climate_vp_venn_png'],
        p_frac = O['climate_vp_frac_png'], p_px   = O['climate_vp_px_png'],
    params:
        clusters_arg = lambda wc, input: input.clusters if CLIMATE_VP_STRUCT == 'qmatrix' else 'NULL',
        lfmm_arg     = lambda wc, input: input.lfmm_pruned if CLIMATE_VP_RESPONSE == 'snps' else 'NULL',
        predictors = PREDICTORS_SELECTED,
        response = CLIMATE_VP_RESPONSE, response_var = CLIMATE_VP_VAR_CUTOFF,
        response_max = CLIMATE_VP_MAX_PCS, response_min = CLIMATE_VP_MIN_PCS,
        pin = CLIMATE_VP_PIN, r2perms = CLIMATE_VP_R2PERMS,
        varpartperms = CLIMATE_VP_PERMS, confound = "TRUE" if CLIMATE_VP_CONFOUND else "FALSE",
        seed = PREGEA_SEED, plot_dir = f"{MOD_CLIMATE}plots/varpart/", inter_dir = INTER,
    log: f"{LOGDIR}climate/varpart.log"
    shell:
        """
        Rscript /pipeline/scripts/pregea_varpart.R \
            {input.projections} {input.eigenvalues} {input.dbmem} {params.clusters_arg} \
            {input.climate} {input.climate_valid} {input.samples_order} {params.predictors} \
            {params.response} {params.response_var} {params.response_max} {params.response_min} \
            {params.lfmm_arg} {params.pin} {params.r2perms} {params.varpartperms} \
            {params.confound} {params.seed} \
            {output.selection} {output.selected} {output.fractions} {output.anova} {output.px} \
            {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """
