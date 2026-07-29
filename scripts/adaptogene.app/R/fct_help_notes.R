#' Central registry of plot help-notes: what a plot shows + which config controls it.
#'
#' One entry per note-ready plot card, feeding help_note() below. Centralizing this
#' prose here (instead of inlining a description + config path at ~24 module call
#' sites) keeps it a single source of truth, per the CLAUDE.md "compute once, reuse
#' downstream" / "single source of truth" rules. Each entry:
#'   - desc:   one-sentence description of what the plot shows
#'   - config: config dot-path(s) that control it, or NULL if none applies
#' @noRd
HELP_NOTES <- list(
    snp_density = list(
        desc   = "SNPs retained per chromosome after all filters.",
        config = NULL
    ),
    pca = list(
        desc   = "Sample scatter on the first two genotype PCs, labeled by sampling site — unsupervised structure overview before any clustering model is fit.",
        config = NULL
    ),
    tracy_widom = list(
        desc   = "Tracy-Widom test statistic per PC — the eigenvalue variance each PC explains, in decreasing order. PCs before the curve flattens capture real structure; flat PCs are noise.",
        config = NULL
    ),
    cross_entropy = list(
        desc   = "sNMF cross-entropy per tested K (lower = better fit) — the standard criterion for choosing the number of ancestral populations.",
        config = "sNMF.k_start, sNMF.k_end, sNMF.repeats"
    ),
    structure_bar = list(
        desc   = "sNMF ancestry proportions per sample at the selected K, one bar per sample.",
        config = "sNMF.k_best"
    ),
    pca_structure = list(
        desc   = "PCA scatter with each sample's genotype PCs pie-chart colored by its sNMF ancestry proportions at the selected K.",
        config = "sNMF.k_best"
    ),
    pop_diff = list(
        desc   = "sNMF population-differentiation p-values per SNP among K clusters (histogram + genome-wide -log10(p) track).",
        config = "sNMF.k_best"
    ),
    pca_structure_k = list(
        desc   = "PCA colored by cluster at the best K.",
        config = "sNMF.k_best"
    ),
    climate_heatmap = list(
        desc   = "Pairwise correlation among climate predictors — collinearity check before GEA.",
        config = "Climate.predictors"
    ),
    climate_density = list(
        desc   = "Distribution of present-day climate values across sampling sites.",
        config = "Climate.predictors, Map.climate_extent"
    ),
    phenotype_density = list(
        desc   = "Distribution of phenotype/trait values across samples.",
        config = NULL
    ),
    mantel = list(
        desc   = "Isolation-by-distance: genetic vs geographic distance correlation.",
        config = "Population.calc_stats, Population.window_size"
    ),
    amova = list(
        desc   = "Variance partition among vs within populations.",
        config = "Population.calc_stats"
    ),
    qq_plot_gea = list(
        desc   = "Observed vs expected p-value quantiles — calibration check for this method/trait.",
        config = "GEA.configs"
    ),
    qq_plot_gwas = list(
        desc   = "Observed vs expected p-value quantiles — calibration check for this method/trait.",
        config = "GWAS.configs"
    ),
    enrich_dotplot = list(
        desc   = "Top enriched GO terms for this region's genes.",
        config = "Enrichment.top_terms, GFF.go_field"
    ),
    enrich_emapplot = list(
        desc   = "Network of related enriched GO terms (needs ≥ 2 terms).",
        config = "Enrichment.top_terms, GFF.go_field"
    ),
    regionplot_img = list(
        desc   = "Zoomed regional Manhattan with gene track.",
        config = "Filter bar → Clumping distance (bp); fixed 100 kb default, tune interactively"
    ),
    hap_clustree_mg = list(
        desc   = "Marker-group cluster stability across scanned epsilon values.",
        config = "haplotype.scan.epsilon_range, .min_group_size, .min_snps"
    ),
    hap_clustree_hap = list(
        desc   = "Haplotype-group cluster stability across epsilon.",
        config = "haplotype.scan.epsilon_range, .min_group_size, .min_snps, .min_haplotype_size"
    ),
    hap_crosshap_viz = list(
        desc   = "Marker-group × haplotype × phenotype panel at the selected epsilon.",
        config = "haplotype.epsilon_selected"
    ),
    hap_boxplot = list(
        desc   = "Phenotype distribution per haplotype group.",
        config = "haplotype.epsilon_selected"
    ),
    hap_piemap = list(
        desc   = "Geographic haplotype-group frequency map.",
        config = "haplotype.epsilon_selected"
    ),
    overall_importance = list(
        desc   = "Gradient Forest variable importance ranking.",
        config = "Maladaptation.methods.gradient_forest.ntree, .cor_threshold"
    ),
    cumulative_importance = list(
        desc   = "Gradient Forest cumulative importance / turnover along each predictor.",
        config = "Maladaptation.methods.gradient_forest.ntree, .cor_threshold"
    ),
    offset_piemap = list(
        desc   = "Predicted genetic offset under future climate.",
        config = "Future.ssp, .year, .models, Maladaptation.methods.gradient_forest.spatial_correction"
    ),
    structure_piemap = list(
        desc   = "Geographic map of per-population ancestry/climate pie charts.",
        config = "Piemap.alpha, .show_labels, .pie_scale, Map.zoom_extent"
    ),
    manhattan_gea_combined = list(
        desc   = "Genome-wide GEA associations from all methods combined, significant SNPs highlighted over the full SNP cloud.",
        config = "GEA.configs"
    ),
    manhattan_gea_method = list(
        desc   = "Genome-wide GEA associations for the selected method and trait, over the full SNP cloud.",
        config = "GEA.configs"
    ),
    manhattan_gwas_combined = list(
        desc   = "Genome-wide GWAS associations from all methods combined, significant SNPs over the full SNP cloud.",
        config = "GWAS.configs"
    ),
    manhattan_gwas_method = list(
        desc   = "Genome-wide GWAS associations for the selected method and trait, over the full SNP cloud.",
        config = "GWAS.configs"
    ),
    rda_screeplot = list(
        desc   = "Constrained RDA eigenvalues; retained axes (green) vs dropped (red), per-axis anova.cca(by=\"axis\") p annotated above each bar.",
        config = "GEA.configs (RDA method params: axes, axis_alpha, condition_pcs)"
    ),
    rda_pval_hist = list(
        desc   = "rdadapt p-value histogram — GIF lambda computed with genomic.control implicit in the robust-Mahalanobis calibration (not LFMM-style genomic control). Should be flat with a spike near 0; a U-shape or hump signals miscalibration (see docs/rda_research.md A4).",
        config = "GEA.configs (RDA method params)"
    ),
    rda_biplot = list(
        desc   = "SNP loadings on the first two constrained axes (binned density, not raw points — WGS-scale marker counts), with predictor biplot vectors and candidate SNPs colored by their most-correlated predictor.",
        config = "GEA.configs (RDA method params: predictor_set)"
    ),
    miami_plot = list(
        desc   = "GEA associations (upward) mirrored against GWAS associations (downward) to reveal loci shared by climate and phenotype.",
        config = "GEAxGWAS.pairwise.window_size, .min_snps"
    ),
    pairwise_miami = list(
        desc   = "Miami plot for one selected trait pair — the two traits' associations mirrored to compare locus overlap.",
        config = "GEAxGWAS.pairwise.window_size, .min_snps"
    ),
    phenomap = list(
        desc   = "Geographic map of phenotype values across sampling sites as scaled pie charts.",
        config = "Piemap.alpha, .show_labels, .pie_scale, Map.zoom_extent"
    ),
    compare_novelty = list(
        desc   = "ExDet climate-novelty surface (NT1 univariate, NT2 multivariate) — where future climate leaves the training range and offset extrapolation is least reliable.",
        config = "Future.ssp, .year, .models"
    ),
    compare_disagree = list(
        desc   = "Rank disagreement between the two offset models mapped against climate novelty.",
        config = "Future.ssp, .year, .models"
    ),
    compare_concordance = list(
        desc   = "Per-site scatter of offset ranks from the two models with a linear fit — overall ranking agreement.",
        config = NULL
    ),
    compare_stability = list(
        desc   = "The 20 sites with the largest offset-rank difference between the two models.",
        config = NULL
    ),
    compare_nway = list(
        desc   = "Kendall's W concordance across all offset models — a single 0-1 agreement statistic over site rankings.",
        config = NULL
    ),

    # ── PreGEA (hyperparameter-exploration mode) ───────────────────────────
    pregea_screeplot = list(
        desc   = "LD-pruned genotype PCA screeplot (eigenvalue bars + broken-stick null) with the swept LFMM K band marked. Anchors the K range — the structure-axis count is an UPPER BOUND on K, not a target.",
        config = "sNMF.k_best, PreGEA.LFMM.k_offset_low, .k_offset_high"
    ),
    pregea_lfmm_hist = list(
        desc   = "P-value histogram per swept K — flat is well-calibrated, U-shaped/peaked at 0 signals residual structure or overcorrection.",
        config = "PreGEA.LFMM.k_offset_low, .k_offset_high, .genomic_control"
    ),
    pregea_lfmm_qq = list(
        desc   = "Observed vs expected -log10(p) quantiles per swept K — deviation from the diagonal at the same K as the histogram check.",
        config = "PreGEA.LFMM.k_offset_low, .k_offset_high"
    ),
    pregea_lfmm_lambda = list(
        desc   = "Genomic inflation factor (lambda_GC) vs K, fit with genomic_control OFF so lambda is actually informative. Recommended K picks the flattest p-value histogram, NOT the lambda closest to 1 — see the note on the histogram panel.",
        config = "PreGEA.LFMM.genomic_control"
    ),
    pregea_lfmm_hits = list(
        desc   = "Significant-hit count vs K at two thresholds (qvalue FDR and Bonferroni) — a secondary cross-check against the flatness-based K recommendation.",
        config = "PreGEA.LFMM.fdr, .bonf_alpha"
    ),
    pregea_emmax_hist = list(
        desc   = "P-value histogram per swept #PCs (kinship matrix held fixed) — flat is well-calibrated.",
        config = "PreGEA.EMMAX.n_pcs_min, .n_pcs_max, .kinship"
    ),
    pregea_emmax_qq = list(
        desc   = "Observed vs expected -log10(p) quantiles per swept #PCs.",
        config = "PreGEA.EMMAX.n_pcs_min, .n_pcs_max"
    ),
    pregea_emmax_lambda = list(
        desc   = "Genomic inflation factor (lambda_GC) vs #PCs — the recommended rung is the smallest #PCs landing inside the lambda tolerance band, i.e. the cheapest covariate count that still calibrates.",
        config = "PreGEA.EMMAX.lambda_tol, .deflation_floor"
    ),
    pregea_emmax_hits = list(
        desc   = "Significant-hit count vs #PCs at two thresholds — cross-check against the lambda-based recommendation.",
        config = "PreGEA.EMMAX.fdr, .bonf_alpha"
    ),
    pregea_rda_collinearity = list(
        desc   = "Predictor-pair |r| pre-screen and post-fit VIF — which climate predictors were dropped before the RDA fit, and why.",
        config = "PreGEA.RDA.collinearity_r, .vif_max, .min_predictors"
    ),
    pregea_rda_screeplot = list(
        desc   = "Constrained-axis eigenvalue screeplot at the recommended Condition()-PC rung.",
        config = "PreGEA.RDA.axis_alpha"
    ),
    pregea_rda_condition_ladder = list(
        desc   = "R2adj and rdadapt() candidate-SNP count vs number of population-structure PCs partialled out via Condition() — more PCs partialled out trades sensitivity for lower false-positive risk from shared ancestry.",
        config = "PreGEA.RDA.condition_pcs_min, .condition_pcs_max"
    ),
    pregea_rda_ordir2step = list(
        desc   = "ordiR2step forward-selection path (Blanchet double stopping rule) — which predictors entered the model, in order, and the cumulative R2adj at each step vs the full-model ceiling.",
        config = "PreGEA.Varpart.ordir2step_pin, .ordir2step_permutations, .r2_permutations"
    ),
    pregea_rda_biplot = list(
        desc   = "RDA ordination biplot at the recommended Condition()-PC rung — sample scores, predictor arrows, and candidate SNP loadings.",
        config = "PreGEA.RDA.condition_pcs_min, .condition_pcs_max"
    ),
    pregea_dbmem_screeplot = list(
        desc   = "dbMEM spatial eigenvector screeplot (Moran's I per MEM) — positive-autocorrelation MEMs are kept as the spatial predictor set for varpart.",
        config = "PreGEA.Varpart.spatial_level, .min_sites"
    ),
    pregea_dbmem_selection_path = list(
        desc   = "Forward-selected MEM path (Blanchet double stopping rule) — which spatial eigenvectors entered the varpart geography term, in order.",
        config = "PreGEA.Varpart.ordir2step_pin, .ordir2step_permutations, .r2_permutations"
    ),
    pregea_varpart_venn = list(
        desc   = "Climate / population-structure / geography variance-partition Venn (adjusted R2 per fraction) against the genomic-PC response matrix.",
        config = "PreGEA.Varpart.response, .structure_table, .confounding_flag"
    ),
    pregea_px_barplot = list(
        desc   = "Lasky et al. 2012 Px metric per climate predictor — each predictor's share of the RDA-explained variance, ranked, for the unconditional and geography-partialled models.",
        config = "PreGEA.Varpart.response"
    ),
    pregea_transfer_guard = list(
        desc   = "LFMM/EMMAX lambda re-estimated on the FULL marker set at the recommended hyperparameters, compared against the LD-pruned ladder value. Only the TREND is expected to transfer — lambda is not scale-free (Yang et al. 2011), so a large jump does not necessarily mean the recommended hyperparameter is wrong.",
        config = "PreGEA.TransferGuard.enabled, .lfmm_k, .emmax_n_pcs"
    )
)

#' Build a neutral hover-note badge for a plot: what it shows + results + config link.
#'
#' Icon-only trigger (no visible label) — this is always-present reference info,
#' not a status alert, so it stays visually quiet next to the card title. Reveals
#' a structured tooltip body on hover via filter_note() (R/utils_ui.R):
#'   - "What it shows": entry$desc (always present)
#'   - "Results": caller-supplied, data-derived from tables already on disk
#'     (e.g. best K, N samples, pops-per-K) — omitted when NULL, never invented
#'   - "Config": entry$config, the dot-path(s) that control the plot — omitted
#'     when NULL
#'
#' Returns NULL silently for an unknown id — a missing help-note should never
#' break a card's rendering.
#'
#' @param id key into HELP_NOTES
#' @param results optional data-derived results line (string or htmltools tag)
#' @param extra optional extra line appended to the tooltip body
#' @param label optional visible badge label (e.g. "K = 3") — same convention as
#'   the MAF/missingness filter badges, which show the live value on the badge
#'   itself, not just in the tooltip. Pass a value the caller already has —
#'   don't recompute it here.
#' @noRd
help_note <- function(id, results = NULL, extra = NULL, label = "") {
    entry <- HELP_NOTES[[id]]
    if (is.null(entry)) return(NULL)

    body <- htmltools::tagList(
        htmltools::p(htmltools::strong("What it shows: "), entry$desc),
        if (!is.null(results))
            htmltools::p(htmltools::strong("Results: "), results),
        if (!is.null(entry$config))
            htmltools::p(class = "text-muted small",
                htmltools::strong("Config: "), htmltools::code(entry$config)),
        if (!is.null(extra)) htmltools::p(extra)
    )

    filter_note(label = label, body = body, class = "bg-secondary")
}
