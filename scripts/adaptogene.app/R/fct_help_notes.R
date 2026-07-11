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
        desc   = "Sample scatter on genotype PCs — unsupervised structure overview.",
        config = NULL
    ),
    structure_bar = list(
        desc   = "sNMF ancestry proportions per sample at the selected K.",
        config = "sNMF.k_best"
    ),
    pca_structure = list(
        desc   = "PCA scatter colored by K-cluster assignment.",
        config = "sNMF.k_best"
    ),
    pop_diff = list(
        desc   = "Population differentiation among K clusters.",
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
        config = "GEA.snp_clumping_distance / GWAS.snp_clumping_distance"
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
    )
)

#' Build a neutral hover-note badge for a plot: what it shows + its config link.
#'
#' Icon-only trigger (no visible label) — this is always-present reference info,
#' not a status alert, so it stays visually quiet next to the card title. Reveals
#' the description + config path on hover via filter_note() (R/utils_ui.R).
#'
#' Returns NULL silently for an unknown id — a missing help-note should never
#' break a card's rendering.
#'
#' @param id key into HELP_NOTES
#' @param extra optional extra line appended to the tooltip body
#' @param label optional visible badge label (e.g. "K = 3") — same convention as
#'   the MAF/missingness filter badges, which show the live value on the badge
#'   itself, not just in the tooltip. Pass a value the caller already has —
#'   don't recompute it here.
#' @noRd
help_note <- function(id, extra = NULL, label = "") {
    entry <- HELP_NOTES[[id]]
    if (is.null(entry)) return(NULL)

    body <- htmltools::tagList(
        htmltools::p(entry$desc),
        if (!is.null(entry$config))
            htmltools::p(
                htmltools::span(class = "text-muted small",
                    "Config: ", htmltools::code(entry$config))
            ),
        if (!is.null(extra)) htmltools::p(extra)
    )

    filter_note(label = label, body = body, class = "bg-secondary")
}
