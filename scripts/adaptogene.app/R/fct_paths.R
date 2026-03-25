# Module path constants — match workflow/rules/common.smk MOD_* values
MOD_PROC    <- "processing"
MOD_STRUCT  <- "structure"
MOD_CLIMATE <- "climate"
MOD_SK      <- "structure_k"
MOD_ASSOC   <- "association"
MOD_PHENO   <- "phenotype_association"
MOD_OVERLAP <- "overlapping"
MOD_REGPLOT <- "regionplot"
MOD_MALAD   <- "maladaptation"
MOD_HAPSCAN <- "haplotype_scan"
MOD_HAP     <- "haplotype"
MOD_INTER   <- "_intermediate"

#' Base results directory for a project
#' @noRd
project_base <- function(project, pipeline_path = get_pipeline_path()) {
    file.path(pipeline_path, paste0(project, "_results"))
}

#' Path within a module's output directory
#' @noRd
mod_path <- function(project, module, ...) {
    file.path(project_base(project), module, ...)
}

# ─── Manhattan plots ─────────────────────────────────────────────────────────

#' Per-method Manhattan background PNG path
#' @noRd
manhattan_bg_path <- function(project, module = MOD_ASSOC, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("manhattan_", trait, "_K", k, "_", adjust, "_background.png"))
}

#' Per-method Manhattan regions background PNG path
#' @noRd
manhattan_regions_bg_path <- function(project, module = MOD_ASSOC, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("manhattan_", trait, "_K", k, "_", adjust, "_regions_background.png"))
}

#' Per-method coordinate mapping JSON path
#' @noRd
manhattan_coords_path <- function(project, module = MOD_ASSOC, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("manhattan_", trait, "_K", k, "_", adjust, "_coords.json"))
}

#' Combined Manhattan background PNG path
#' @noRd
combined_manhattan_bg_path <- function(project, module = MOD_ASSOC, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("manhattan_combined_K", k, "_background.png"))
}

#' Combined Manhattan regions background PNG path
#' @noRd
combined_manhattan_regions_bg_path <- function(project, module = MOD_ASSOC, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("manhattan_combined_K", k, "_regions_background.png"))
}

#' Combined Manhattan coordinate mapping JSON path
#' @noRd
combined_manhattan_coords_path <- function(project, module = MOD_ASSOC, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("manhattan_combined_K", k, "_coords.json"))
}

#' QQ plot path (per-method)
#' @noRd
qq_plot_path <- function(project, module = MOD_ASSOC, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("qq_", trait, "_K", k, "_", adjust, ".png"))
}

#' Combined QQ plot path
#' @noRd
qq_combined_path <- function(project, module = MOD_ASSOC, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("qq_combined_K", k, ".png"))
}

# ─── Miami plot ───────────────────────────────────────────────────────────────

#' Miami background PNG path
#' @noRd
miami_bg_path <- function(project, k) {
    mod_path(project, MOD_OVERLAP, "plots",
             paste0("miami_combined_K", k, "_background.png"))
}

#' Miami regions background PNG path
#' @noRd
miami_regions_bg_path <- function(project, k) {
    mod_path(project, MOD_OVERLAP, "plots",
             paste0("miami_combined_K", k, "_regions_background.png"))
}

#' Miami coordinate mapping JSON path
#' @noRd
miami_coords_path <- function(project, k) {
    mod_path(project, MOD_OVERLAP, "plots",
             paste0("miami_combined_K", k, "_coords.json"))
}

# ─── Association tables ───────────────────────────────────────────────────────

#' Glob all per-method sig SNPs files for a module (across all methods and adjustments)
#' @noRd
find_method_sigsnps_files <- function(project, module = MOD_ASSOC) {
    Sys.glob(mod_path(project, module, "tables", "methods", "*", "*_sig_snps_*.tsv"))
}

#' Pairwise trait overlap table (pipeline-computed)
#' @noRd
pairwise_table_path <- function(project) {
    mod_path(project, MOD_OVERLAP, "tables", "pairwise_overlap_table.tsv")
}

#' Pairwise collapsed sig SNPs (long format, one row per SNP per trait)
#' @noRd
pairwise_collapsed_path <- function(project) {
    mod_path(project, MOD_OVERLAP, "tables", "pairwise_collapsed_snps.tsv")
}

#' Selected SNPs table (combined all methods)
#' @noRd
selected_snps_path <- function(project, module = MOD_ASSOC) {
    fname <- if (module == MOD_OVERLAP) "selected_snps_all.tsv" else "selected_snps.tsv"
    mod_path(project, module, "tables", fname)
}

#' Per-method sig SNPs table
#' @noRd
method_sigsnps_path <- function(project, module = MOD_ASSOC, method, adjust) {
    mod_path(project, module, "tables", "methods", method,
             paste0(method, "_pvalues_K*_sig_snps_", adjust, ".tsv"))
}

#' Per-trait regions table
#' @noRd
regions_per_trait_path <- function(project, module = MOD_ASSOC) {
    fname <- if (module == MOD_OVERLAP) "regions_per_trait_all.tsv" else "regions_per_trait.tsv"
    mod_path(project, module, "tables", fname)
}

#' Combined regions table
#' @noRd
regions_combined_path <- function(project, module = MOD_ASSOC) {
    mod_path(project, module, "tables", "regions_combined.tsv")
}

#' Genes per region table
#' @noRd
genes_per_region_path <- function(project, module = MOD_ASSOC) {
    mod_path(project, module, "tables", "genes_per_region_collapsed.tsv")
}

#' GO enrichment table for a specific region and trait
#' @noRd
enrichment_table_path <- function(project, module = MOD_ASSOC, trait, region_id) {
    mod_path(project, module, "tables", "enrichment", trait,
             paste0("Region_", region_id, "_", trait, "_enrichment.tsv"))
}

# ─── Enrichment plots ─────────────────────────────────────────────────────────

#' Enrichment plot path
#' @noRd
enrichment_plot_path <- function(project, module = MOD_ASSOC, trait, region_id, plot_type) {
    mod_path(project, module, "plots", "enrichment", trait,
             paste0("Region_", region_id, "_", trait, "_", plot_type, ".png"))
}

# ─── Structure ────────────────────────────────────────────────────────────────

#' PCA plot path
#' @noRd
pca_path <- function(project) {
    mod_path(project, MOD_STRUCT, "plots", "pca.png")
}

#' Tracy-Widom path
#' @noRd
tracy_widom_path <- function(project) {
    mod_path(project, MOD_STRUCT, "plots", "tracy_widom.png")
}

#' Cross-entropy plot path
#' @noRd
cross_entropy_path <- function(project, k_start, k_end) {
    mod_path(project, MOD_STRUCT, "plots",
             paste0("cross_entropy_K", k_start, "-", k_end, ".png"))
}

#' Structure barplot path for a given K
#' @noRd
structure_k_path <- function(project, k) {
    mod_path(project, MOD_STRUCT, "plots", paste0("K", k),
             paste0("structure_K", k, ".png"))
}

#' PCA structure plot path for a given K
#' @noRd
pca_structure_path <- function(project, k) {
    mod_path(project, MOD_STRUCT, "plots", paste0("K", k),
             paste0("pca_structure_K", k, ".png"))
}

#' Pop differentiation path for a given K
#' @noRd
pop_diff_path <- function(project, k) {
    mod_path(project, MOD_STRUCT, "plots", paste0("K", k),
             paste0("pop_diff_K", k, ".png"))
}

# ─── Structure K ─────────────────────────────────────────────────────────────

#' Piemap path
#' @noRd
piemap_path <- function(project, bio, metric = NULL, zoom = NULL) {
    dir <- mod_path(project, MOD_SK, "plots", "piemap")
    if (!is.null(zoom) && nzchar(zoom) && zoom != "none") {
        dir <- file.path(dir, "zoom", zoom)
        fname <- paste0("piemap_bio_", bio,
                        if (!is.null(metric) && nzchar(metric) && metric != "none") paste0("_", metric) else "", ".png")
    } else if (!is.null(metric) && nzchar(metric) && metric != "none") {
        fname <- paste0("piemap_bio_", bio, "_", metric, ".png")
    } else {
        fname <- paste0("piemap_bio_", bio, ".png")
    }
    file.path(dir, fname)
}

#' Climate correlation heatmap path
#' @noRd
climate_heatmap_path <- function(project) {
    mod_path(project, MOD_CLIMATE, "plots", "correlation_heatmap.png")
}

#' Climate density plot path
#' @noRd
climate_density_path <- function(project, bio = NULL) {
    if (is.null(bio)) {
        mod_path(project, MOD_CLIMATE, "plots", "density_plot_present.png")
    } else {
        mod_path(project, MOD_CLIMATE, "plots", paste0("density_plot_present_bio", bio, ".png"))
    }
}

#' LD decay plot path
#' @noRd
ld_decay_path <- function(project, per_chr = FALSE) {
    fname <- if (per_chr) "ld_decay_per_chromosome.png" else "ld_decay_genome_wide.png"
    mod_path(project, MOD_SK, "plots", "ld_decay", fname)
}

# ─── Maladaptation ───────────────────────────────────────────────────────────

#' GF importance plot path
#' @noRd
gf_importance_path <- function(project, suffix, type = "overall") {
    mod_path(project, MOD_MALAD, "plots", suffix,
             paste0(type, "_importance.png"))
}

#' GF genetic offset piemap path
#' @noRd
gf_offset_piemap_path <- function(project, suffix, variant = "base") {
    # "base" maps to no suffix (genetic_offset_piemap.png); other variants appended
    fname <- if (is.null(variant) || variant == "base" || !nzchar(variant)) {
        "genetic_offset_piemap.png"
    } else {
        paste0("genetic_offset_piemap_", variant, ".png")
    }
    mod_path(project, MOD_MALAD, "plots", suffix, fname)
}

#' GF site offset table path
#' @noRd
gf_site_table_path <- function(project, suffix) {
    mod_path(project, MOD_MALAD, "tables", suffix, "genetic_offset_site.tsv")
}

# ─── Phenotype piemap ────────────────────────────────────────────────────────

#' Phenotype association piemap (phenomap) path
#' @noRd
pheno_piemap_path <- function(project, trait) {
    mod_path(project, MOD_PHENO, "plots", "piemap",
             paste0("phenomap_", trait, ".png"))
}

# ─── Regional Manhattan (topr) ────────────────────────────────────────────────

#' Regionplot image path
#' @noRd
regionplot_path <- function(project, region_id, trait = NULL) {
    # Sanitize region_id: pipeline replaces colons and hyphens with underscores
    safe_id <- gsub("[:\\-]", "_", region_id)
    fname <- if (is.null(trait)) {
        paste0("regionplot_", safe_id, ".png")
    } else {
        paste0("regionplot_", safe_id, "_", trait, ".png")
    }
    mod_path(project, MOD_REGPLOT, fname)
}

# ─── Haplotype ────────────────────────────────────────────────────────────────

#' Haplotype clustree path (per-region)
#' @noRd
hap_clustree_path <- function(project, tag, region_id, type = "MG") {
    mod_path(project, MOD_HAPSCAN, tag, "plots", "clustree",
             paste0("Region_", region_id, "_clustree_", type, ".png"))
}

#' Haplotype crosshap_viz path (per-trait)
#' @noRd
crosshap_viz_path <- function(project, tag, region_id, trait) {
    mod_path(project, MOD_HAP, tag, "plots",
             paste0("Region_", region_id, "_crosshap_viz_", trait, ".png"))
}

#' Haplotype boxplot path (per-trait)
#' @noRd
hap_boxplot_path <- function(project, tag, region_id, trait) {
    mod_path(project, MOD_HAP, tag, "plots",
             paste0("Region_", region_id, "_boxplot_", trait, ".png"))
}

#' Haplotype piemap path (per-trait)
#' @noRd
hap_piemap_path <- function(project, tag, region_id, trait) {
    mod_path(project, MOD_HAP, tag, "plots",
             paste0("haplotype_piemap_region_", region_id, "_", trait, ".png"))
}

#' Haplotype scan status path
#' @noRd
hap_scan_status_path <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "scan_status.tsv")
}

#' Haplotype selected regions path
#' @noRd
hap_selected_regions_path <- function(project, tag) {
    mod_path(project, MOD_HAPSCAN, tag, "tables", "selected_regions.tsv")
}

#' Haplotype assignments table path
#' @noRd
hap_assignments_path <- function(project, tag) {
    mod_path(project, MOD_HAP, tag, "tables", "Haplotype_assignments.tsv")
}

#' Haplotype frequencies table path
#' @noRd
hap_frequencies_path <- function(project, tag) {
    mod_path(project, MOD_HAP, tag, "tables", "Haplotype_frequencies.tsv")
}

#' Pipeline summary TSV path
#' @noRd
pipeline_summary_path <- function(project) {
    file.path(project_base(project), "pipeline_summary.tsv")
}
