# Module path constants — match workflow/rules/common.smk MOD_* values
MOD_PROC      <- "Processing"
MOD_PRESTRUCT <- "PreStructure"
MOD_STRUCT    <- "Structure"
MOD_CLIMATE   <- "climate"
MOD_GEA       <- "GEA"
MOD_GWAS      <- "GWAS"
MOD_GEAXGWAS  <- "GEAxGWAS"
MOD_MALAD     <- "Maladaptation"
MOD_INTER     <- "_intermediate"

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
manhattan_bg_path <- function(project, module = MOD_GEA, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("manhattan_", trait, "_K", k, "_", adjust, "_background.png"))
}

#' Per-method coordinate mapping JSON path
#' @noRd
manhattan_coords_path <- function(project, module = MOD_GEA, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("manhattan_", trait, "_K", k, "_", adjust, "_coords.json"))
}

#' Combined Manhattan background PNG path
#' @noRd
combined_manhattan_bg_path <- function(project, module = MOD_GEA, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("manhattan_combined_K", k, "_background.png"))
}

#' Combined Manhattan coordinate mapping JSON path
#' @noRd
combined_manhattan_coords_path <- function(project, module = MOD_GEA, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("manhattan_combined_K", k, "_coords.json"))
}

#' QQ plot path (per-method)
#' @noRd
qq_plot_path <- function(project, module = MOD_GEA, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("qq_", trait, "_K", k, "_", adjust, ".png"))
}

#' Combined QQ plot path
#' @noRd
qq_combined_path <- function(project, module = MOD_GEA, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("qq_combined_K", k, ".png"))
}

# ─── Miami plot ───────────────────────────────────────────────────────────────

#' Miami background PNG path
#' @noRd
miami_bg_path <- function(project, k) {
    mod_path(project, MOD_GEAXGWAS, "plots",
             paste0("miami_combined_K", k, "_background.png"))
}

#' Miami coordinate mapping JSON path
#' @noRd
miami_coords_path <- function(project, k) {
    mod_path(project, MOD_GEAXGWAS, "plots",
             paste0("miami_combined_K", k, "_coords.json"))
}

# ─── Association tables ───────────────────────────────────────────────────────

#' Glob all per-method sig SNPs files for a module (across all methods and adjustments)
#' @noRd
find_method_sigsnps_files <- function(project, module = MOD_GEA) {
    Sys.glob(mod_path(project, module, "tables", "methods", "*", "*_sig_snps_*.tsv"))
}

#' Pairwise trait overlap table (pipeline-computed)
#' @noRd
pairwise_table_path <- function(project) {
    mod_path(project, MOD_GEAXGWAS, "tables", "pairwise_overlap_table.tsv")
}

#' Pairwise collapsed sig SNPs (long format, one row per SNP per trait)
#' @noRd
pairwise_collapsed_path <- function(project) {
    mod_path(project, MOD_GEAXGWAS, "tables", "pairwise_collapsed_snps.tsv")
}

#' Selected SNPs table (combined all methods)
#' @noRd
selected_snps_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "selected_snps.tsv")
}

#' Per-method sig SNPs table (glob pattern for discovery)
#' @noRd
method_sigsnps_path <- function(project, module = MOD_GEA, method, adjust) {
    mod_path(project, module, "tables", "methods", method,
             paste0(method, "_pvalues_K*_sig_snps_", adjust, ".tsv"))
}

#' Per-method sig SNPs table (exact K — for loading in Shiny overlay)
#' @noRd
method_sigsnps_direct_path <- function(project, module = MOD_GEA, method, k, adjust) {
    mod_path(project, module, "tables", "methods", method,
             paste0(method, "_pvalues_K", k, "_sig_snps_", adjust, ".tsv"))
}

#' Per-trait regions table
#' @noRd
regions_per_trait_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "regions_per_trait.tsv")
}

#' Combined regions table
#' @noRd
regions_combined_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "regions_combined.tsv")
}

#' Genes per region table
#' @noRd
genes_per_region_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "genes_per_region_collapsed.tsv")
}

#' GO enrichment table for a specific region and trait
#' @noRd
enrichment_table_path <- function(project, module = MOD_GEA, trait, region_id) {
    mod_path(project, module, "tables", "enrichment", trait,
             paste0("Region_", region_id, "_", trait, "_enrichment.tsv"))
}

# ─── Enrichment plots ─────────────────────────────────────────────────────────

#' Enrichment plot path
#' @noRd
enrichment_plot_path <- function(project, module = MOD_GEA, trait, region_id, plot_type) {
    mod_path(project, module, "plots", "enrichment", trait,
             paste0("Region_", region_id, "_", trait, "_", plot_type, ".png"))
}

# ─── PreStructure ─────────────────────────────────────────────────────────────

#' PCA plot path
#' @noRd
pca_path <- function(project) {
    mod_path(project, MOD_PRESTRUCT, "plots", "pca.png")
}

#' Tracy-Widom path
#' @noRd
tracy_widom_path <- function(project) {
    mod_path(project, MOD_PRESTRUCT, "plots", "tracy_widom.png")
}

#' Cross-entropy plot path
#' @noRd
cross_entropy_path <- function(project, k_start, k_end) {
    mod_path(project, MOD_PRESTRUCT, "plots",
             paste0("cross_entropy_K", k_start, "-", k_end, ".png"))
}

#' Structure barplot path for a given K
#' @noRd
structure_k_path <- function(project, k) {
    mod_path(project, MOD_PRESTRUCT, "plots", paste0("K", k),
             paste0("structure_K", k, ".png"))
}

#' PCA structure plot path for a given K
#' @noRd
pca_structure_path <- function(project, k) {
    mod_path(project, MOD_PRESTRUCT, "plots", paste0("K", k),
             paste0("pca_structure_K", k, ".png"))
}

#' Pop differentiation path for a given K
#' @noRd
pop_diff_path <- function(project, k) {
    mod_path(project, MOD_PRESTRUCT, "plots", paste0("K", k),
             paste0("pop_diff_K", k, ".png"))
}

# ─── Structure ────────────────────────────────────────────────────────────────

#' Piemap path
#' @noRd
piemap_path <- function(project, bio, metric = NULL, zoom = NULL) {
    dir <- mod_path(project, MOD_STRUCT, "plots", "piemap")
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

#' Phenotype density plot path
#' @noRd
phenotype_density_path <- function(project) {
    mod_path(project, MOD_CLIMATE, "plots", "density_plot_phenotypes.png")
}

#' Climate invariant predictors file path
#' @noRd
climate_invariant_path <- function(project) {
    mod_path(project, MOD_CLIMATE, "tables", "present", "climate_invariant_predictors.tsv")
}

#' Phenotype missing summary path (GWAS mode)
#' @noRd
pheno_missing_summary_path <- function(project) {
    mod_path(project, MOD_GWAS, "tables", "phenotype_missing_summary.tsv")
}

#' LEA .removed file path — glob from _work/ (LD-pruned variant)
#' Returns empty string if not found.
#' @noRd
removed_snps_path <- function(project) {
    base  <- project_base(project)
    files <- Sys.glob(file.path(base, "_work", "maf*", paste0(project, ".removed")))
    if (length(files) == 0) "" else files[1]
}

#' Glob first per-SNP pvalue TSV for a module/method (any K, any adjust).
#' Used to extract trait column names. Excludes snp-mode `_sig_snps_` and
#' @noRd
find_pvalue_tsv <- function(project, module = MOD_GEA) {
    files <- Sys.glob(mod_path(project, module, "tables", "methods", "*", "*_pvalues_K*.tsv"))
    files <- grep("_sig_snps_", files, value = TRUE, invert = TRUE)
    if (length(files) == 0) "" else files[1]
}

#' LD decay plot path
#' @noRd
ld_decay_path <- function(project, per_chr = FALSE) {
    fname <- if (per_chr) "ld_decay_per_chromosome.png" else "ld_decay_genome_wide.png"
    mod_path(project, MOD_STRUCT, "plots", "ld_decay", fname)
}

#' LD decay half-distances TSV path
#' @noRd
ld_decay_table_path <- function(project) {
    mod_path(project, MOD_STRUCT, "tables", "ld_decay_half_distances.tsv")
}

# ─── Maladaptation ───────────────────────────────────────────────────────────

#' GF selected SNPs path (intermediate file under _intermediate/{method}/{suffix}/)
#' @noRd
gf_selected_snps_path <- function(project, suffix, method = "gradient_forest") {
    mod_path(project, MOD_INTER, method, suffix, "selected_snps.tsv")
}

#' GF importance plot path
#' @noRd
gf_importance_path <- function(project, suffix, type = "overall", method = "gradient_forest") {
    mod_path(project, MOD_MALAD, "plots", method, suffix,
             paste0(type, "_importance.png"))
}

#' GF genetic offset piemap path
#' @noRd
gf_offset_piemap_path <- function(project, suffix, variant = "base", method = "gradient_forest") {
    # "base" maps to "notrait" suffix; other variants used as-is
    piemap_variant <- if (is.null(variant) || variant == "base" || !nzchar(variant)) "notrait" else variant
    fname <- paste0("genetic_offset_piemap_", piemap_variant, ".png")
    mod_path(project, MOD_MALAD, "plots", method, suffix, fname)
}

#' GF site offset table path
#' @noRd
gf_site_table_path <- function(project, suffix, method = "gradient_forest") {
    mod_path(project, MOD_MALAD, "tables", method, suffix, "genetic_offset_site.tsv")
}

# ─── Phenotype piemap ────────────────────────────────────────────────────────

#' Phenotype association piemap (phenomap) path
#' @noRd
pheno_piemap_path <- function(project, trait) {
    mod_path(project, MOD_GWAS, "plots", "piemap",
             paste0("phenomap_", trait, ".png"))
}

# ─── Input data paths ─────────────────────────────────────────────────────────

#' Resolve full path to the project's GFF annotation file
#' @noRd
gff_path <- function(project, pipeline_path = get_pipeline_path()) {
    config  <- read_project_config(project, pipeline_path)
    inp_dir <- config_get(config, "Input", "dir", default = "data")
    gff_rel <- config_get(config, "Input", "gff", default = "")
    if (is.null(gff_rel) || gff_rel == "") return("")
    file.path(pipeline_path, inp_dir, gff_rel)
}

# ─── Regional Manhattan (topr) ────────────────────────────────────────────────

#' Directory for on-demand regionplots persisted to disk.
#' Subdirectory per combo_hash keeps files for different filter combinations separate.
#' @noRd
ondemand_regionplot_dir <- function(project, module = MOD_GEA, combo_hash) {
    mod_path(project, "regionplot", module, combo_hash)
}

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
    mod_path(project, "regionplot", fname)
}

# ─── Haplotype ────────────────────────────────────────────────────────────────
# All haplotype outputs live under _intermediate/haplotype/{tag}/ (Shiny on-demand only;
# no pipeline haplotype modes since Phase 3).

#' Haplotype clustree path (per-region)
#' @noRd
hap_clustree_path <- function(project, tag, region_id, type = "MG") {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots", "clustree",
             paste0("Region_", region_id, "_clustree_", type, ".png"))
}

#' Haplotype crosshap_viz path (per-trait)
#' @noRd
crosshap_viz_path <- function(project, tag, region_id, trait) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots",
             paste0("Region_", region_id, "_crosshap_viz_", trait, ".png"))
}

#' Haplotype boxplot path (per-trait)
#' @noRd
hap_boxplot_path <- function(project, tag, region_id, trait) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots",
             paste0("Region_", region_id, "_boxplot_", trait, ".png"))
}

#' Haplotype piemap path (per-trait)
#' @noRd
hap_piemap_path <- function(project, tag, region_id, trait) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots",
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
    mod_path(project, MOD_INTER, "haplotype", tag, "tables", "selected_regions.tsv")
}

#' Haplotype assignments table path
#' @noRd
hap_assignments_path <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "tables", "Haplotype_assignments.tsv")
}

#' Haplotype frequencies table path
#' @noRd
hap_frequencies_path <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "tables", "Haplotype_frequencies.tsv")
}

#' Haplotype intermediate directory (HapObject .qs files)
#' @noRd
hap_intermediate_dir <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag)
}

#' Haplotype scan clustree plots output directory
#' @noRd
hap_scan_plots_dir <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots", "clustree")
}

#' Haplotype viz plots output directory
#' @noRd
hap_viz_plots_dir <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots")
}

#' Haplotype viz tables output directory
#' @noRd
hap_viz_tables_dir <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "tables")
}

#' Derive VCF basename from config (mirrors pipeline's VCF_BASE).
#' Falls back to `project` if config is missing an Input.vcf entry.
#' @noRd
.vcf_base <- function(project, config = NULL) {
    vcf <- if (!is.null(config)) config_get(config, "Input", "vcf", default = NULL) else NULL
    if (is.null(vcf) || !nzchar(vcf)) return(project)
    name <- basename(vcf)
    name <- sub("\\.vcf\\.gz$", "", name)
    name <- sub("\\.vcf$",     "", name)
    name
}

#' Filtered VCF path (glob from _work, returns first match)
#' @noRd
hap_filtered_vcf_path <- function(project, config = NULL) {
    base <- project_base(project)
    vb <- .vcf_base(project, config)
    vcfs <- Sys.glob(file.path(base, "_work", "maf*", paste0(vb, ".vcf")))
    if (length(vcfs) == 0) "" else vcfs[1]
}

#' Imputed VCF path for LD computation (glob from _work, returns first match)
#' @noRd
hap_imputed_vcf_path <- function(project, k_best, config = NULL) {
    base <- project_base(project)
    vb <- .vcf_base(project, config)
    vcfs <- Sys.glob(file.path(base, "_work", "maf*",
                               paste0(vb, "_K", k_best, "imp.vcf")))
    if (length(vcfs) == 0) "" else vcfs[1]
}

#' Aligned metadata path
#' @noRd
hap_metadata_path <- function(project) {
    mod_path(project, MOD_PROC, "tables", "metadata.tsv")
}

#' Clusters table path for k_best
#' @noRd
hap_clusters_path <- function(project, k_best) {
    mod_path(project, MOD_PRESTRUCT, "tables", paste0("K", k_best),
             paste0("clusters_K", k_best, ".tsv"))
}

# ── Processing QC paths ──────────────────────────────────────────────────────

#' QC plot path (Processing/plots/<filename>)
#' @noRd
qc_plot_path <- function(project, filename) {
    mod_path(project, MOD_PROC, "plots", filename)
}

#' QC table path (Processing/tables/<filename>)
#' @noRd
qc_table_path <- function(project, filename) {
    mod_path(project, MOD_PROC, "tables", filename)
}

#' Pipeline summary TSV path
#' @noRd
pipeline_summary_path <- function(project) {
    file.path(project_base(project), "pipeline_summary.tsv")
}
