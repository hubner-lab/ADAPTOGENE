#' Discover all available projects from the pipeline path
#' @noRd
find_projects <- function(pipeline_path = get_pipeline_path()) {
    dirs <- list.dirs(pipeline_path, full.names = FALSE, recursive = FALSE)
    projects <- gsub("_results$", "", dirs[grepl("_results$", dirs)])
    sort(projects)
}

#' Find K values for a project (from structure/plots/K{k}/ directories)
#' @noRd
find_k_values <- function(project) {
    struct_plots <- mod_path(project, MOD_STRUCT, "plots")
    if (!dir.exists(struct_plots)) return(integer(0))
    k_dirs <- list.dirs(struct_plots, full.names = FALSE, recursive = FALSE)
    k_nums <- as.integer(gsub("^K", "", k_dirs[grepl("^K\\d+$", k_dirs)]))
    sort(k_nums)
}

#' Find K_start and K_end from cross_entropy plot filenames
#' @noRd
find_k_range <- function(project) {
    struct_plots <- mod_path(project, MOD_STRUCT, "plots")
    files <- list.files(struct_plots, pattern = "^cross_entropy_K.*\\.png$",
                        full.names = TRUE)
    if (length(files) == 0) return(list(k_start = NA, k_end = NA, path = NULL))
    # Match floats or integers: K2-10 or K2.0-10.0
    m <- regmatches(basename(files[1]),
                    regexpr("K([0-9.]+)-([0-9.]+)", basename(files[1])))
    if (length(m) == 0) return(list(k_start = NA, k_end = NA, path = files[1]))
    parts <- strsplit(gsub("^K", "", m), "-")[[1]]
    list(k_start = as.integer(as.numeric(parts[1])),
         k_end   = as.integer(as.numeric(parts[2])),
         path    = files[1])
}

#' Find all association methods available for a module
#' @noRd
find_assoc_methods <- function(project, module = MOD_ASSOC) {
    methods_dir <- mod_path(project, module, "tables", "methods")
    if (!dir.exists(methods_dir)) return(character(0))
    list.dirs(methods_dir, full.names = FALSE, recursive = FALSE)
}

#' Find traits from regions_per_trait.tsv (unique trait column values)
#' @noRd
find_assoc_traits <- function(project, module = MOD_ASSOC) {
    rpt <- regions_per_trait_path(project, module)
    if (!file.exists(rpt)) return(character(0))
    tryCatch({
        dt <- data.table::fread(rpt, select = "trait")
        sort(unique(dt$trait))
    }, error = function(e) character(0))
}

#' Find bio variable numbers from piemap filenames
#' @noRd
find_bio_values <- function(project) {
    piemap_dir <- mod_path(project, MOD_SK, "plots", "piemap")
    if (!dir.exists(piemap_dir)) return(integer(0))
    files <- list.files(piemap_dir, pattern = "^piemap_bio_\\d+\\.png$")
    sort(as.integer(gsub("piemap_bio_(\\d+)\\.png", "\\1", files)))
}

#' Find GF run suffixes from maladaptation output directories
#' @noRd
find_gf_suffixes <- function(project) {
    malad_dir <- mod_path(project, MOD_MALAD, "plots")
    if (!dir.exists(malad_dir)) return(character(0))
    dirs <- list.dirs(malad_dir, full.names = FALSE, recursive = FALSE)
    dirs[dirs != "zoom"]
}

#' Find haplotype tags (validated: must have scan_done.flag + HapObject files)
#' @noRd
find_haplotype_tags <- function(project) {
    hap_inter <- mod_path(project, MOD_INTER, "haplotype")
    if (!dir.exists(hap_inter)) return(character(0))
    tags <- list.dirs(hap_inter, full.names = FALSE, recursive = FALSE)
    valid <- vapply(tags, function(tag) {
        flag <- mod_path(project, MOD_INTER, "flags",
                         paste0("haplotype_", tag, "_scan_done.flag"))
        hap_files <- list.files(
            mod_path(project, MOD_INTER, "haplotype", tag),
            pattern = "HapObject\\.qs$"
        )
        file.exists(flag) && length(hap_files) > 0
    }, logical(1))
    tags[valid]
}

#' Check which pipeline modules have completed output directories
#' @noRd
check_module_status <- function(project) {
    modules <- c(
        processing     = MOD_PROC,
        structure      = MOD_STRUCT,
        structure_k    = MOD_SK,
        association    = MOD_ASSOC,
        phenotype      = MOD_PHENO,
        overlapping    = MOD_OVERLAP,
        maladaptation  = MOD_MALAD,
        haplotype      = MOD_HAP
    )
    vapply(modules, function(m) {
        dir.exists(mod_path(project, m))
    }, logical(1))
}

#' Find zoom region tags for a GF suffix
#' @noRd
find_gf_zooms <- function(project, suffix) {
    zoom_dir <- mod_path(project, MOD_MALAD, "plots", suffix, "zoom")
    if (!dir.exists(zoom_dir)) return(character(0))
    files <- list.files(zoom_dir, pattern = "\\.png$")
    gsub("\\.png$", "", files)
}

#' Find haplotype regions for a tag
#' @noRd
find_haplotype_regions <- function(project, tag) {
    sel_file <- hap_selected_regions_path(project, tag)
    if (!file.exists(sel_file)) return(character(0))
    tryCatch({
        dt <- data.table::fread(sel_file, select = "region_id")
        as.character(dt$region_id)
    }, error = function(e) character(0))
}

#' Find haplotype regions that have clustree plots (scan completed)
#' @noRd
find_haplotype_scan_regions <- function(project, tag) {
    clustree_dir <- mod_path(project, MOD_HAPSCAN, tag, "plots", "clustree")
    if (!dir.exists(clustree_dir)) return(character(0))
    files <- list.files(clustree_dir, pattern = "^Region_(.*)_clustree_MG\\.png$")
    gsub("^Region_(.*)_clustree_MG\\.png$", "\\1", files)
}

#' Find traits available for a haplotype region (from boxplot files)
#' @noRd
find_haplotype_traits <- function(project, tag, region_id) {
    plots_dir <- mod_path(project, MOD_HAP, tag, "plots")
    if (!dir.exists(plots_dir)) return(character(0))
    pat <- paste0("^Region_", region_id, "_boxplot_(.*)\\.png$")
    files <- list.files(plots_dir, pattern = pat)
    sort(gsub(pat, "\\1", files))
}
