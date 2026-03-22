#' Load pipeline summary TSV
#' @noRd
load_pipeline_summary <- function(project) {
    p <- pipeline_summary_path(project)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE),
        error = function(e) data.table::data.table()
    )
}

#' Load regions per trait
#' @noRd
load_regions_per_trait <- function(project, module = MOD_ASSOC) {
    p <- regions_per_trait_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE,
                                colClasses = c(chr = "character", region_id = "character"))
        dt
    }, error = function(e) data.table::data.table())
}

#' Load combined regions
#' @noRd
load_regions_combined <- function(project, module = MOD_ASSOC) {
    p <- regions_combined_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(chr = "character", region_id = "character"))
    }, error = function(e) data.table::data.table())
}

#' Load selected SNPs table (all sig SNPs with region assignments)
#' @noRd
load_selected_snps <- function(project, module = MOD_ASSOC) {
    p <- selected_snps_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(chr = "character", region_id = "character",
                                         SNPID = "character"))
    }, error = function(e) data.table::data.table())
}

#' Load genes per region (collapsed)
#' @noRd
load_genes <- function(project, module = MOD_ASSOC) {
    p <- genes_per_region_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(chr = "character", region_id = "character"))
    }, error = function(e) data.table::data.table())
}

#' Load GO enrichment table for a specific region and trait
#' @noRd
load_enrichment_table <- function(project, module = MOD_ASSOC, trait, region_id) {
    p <- enrichment_table_path(project, module, trait, region_id)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE)
    }, error = function(e) data.table::data.table())
}

#' Load haplotype scan status table
#' @noRd
load_haplotype_status <- function(project, tag) {
    p <- mod_path(project, MOD_HAPSCAN, tag, "scan_status.tsv")
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(region_id = "character", chr = "character"))
    }, error = function(e) data.table::data.table())
}

#' Load GF site offset table
#' @noRd
load_gf_site_table <- function(project, suffix) {
    p <- gf_site_table_path(project, suffix)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(site = "character"))
    }, error = function(e) data.table::data.table())
}

#' Load overlap summary table
#' @noRd
load_overlap_summary <- function(project) {
    p <- mod_path(project, MOD_OVERLAP, "tables", "overlap_summary.tsv")
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(chr = "character", region_id = "character"))
    }, error = function(e) data.table::data.table())
}

#' Cached data loader using cachem
#' Creates a session-level cache automatically on first call.
#' @noRd
.cache_env <- new.env(parent = emptyenv())

get_app_cache <- function() {
    if (is.null(.cache_env$cache)) {
        .cache_env$cache <- cachem::cache_mem(max_size = 200 * 1024^2)  # 200 MB
    }
    .cache_env$cache
}

#' Load with session cache (avoids re-reading the same file across tab switches)
#' @noRd
load_cached <- function(cache_key, loader_fn) {
    cache <- get_app_cache()
    safe_key <- gsub("[^a-z0-9]", "", tolower(cache_key))
    cached <- cache$get(safe_key)
    if (!cachem::is.key_missing(cached)) return(cached)
    val <- loader_fn()
    cache$set(safe_key, val)
    val
}
