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

#' Load selected SNPs table (all sig SNPs) and pivot wide→long
#'
#' Pipeline outputs wide format: SNPID, chr, pos, {METHOD1, METHOD2, ...}, min_pvalue
#' where method columns contain the trait name (or "" if not significant for that method).
#' Returns long format: SNPID, chr, pos, pvalue, method, trait, region_id
#' @noRd
load_selected_snps <- function(project, module = MOD_ASSOC) {
    p <- selected_snps_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE,
                                colClasses = c(chr = "character", SNPID = "character"))

        # Identify method columns (everything except the fixed columns)
        fixed_cols  <- c("SNPID", "chr", "pos", "min_pvalue")
        method_cols <- setdiff(names(dt), fixed_cols)

        if (length(method_cols) == 0) return(data.table::data.table())

        # Pivot wide→long: method column name → method, value → trait
        long <- data.table::melt(dt,
            id.vars       = c("SNPID", "chr", "pos", "min_pvalue"),
            measure.vars  = method_cols,
            variable.name = "method",
            value.name    = "trait",
            variable.factor = FALSE
        )

        # Keep only rows where this method flagged the SNP (trait != "")
        long <- long[trait != ""]

        # Some methods (e.g. LFMM) store comma-separated traits when a SNP is
        # significant for multiple traits. Split those into one row per trait.
        if (any(grepl(",", long$trait, fixed = TRUE))) {
            long <- long[, {
                trs <- unlist(strsplit(trait, ",", fixed = TRUE))
                .(SNPID = SNPID, chr = chr, pos = pos, min_pvalue = min_pvalue,
                  method = method, trait = trimws(trs))
            }, by = seq_len(nrow(long))][, seq_len := NULL]
        }

        # Use min_pvalue as pvalue (background PNG handles per-method rendering)
        long[, pvalue := min_pvalue]
        long[, min_pvalue := NULL]

        # Assign region_id using coordinate-based lookup against regions_combined.
        # This ensures region_ids match the dropdown (combined format without trait suffix).
        reg_combined_path <- regions_combined_path(project, module)
        if (file.exists(reg_combined_path)) {
            regs <- data.table::fread(reg_combined_path, sep = "\t", header = TRUE,
                                      colClasses = c(chr = "character", region_id = "character"))
            if (nrow(regs) > 0 && all(c("chr", "start", "end", "region_id") %in% names(regs))) {
                regs_ov <- data.table::data.table(
                    chr       = as.character(regs$chr),
                    start     = as.integer(regs$start),
                    end       = as.integer(regs$end),
                    region_id = as.character(regs$region_id)
                )
                # One row per unique SNP position for foverlaps
                snp_pos <- unique(long[, .(SNPID, chr = as.character(chr), pos = as.integer(pos))])
                snp_ov  <- data.table::data.table(
                    SNPID = snp_pos$SNPID,
                    chr   = snp_pos$chr,
                    start = snp_pos$pos,
                    end   = snp_pos$pos
                )
                data.table::setkey(regs_ov, chr, start, end)
                data.table::setkey(snp_ov,  chr, start, end)
                ov <- data.table::foverlaps(snp_ov, regs_ov,
                    by.x = c("chr", "start", "end"),
                    by.y = c("chr", "start", "end"),
                    type = "within", mult = "first", nomatch = NA)
                rid_map <- ov[, .(SNPID, region_id)]
                long <- rid_map[long, on = "SNPID"]
            }
        }

        if (!"region_id" %in% names(long)) long[, region_id := NA_character_]
        long[, region_id := as.character(region_id)]

        long
    }, error = function(e) data.table::data.table())
}

#' Load genes per region (collapsed)
#' @noRd
load_genes <- function(project, module = MOD_ASSOC) {
    p <- genes_per_region_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE,
                                colClasses = c(chr = "character", region_id = "character"))
        # genes_per_region_collapsed uses per-trait region_ids ({chr}_{start}-{end}_{trait}).
        # Normalize to combined format ({chr}_{start}-{end}) to match region dropdown.
        if ("region_id" %in% names(dt) && nrow(dt) > 0) {
            dt[, region_id := sub("_[A-Za-z].*$", "", region_id)]
            dt <- unique(dt, by = intersect(c("gene_id", "region_id"), names(dt)))
        }
        dt
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
    p <- hap_scan_status_path(project, tag)
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
                          colClasses = c(chr = "character",
                                         gea_region_id = "character",
                                         gwas_region_id = "character"))
    }, error = function(e) data.table::data.table())
}

#' Load all per-method sig SNPs from both association and phenotype_association modules
#' Load pairwise trait overlap table (pipeline-computed)
#' @noRd
load_pairwise_table <- function(project) {
    p <- pairwise_table_path(project)
    if (!file_ok(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
            colClasses = c(trait_a = "character", source_a = "character",
                           trait_b = "character", source_b = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load pairwise collapsed sig SNPs (long format: one row per SNP per trait)
#' @noRd
load_pairwise_collapsed <- function(project) {
    p <- pairwise_collapsed_path(project)
    if (!file_ok(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
            colClasses = c(chr = "character", SNPID = "character",
                           trait = "character", source = "character")),
        error = function(e) data.table::data.table()
    )
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
