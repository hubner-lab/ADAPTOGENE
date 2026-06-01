# Source shared pipeline utilities for interactive threshold computation.
# compute_pval_threshold() is used directly; lib/sig_snps.R is not needed
# here (we implement a simpler single-threaded version inline).
if (file.exists("/pipeline/scripts/R/utils/pval_threshold.R")) {
    source("/pipeline/scripts/R/utils/pval_threshold.R")
}

# ── Processing QC loaders ────────────────────────────────────────────────────

#' Load filtering summary table
#' @noRd
load_filtering_summary <- function(project) {
    p <- qc_table_path(project, "filtering_summary.tsv")
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE),
        error = function(e) data.table::data.table()
    )
}

#' Load sample heterozygosity table
#' @noRd
load_sample_heterozygosity <- function(project) {
    p <- qc_table_path(project, "sample_heterozygosity.tsv")
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(sample = "character", site = "character")),
        error = function(e) data.table::data.table()
    )
}

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
load_regions_per_trait <- function(project, module = MOD_GEA) {
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
load_regions_combined <- function(project, module = MOD_GEA) {
    p <- regions_combined_path(project, module)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(chr = "character", region_id = "character"))
    }, error = function(e) data.table::data.table())
}

#' Assign region_id to a SNP data.table via coordinate-based foverlaps lookup.
#'
#' Uses regions_combined.tsv to map each SNP position to its enclosing region.
#' Region IDs are in combined format (chr_start-end, no trait suffix).
#' @param snps_dt data.table with columns SNPID, chr, pos
#' @param project project name
#' @param module pipeline module (MOD_GEA, MOD_GWAS, MOD_GEAXGWAS)
#' @return snps_dt with region_id column added (NA for SNPs outside any region)
#' @noRd
assign_region_ids <- function(snps_dt, project, module) {
    reg_combined_path <- regions_combined_path(project, module)
    snps_dt[, region_id := NA_character_]
    if (!file.exists(reg_combined_path)) return(snps_dt)
    regs <- tryCatch(
        data.table::fread(reg_combined_path, sep = "\t", header = TRUE,
                          colClasses = c(chr = "character", region_id = "character")),
        error = function(e) data.table::data.table()
    )
    if (nrow(regs) == 0 || !all(c("chr", "start", "end", "region_id") %in% names(regs)))
        return(snps_dt)

    regs_ov <- data.table::data.table(
        chr       = as.character(regs$chr),
        start     = as.integer(regs$start),
        end       = as.integer(regs$end),
        region_id = as.character(regs$region_id)
    )
    snp_pos <- unique(snps_dt[, .(SNPID, chr = as.character(chr), pos = as.integer(pos))])
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
    snps_dt[rid_map, region_id := i.region_id, on = "SNPID"]
    snps_dt[, region_id := as.character(region_id)]
    snps_dt
}

#' Load selected SNPs table (all sig SNPs) and pivot wide→long
#'
#' Pipeline outputs wide format: SNPID, chr, pos, {METHOD1, METHOD2, ...}, min_pvalue
#' where method columns contain the trait name (or "" if not significant for that method).
#' Returns long format: SNPID, chr, pos, pvalue, method, trait, region_id
#' @noRd
load_selected_snps <- function(project, module = MOD_GEA, k_best = NA) {
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

        # Join actual per-method p-values from sig_snps files.
        # min_pvalue is kept for combined Manhattan y-range but NOT shown in tables.
        method_dts <- load_all_method_sigsnps(project, module, k = k_best)
        if (length(method_dts) > 0) {
            method_pvals <- data.table::rbindlist(method_dts, use.names = TRUE)
            method_pvals <- method_pvals[, .(SNPID, method, trait, pvalue)]
            long <- merge(long, method_pvals,
                          by = c("SNPID", "method", "trait"), all.x = TRUE)
            # Rows without a match keep min_pvalue as fallback
            long[is.na(pvalue), pvalue := min_pvalue]
        } else {
            long[, pvalue := min_pvalue]
        }
        # Keep min_pvalue for combined Manhattan overlay y-placement
        # (combined background PNG y-axis is calibrated to min across all methods)

        long <- assign_region_ids(long, project, module)
        long
    }, error = function(e) data.table::data.table())
}

#' Load per-method sig SNPs (with actual per-method p-values)
#'
#' Reads the per-method sig_snps TSV produced by find_sig_snps.R, which contains
#' the actual p-value for that specific method (not min across all methods).
#' This is required for correct y-axis placement in per-method Manhattan overlays,
#' because each method's background PNG has a y-range calibrated to its own p-values.
#' Returns long format: SNPID, chr, pos, pvalue, method, trait, region_id
#' @noRd
load_method_sigsnps <- function(project, module = MOD_GEA, method, k, adjust) {
    p <- method_sigsnps_direct_path(project, module, method, k, adjust)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch({
        dt <- data.table::fread(p, sep = "\t", header = TRUE,
                                colClasses = c(chr = "character", SNPID = "character",
                                               method = "character", trait = "character"))
        if (nrow(dt) == 0) return(data.table::data.table())

        # Keep only the columns needed for the overlay
        keep <- intersect(c("SNPID", "chr", "pos", "pvalue", "method", "trait"), names(dt))
        dt <- dt[, ..keep]

        dt <- assign_region_ids(dt, project, module)
        dt
    }, error = function(e) data.table::data.table())
}

#' Load genes per region (collapsed)
#' @noRd
load_genes <- function(project, module = MOD_GEA) {
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
load_enrichment_table <- function(project, module = MOD_GEA, trait, region_id) {
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

#' Parse GFF9 attribute string into a data.table of key-value columns
#'
#' Each row's attributes are parsed independently, so variable key order and
#' missing keys (filled with NA) are handled correctly.
#'
#' @param attr_strings character vector of raw GFF attribute strings
#' @return data.table with one column per attribute key
#' @noRd
parse_gff_attributes <- function(attr_strings) {
    parsed <- lapply(attr_strings, function(s) {
        if (is.na(s) || !nzchar(s)) return(list())
        pairs <- strsplit(s, ";", fixed = TRUE)[[1]]
        pairs <- pairs[nzchar(pairs)]
        keys  <- sub("=.*",     "", pairs)
        vals  <- sub("^[^=]+=", "", pairs)
        stats::setNames(as.list(vals), keys)
    })
    data.table::rbindlist(parsed, fill = TRUE)
}

#' Load and parse the project GFF as a gene coordinate table
#'
#' Reads the GFF file, filters to the configured feature type (default "mRNA"),
#' extracts gene_id from the Parent= or ID= attribute, and parses all attribute
#' fields into columns. Cached per project to avoid re-loading the full GFF on
#' every reactive trigger.
#'
#' @param project project name
#' @param config parsed YAML config (from project_data()$config)
#' @return data.table with gene_id, chr, start, end, + attribute columns.
#'   Empty data.table if GFF not found or not readable.
#' @noRd
load_gff_genes <- function(project, config) {
    cache_key <- paste0("gff_genes_", project)
    load_cached(cache_key, function() {
        inp_dir  <- config_get(config, "Input", "dir",     default = "data")
        gff_rel  <- config_get(config, "Input", "gff",     default = "")
        feat     <- config_get(config, "GFF",   "feature", default = "mRNA")

        if (is.null(gff_rel) || gff_rel == "") return(data.table::data.table())

        gff_abs <- file.path(get_pipeline_path(), inp_dir, gff_rel)
        if (!file.exists(gff_abs)) return(data.table::data.table())

        tryCatch({
            # Read GFF, skip comment lines, filter to target feature
            dt <- data.table::fread(
                cmd    = paste0("grep -v '^#' ", shQuote(gff_abs)),
                header = FALSE,
                col.names = c("chr", "source", "feature", "start", "end",
                              "score", "strand", "phase", "attributes"),
                colClasses = c("character", "character", "character",
                               "integer", "integer",
                               "character", "character", "character", "character")
            )
            dt <- dt[feature == feat]
            if (nrow(dt) == 0) return(data.table::data.table())

            dt[, chr   := sub("^chr", "", as.character(chr))]
            dt[, start := as.integer(start)]
            dt[, end   := as.integer(end)]

            # Extract gene_id: prefer Parent=, fall back to ID=
            dt[, gene_id := {
                id <- regmatches(attributes,
                    regexpr("(?<=Parent=)[^;]+", attributes, perl = TRUE))
                na_idx <- !nzchar(id)
                if (any(na_idx)) {
                    id2 <- regmatches(attributes,
                        regexpr("(?<=ID=)[^;]+", attributes, perl = TRUE))
                    id[na_idx] <- id2[na_idx]
                }
                # Strip isoform suffixes (.1, .2, _exon_1 etc.)
                sub("\\.[0-9]+(_exon_[0-9]+)?$", "", id)
            }]

            # Parse attributes into individual columns
            attr_dt <- parse_gff_attributes(dt$attributes)
            # Drop ID/Parent — redundant with gene_id
            drop_cols <- intersect(c("ID", "Parent"), names(attr_dt))
            if (length(drop_cols) > 0) attr_dt[, (drop_cols) := NULL]
            cbind(dt[, .(gene_id, chr, start, end)], attr_dt)
        }, error = function(e) data.table::data.table())
    })
}

#' Load all per-method sig SNPs for a module (for interactive Shiny combine)
#'
#' Reads all {method}_pvalues_K*_sig_snps_*.tsv files found under
#' {module}/tables/methods/, returning a named list (method → data.table).
#' Each data.table has columns: SNPID, chr, pos, pvalue, method, trait.
#' Returns an empty list when no files exist (old pipeline run / fallback).
#' @noRd
load_all_method_sigsnps <- function(project, module = MOD_GEA, k = NA) {
    files <- find_method_sigsnps_files(project, module)
    if (length(files) == 0) return(list())
    if (!is.na(k)) {
        k_pattern <- paste0("_K", k, "_")
        files <- files[grepl(k_pattern, basename(files), fixed = TRUE)]
    }
    if (length(files) == 0) return(list())

    result <- list()
    for (f in files) {
        # Extract method from path: .../tables/methods/{method}/...
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        methods_idx <- which(parts == "methods")
        if (length(methods_idx) == 0) next
        method_name <- parts[methods_idx + 1L]

        dt <- tryCatch(
            data.table::fread(f, sep = "\t", header = TRUE,
                              colClasses = c(chr    = "character",
                                             SNPID  = "character",
                                             method = "character",
                                             trait  = "character")),
            error = function(e) data.table::data.table()
        )
        if (nrow(dt) == 0) next
        keep <- intersect(c("SNPID", "chr", "pos", "pvalue", "method", "trait"), names(dt))
        if (length(keep) < 5L) next
        result[[method_name]] <- dt[, ..keep]
    }
    result
}

#' Load all per-method WZA sig-windows (WZA regime counterpart of load_all_method_sigsnps)
#' @noRd
load_all_method_wza_sigwindows <- function(project, module = MOD_GEA, k = NA) {
    glob <- mod_path(project, module, "tables", "methods", "*",
                     "*_wza_K*_sig_windows_*.tsv")
    files <- Sys.glob(glob)
    if (length(files) == 0) return(list())
    if (!is.na(k)) {
        k_pattern <- paste0("_K", k, "_")
        files <- files[grepl(k_pattern, basename(files), fixed = TRUE)]
    }
    if (length(files) == 0) return(list())

    result <- list()
    for (f in files) {
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        methods_idx <- which(parts == "methods")
        if (length(methods_idx) == 0) next
        method_name <- parts[methods_idx + 1L]

        dt <- tryCatch(
            data.table::fread(f, sep = "\t", header = TRUE,
                              colClasses = c(chr    = "character",
                                             SNPID  = "character",
                                             method = "character",
                                             trait  = "character")),
            error = function(e) data.table::data.table()
        )
        if (nrow(dt) == 0) next
        keep <- intersect(c("SNPID", "chr", "pos", "pvalue", "method", "trait"), names(dt))
        if (length(keep) < 5L) next
        result[[method_name]] <- dt[, ..keep]
    }
    result
}

#' Load phenotype missing summary (GWAS mode)
#' @noRd
load_pheno_missing_summary <- function(project) {
    p <- pheno_missing_summary_path(project)
    if (!file.exists(p)) return(data.table::data.table())
    tryCatch(
        data.table::fread(p, sep = "\t", header = TRUE,
                          colClasses = c(trait = "character")),
        error = function(e) data.table::data.table()
    )
}

#' Load trait names from any pvalue TSV (column names after chr/pos/SNPID)
#' Used to detect traits that have pvalue output but zero significant SNPs.
#' Returns character(0) when no pvalue files exist.
#' @noRd
load_all_trait_names <- function(project, module = MOD_GEA) {
    f <- find_pvalue_tsv(project, module)
    if (!nzchar(f) || !file.exists(f)) return(character(0))
    tryCatch({
        hdr <- data.table::fread(f, nrows = 0L)
        fixed <- c("chr", "pos", "SNPID")
        setdiff(names(hdr), fixed)
    }, error = function(e) character(0))
}

#' Load all per-method sig SNPs for the Overlapping tab
#'
#' Combines per-method sig SNPs from both MOD_GEA and MOD_GWAS.
#' GEA traits (bio_*) and GWAS traits (phenotype columns) are naturally distinct,
#' so same method names across modules are merged without prefixing.
#'
#' @param project Character: project name
#' @param k       Optional K value to filter by (NA = all)
#' @return Named list of data.tables (method -> dt), same format as load_all_method_sigsnps()
#' @noRd
load_all_overlap_method_sigsnps <- function(project, k = NA) {
    gea  <- load_all_method_sigsnps(project, MOD_GEA,  k)
    gwas <- load_all_method_sigsnps(project, MOD_GWAS, k)

    all_methods <- union(names(gea), names(gwas))
    result <- list()
    for (m in all_methods) {
        parts <- Filter(
            function(x) !is.null(x) && nrow(x) > 0,
            list(gea[[m]], gwas[[m]])
        )
        if (length(parts) > 0)
            result[[m]] <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
    }
    result
}

#' Load full p-value TSVs for all methods (wide format: one column per trait).
#'
#' Returns a named list (method → data.table with SNPID, chr, pos, trait1..N).
#' Excludes _sig_snps_, _wza_, and _block_ files.
#' @noRd
load_all_method_pvalues <- function(project, module = MOD_GEA, k = NA) {
    glob <- mod_path(project, module, "tables", "methods", "*", "*_pvalues_K*.tsv")
    files <- Sys.glob(glob)
    # Exclude downstream outputs (sig_snps, wza, block)
    files <- grep("_sig_snps_|_wza_|_block_", files, value = TRUE, invert = TRUE)
    if (length(files) == 0) return(list())
    if (!is.na(k)) {
        k_pattern <- paste0("_K", k, ".tsv")
        files <- files[endsWith(basename(files), k_pattern)]
    }
    if (length(files) == 0) return(list())

    result <- list()
    for (f in files) {
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        methods_idx <- which(parts == "methods")
        if (length(methods_idx) == 0) next
        method_name <- parts[methods_idx + 1L]

        dt <- tryCatch(
            data.table::fread(f, sep = "\t", header = TRUE,
                              colClasses = c(chr = "character", SNPID = "character")),
            error = function(e) data.table::data.table()
        )
        if (nrow(dt) == 0) next
        result[[method_name]] <- dt
    }
    result
}

#' Load full WZA p-value TSVs for all methods (wide format: one column per trait).
#' @noRd
load_all_method_wza_pvalues <- function(project, module = MOD_GEA, k = NA) {
    glob <- mod_path(project, module, "tables", "methods", "*", "*_wza_K*.tsv")
    files <- Sys.glob(glob)
    files <- grep("_sig_windows_", files, value = TRUE, invert = TRUE)
    if (length(files) == 0) return(list())
    if (!is.na(k)) {
        k_pattern <- paste0("_K", k, ".tsv")
        files <- files[endsWith(basename(files), k_pattern)]
    }
    if (length(files) == 0) return(list())

    result <- list()
    for (f in files) {
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        methods_idx <- which(parts == "methods")
        if (length(methods_idx) == 0) next
        method_name <- parts[methods_idx + 1L]

        dt <- tryCatch(
            data.table::fread(f, sep = "\t", header = TRUE,
                              colClasses = c(chr = "character", SNPID = "character")),
            error = function(e) data.table::data.table()
        )
        if (nrow(dt) == 0) next
        result[[method_name]] <- dt
    }
    result
}

#' Compute a fingerprint string for a set of files (mtime + size).
#'
#' Used to detect when source p-value files have been regenerated by the pipeline,
#' so disk caches keyed on this fingerprint automatically miss after a Run.
#' @noRd
pvalues_file_fingerprint <- function(project, module, k, regime = "snp") {
    suffix_pat <- if (regime == "wza") "*_wza_K*.tsv" else "*_pvalues_K*.tsv"
    glob_pat   <- mod_path(project, module, "tables", "methods", "*", suffix_pat)
    files      <- Sys.glob(glob_pat)
    files      <- grep("_sig_snps_|_block_|_sig_windows_", files,
                       value = TRUE, invert = TRUE)
    if (!is.na(k)) {
        k_pat  <- paste0("_K", k, ".tsv")
        files  <- files[endsWith(basename(files), k_pat)]
    }
    if (length(files) == 0) return("empty")
    fi <- file.info(files[order(files)])  # stable sort so fingerprint is deterministic
    digest::digest(list(mtime = fi$mtime, size = fi$size), algo = "md5")
}

#' Compute significant SNPs interactively using a user-chosen threshold.
#'
#' Wraps compute_pval_threshold() (sourced from pipeline utils) per method × trait.
#' Disk-caches results keyed by (module, k, regime, type, value, file_fingerprint)
#' so qval (slow) is paid once per unique combination; bonf/top/custom are sub-second.
#' Including the file fingerprint ensures the cache is automatically bypassed when the
#' pipeline regenerates p-value files with a different trait set.
#'
#' @param pvalues_list Named list of wide data.tables (method → dt) from
#'   load_all_method_pvalues() or load_all_method_wza_pvalues().
#' @param type   Character: "bonf" | "qval" | "top" | "custom"
#' @param value  Numeric threshold value (alpha for bonf/qval, N for top, raw p for custom)
#' @param k      K-best value (for cache key)
#' @param regime Character: "snp" or "wza"
#' @param project Project name (for cache dir)
#' @param module  Pipeline module (for cache dir)
#' @return Named list (method → data.table with SNPID/chr/pos/pvalue/method/trait),
#'   same shape as load_all_method_sigsnps().
#' @noRd
compute_method_sigsnps_cached <- function(pvalues_list, type, value,
                                           k, regime, project, module) {
    if (length(pvalues_list) == 0 || is.null(pvalues_list)) return(list())

    # Validate that compute_pval_threshold is available (sourced above)
    if (!exists("compute_pval_threshold", mode = "function")) {
        warning("compute_pval_threshold not available — returning empty sig set")
        return(list())
    }

    # Disk-cache key: (module, k, regime, type, value) + source file fingerprint.
    # The fingerprint (md5 of mtime+size of source pvalue files) ensures that when the
    # pipeline regenerates files with a new trait set, the stale cache entry is bypassed.
    fp <- pvalues_file_fingerprint(project, module, k, regime)
    cache_key  <- list(module = module, k = k, regime = regime,
                       type = type, value = as.numeric(value), fp = fp)
    cache_hash <- digest::digest(cache_key, algo = "md5")
    cache_dir  <- interactive_sigsnps_dir(project, module)
    cache_file <- file.path(cache_dir, paste0(cache_hash, ".tsv"))

    if (file.exists(cache_file)) {
        flat <- tryCatch(
            data.table::fread(cache_file, sep = "\t", header = TRUE,
                              colClasses = c(chr    = "character",
                                             SNPID  = "character",
                                             method = "character",
                                             trait  = "character")),
            error = function(e) NULL
        )
        if (!is.null(flat) && nrow(flat) > 0) {
            return(split(flat, by = "method", keep.by = TRUE))
        }
    }

    # Compute: iterate methods and traits
    all_rows <- list()
    for (m in names(pvalues_list)) {
        pv_dt <- pvalues_list[[m]]
        if (nrow(pv_dt) == 0) next

        fixed_cols <- c("SNPID", "chr", "pos")
        # WZA tables also have n_snps, mean_maf — skip non-trait cols
        trait_cols <- setdiff(names(pv_dt), c(fixed_cols, "n_snps", "mean_maf"))
        if (length(trait_cols) == 0) next

        for (tr in trait_cols) {
            if (!(tr %in% names(pv_dt))) next
            pvec <- pv_dt[[tr]]

            thresh_result <- tryCatch(
                compute_pval_threshold(pvec, type, as.numeric(value)),
                error = function(e) list(status = "error", threshold = NA_real_)
            )
            if (thresh_result$status != "ok" || is.na(thresh_result$threshold)) next

            mask <- pvec < thresh_result$threshold & !is.na(pvec)
            if (!any(mask)) next

            all_rows[[length(all_rows) + 1L]] <- pv_dt[mask, .(
                SNPID  = SNPID,
                chr    = chr,
                pos    = pos,
                pvalue = .SD[[tr]],
                method = m,
                trait  = tr
            ), .SDcols = tr]
        }
    }

    if (length(all_rows) == 0) {
        flat <- data.table::data.table(
            SNPID  = character(),
            chr    = character(),
            pos    = integer(),
            pvalue = numeric(),
            method = character(),
            trait  = character()
        )
    } else {
        flat <- data.table::rbindlist(all_rows)
    }

    # Save to disk cache
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(flat, cache_file, sep = "\t")

    if (nrow(flat) == 0) return(list())
    split(flat, by = "method", keep.by = TRUE)
}

#' Compute significance threshold per trait × method (no IO, no cache)
#'
#' Returns a named numeric vector keyed "trait::method" -> threshold p-value.
#' NA when too few tests or compute_pval_threshold status != "ok".
#' Works for all modes:
#'   bonf -> 0.05 / n_tested
#'   qval -> max p still passing FDR
#'   top  -> Nth smallest p
#' For WZA regime, pass the WZA window pvalue list — denominator is n_windows.
#'
#' @param pvalues_list Named list: method -> data.table with SNPID/chr/pos + trait cols
#' @param type         Adjustment type: "bonf", "qval", "top"
#' @param value        Significance level / FDR / top-N (numeric or character)
#' @return Named numeric vector "trait::method" -> threshold (NA if unavailable)
#' @noRd
compute_method_thresholds <- function(pvalues_list, type, value) {
    if (length(pvalues_list) == 0 || is.null(pvalues_list)) return(list())
    fixed_cols <- c("SNPID", "chr", "pos", "n_snps", "mean_maf")
    out <- list()
    for (m in names(pvalues_list)) {
        pv_dt <- pvalues_list[[m]]
        if (is.null(pv_dt) || nrow(pv_dt) == 0) next
        trait_cols <- setdiff(names(pv_dt), fixed_cols)
        if (length(trait_cols) == 0) next
        for (tr in trait_cols) {
            pvec <- pv_dt[[tr]]
            result <- tryCatch(
                compute_pval_threshold(pvec, type, as.numeric(value)),
                error = function(e) list(status = "error", threshold = NA_real_)
            )
            key <- paste0(tr, "::", m)
            out[[key]] <- if (identical(result$status, "ok")) result$threshold else NA_real_
        }
    }
    out
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
#'
#' @param cache_key  Base cache key string.
#' @param loader_fn  Zero-argument function that loads the data.
#' @param fingerprint Optional string appended to the key. Pass a file's mtime or
#'   a digest of mtime+size so the cache is bypassed when the file is regenerated
#'   by the pipeline without changing its path.
#' @noRd
load_cached <- function(cache_key, loader_fn, fingerprint = NULL) {
    cache <- get_app_cache()
    full_key <- if (!is.null(fingerprint)) paste0(cache_key, "_", fingerprint) else cache_key
    safe_key <- gsub("[^a-z0-9]", "", tolower(full_key))
    cached <- cache$get(safe_key)
    if (!cachem::is.key_missing(cached)) return(cached)
    val <- loader_fn()
    cache$set(safe_key, val)
    val
}
