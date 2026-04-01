# On-the-fly region computation for interactive Shiny exploration.
# Ports single-linkage clustering from scripts/create_regions.R using data.table
# (avoids GenomicRanges dependency in Shiny).

#' Compute all regions from sig SNPs using single-linkage clustering
#'
#' Two SNPs merge into the same region when their genomic gap <= distance.
#' Equivalent to create_regions.R's extend-by-distance/2 + GenomicRanges::reduce().
#'
#' @param sig_snps data.table with columns SNPID, chr, pos, pvalue, method, trait
#' @param distance integer. Maximum gap between neighboring SNPs to be in the
#'   same region (i.e. REGION_DISTANCE config value).
#' @return data.table with region_id, chr, start, end, length, snp_count,
#'   snp_ids, traits, methods, min_pvalue. Empty data.table if no SNPs.
#' @noRd
compute_all_regions <- function(sig_snps, distance = 2000000L) {
    distance <- as.integer(distance)

    if (is.null(sig_snps) || nrow(sig_snps) == 0) return(.empty_regions())

    # Unique positions (one row per SNP)
    unique_snps <- unique(sig_snps[, .(SNPID, chr = as.character(chr),
                                       pos = as.integer(pos),
                                       min_pvalue = pvalue)])
    unique_snps[, min_pvalue := min(min_pvalue, na.rm = TRUE), by = "SNPID"]
    unique_snps <- unique(unique_snps, by = "SNPID")

    data.table::setorder(unique_snps, chr, pos)

    # Single-linkage clustering via gap detection
    unique_snps[, gap     := c(Inf, diff(pos)), by = "chr"]
    unique_snps[, cluster := cumsum(gap > distance), by = "chr"]

    # Aggregate: region boundaries = min(pos) - distance to max(pos) + distance
    regions <- unique_snps[, .(
        chr       = chr[1],
        start     = pmax(1L, min(pos) - distance),
        end       = max(pos) + distance,
        snp_count = .N,
        snp_ids   = paste(SNPID, collapse = ","),
        min_pvalue = min(min_pvalue, na.rm = TRUE)
    ), by = c("chr", "cluster")]

    regions[, region_id := paste0(chr, "_", start, "-", end)]
    regions[, length    := end - start]
    regions[, cluster   := NULL]
    data.table::setorder(regions, chr, start)

    # Join back to full sig_snps to collect traits and methods per region
    # Use foverlaps: match each SNP to its enclosing region
    snp_pos <- unique(sig_snps[, .(SNPID, chr = as.character(chr),
                                    pos = as.integer(pos), method, trait)])
    snp_pos[, s := pos][, e := pos]
    reg_ov <- data.table::copy(regions[, .(region_id, chr, start, end)])
    data.table::setkey(reg_ov, chr, start, end)
    data.table::setkey(snp_pos, chr, s, e)
    ov <- data.table::foverlaps(snp_pos, reg_ov,
        by.x = c("chr", "s", "e"),
        by.y = c("chr", "start", "end"),
        type = "within", mult = "first", nomatch = NA)

    # Aggregate traits and methods per region
    trait_methods <- ov[!is.na(region_id), .(
        traits  = paste(sort(unique(trait)),  collapse = ","),
        methods = paste(sort(unique(method)), collapse = ",")
    ), by = "region_id"]

    regions <- trait_methods[regions, on = "region_id"]

    # Select and order columns
    data.table::setcolorder(regions,
        c("region_id", "chr", "start", "end", "length",
          "snp_count", "snp_ids", "traits", "methods", "min_pvalue"))
    regions
}

#' Find genes overlapping a region (data.table-based, no GenomicRanges)
#'
#' @param gff_genes data.table from load_gff_genes(): gene_id, chr, start, end, ...
#' @param region_row single-row data.table with chr, start, end
#' @param promoter_length integer. Extend region upstream to capture promoters.
#' @return data.table of genes overlapping [start - promoter_length, end].
#'   Empty data.table if gff_genes is NULL or no overlap found.
#' @noRd
find_genes_in_region <- function(gff_genes, region_row, promoter_length = 10000L) {
    if (is.null(gff_genes) || nrow(gff_genes) == 0) return(data.table::data.table())
    if (is.null(region_row) || nrow(region_row) == 0) return(data.table::data.table())

    promoter_length <- as.integer(promoter_length)
    chr_val <- as.character(region_row$chr[1])
    r_start <- pmax(1L, as.integer(region_row$start[1]) - promoter_length)
    r_end   <- as.integer(region_row$end[1])

    # Extend region window
    query <- data.table::data.table(
        chr   = chr_val,
        start = r_start,
        end   = r_end
    )

    genes_chr <- data.table::copy(gff_genes[chr == chr_val])
    if (nrow(genes_chr) == 0) return(data.table::data.table())

    data.table::setkey(genes_chr, chr, start, end)
    data.table::setkey(query, chr, start, end)

    ov <- data.table::foverlaps(query, genes_chr,
        by.x = c("chr", "start", "end"),
        by.y = c("chr", "start", "end"),
        type = "any", mult = "all", nomatch = NULL)

    if (is.null(ov) || nrow(ov) == 0) return(data.table::data.table())

    # Remove query columns (i.start, i.end), keep gene columns
    drop_cols <- intersect(c("i.start", "i.end"), names(ov))
    if (length(drop_cols) > 0) ov[, (drop_cols) := NULL]
    ov
}

#' Run GO enrichment for a single region via processx subprocess
#'
#' Calls scripts/run_enrichment.R then scripts/plot_enrichment.R.
#' Results written to a session-specific temp directory.
#'
#' @param genes_dt data.table of genes for this region (from find_genes_in_region())
#' @param region_id character combined region_id (e.g. "1_1000000-2000000")
#' @param trait character trait name (for per-trait enrichment)
#' @param project_data reactive project data bundle (name, config)
#' @return list(table, dotplot_path, emapplot_path) or NULL on failure
#' @noRd
run_enrichment_subprocess <- function(genes_dt, region_id, trait, project_data) {
    pd <- project_data

    go_field   <- config_get(pd$config, "association", "go_field",  default = NULL)
    gff_feat   <- config_get(pd$config, "gff", "feature",           default = "mRNA")
    top_terms  <- config_get(pd$config, "enrichment", "top_terms",  default = 10L)
    plot_w     <- config_get(pd$config, "enrichment", "plot_width",  default = 10L)
    plot_h     <- config_get(pd$config, "enrichment", "plot_height", default = 8L)
    cnet_label <- config_get(pd$config, "enrichment", "cnet_label", default = "gene_id")

    if (is.null(go_field) || go_field == "" || go_field == "NULL") {
        return(list(error = "GO field not configured (association.go_field)"))
    }

    gff <- .resolve_gff_path(pd)
    if (!file_ok(gff)) {
        return(list(error = paste0("GFF not found: ", gff)))
    }

    if (is.null(genes_dt) || nrow(genes_dt) == 0) {
        return(list(error = "No genes in region"))
    }

    pipeline_path <- get_pipeline_path()

    # Temp session directory for this region's enrichment output
    tmp_base  <- file.path(tempdir(), "adaptogene_enrichment",
                           paste0(pd$name, "_", gsub("[^a-z0-9]", "", tolower(region_id)),
                                  "_", trait))
    tmp_genes   <- file.path(tmp_base, "genes.tsv")
    tmp_tables  <- file.path(tmp_base, "tables")
    tmp_inter   <- file.path(tmp_base, "intermediate")
    tmp_plots   <- file.path(tmp_base, "plots")

    dir.create(tmp_tables, recursive = TRUE, showWarnings = FALSE)
    dir.create(tmp_inter,  recursive = TRUE, showWarnings = FALSE)
    dir.create(tmp_plots,  recursive = TRUE, showWarnings = FALSE)

    # Write per-trait genes file (run_enrichment.R expects region_id with trait suffix)
    per_trait_id <- paste0(region_id, "_", trait)
    genes_out <- data.table::copy(genes_dt)
    genes_out[, region_id := per_trait_id]
    data.table::fwrite(genes_out, tmp_genes, sep = "\t")

    # ── Step 1: Run enrichment ─────────────────────────────────────────────────
    res_enrich <- tryCatch(
        processx::run(
            command = "Rscript",
            args    = c(
                file.path(pipeline_path, "scripts", "run_enrichment.R"),
                tmp_genes,
                gff,
                go_field,
                gff_feat,
                tmp_tables,
                tmp_inter
            ),
            wd      = pipeline_path,
            timeout = 120
        ),
        error = function(e) list(status = 1L, stderr = conditionMessage(e))
    )

    if (res_enrich$status != 0L) {
        return(list(error = paste0("Enrichment failed: ",
                                   utils::tail(strsplit(res_enrich$stderr, "\n")[[1]], 3))))
    }

    # Check if TSV was produced
    tsv_path <- file.path(tmp_tables, trait,
                          paste0("Region_", per_trait_id, "_enrichment.tsv"))
    if (!file_ok(tsv_path)) {
        return(list(error = "No enrichment results (no significant GO terms)"))
    }

    enrich_tbl <- tryCatch(
        data.table::fread(tsv_path, sep = "\t", header = TRUE),
        error = function(e) NULL
    )
    if (is.null(enrich_tbl) || nrow(enrich_tbl) == 0) {
        return(list(error = "No enrichment results (no significant GO terms)"))
    }

    # ── Step 2: Generate plots ─────────────────────────────────────────────────
    res_plots <- tryCatch(
        processx::run(
            command = "Rscript",
            args    = c(
                file.path(pipeline_path, "scripts", "plot_enrichment.R"),
                tmp_inter,
                tmp_plots,
                as.character(as.integer(top_terms)),
                as.character(as.numeric(plot_w)),
                as.character(as.numeric(plot_h)),
                cnet_label,
                tmp_genes,
                "0"  # top_plot_regions = 0 (all)
            ),
            wd      = pipeline_path,
            timeout = 120
        ),
        error = function(e) list(status = 1L, stderr = conditionMessage(e))
    )

    dotplot_path  <- file.path(tmp_plots, trait,
                               paste0("Region_", per_trait_id, "_dotplot.png"))
    emapplot_path <- file.path(tmp_plots, trait,
                               paste0("Region_", per_trait_id, "_emapplot.png"))

    list(
        table         = enrich_tbl,
        dotplot_path  = if (file_ok(dotplot_path))  dotplot_path  else NULL,
        emapplot_path = if (file_ok(emapplot_path)) emapplot_path else NULL
    )
}

#' Run regionplot for a single region via processx subprocess
#'
#' Calls scripts/plot_regionplot.R with CUSTOM_REGION set.
#'
#' @param region_row single-row data.table with chr, start, end
#' @param trait character trait name
#' @param project_data project data bundle
#' @param module character pipeline module (MOD_ASSOC or MOD_PHENO)
#' @return character path to generated PNG, or NULL on failure
#' @noRd
run_regionplot_subprocess <- function(region_row, trait, project_data, module = MOD_ASSOC) {
    pd <- project_data

    # Convert GFF on-the-fly if the topr annotation doesn't exist yet
    gff_topr_p <- .ensure_gff_topr(pd)
    if (!file_ok(gff_topr_p)) {
        return(list(error = "Could not find or generate topr gene annotation. Check GFF config."))
    }

    # Build ASSOC_TABLES: method:adjust:filepath for each available method sig SNPs file
    method_files <- find_method_sigsnps_files(pd$name, module)
    if (length(method_files) == 0) {
        return(list(error = "No method result files found. Run the association mode first."))
    }

    # Build method:adjust_str:pvalues_file triples from the paths
    assoc_tables <- .build_assoc_tables_arg(method_files)
    if (is.null(assoc_tables)) {
        return(list(error = "Could not build association table arguments."))
    }

    custom_region <- paste0(region_row$chr[1], ":",
                             region_row$start[1], "-",
                             region_row$end[1])

    genes_highlight <- config_get(pd$config, "regionplot", "genes",
                                  default = "NULL") %||% "NULL"
    if (length(genes_highlight) > 1) genes_highlight <- paste(genes_highlight, collapse = "|")

    pipeline_path <- get_pipeline_path()

    tmp_base <- file.path(tempdir(), "adaptogene_regionplot",
                          paste0(pd$name, "_",
                                 gsub("[^a-z0-9]", "", tolower(custom_region)),
                                 "_", trait))
    dir.create(tmp_base, recursive = TRUE, showWarnings = FALSE)

    res <- tryCatch(
        processx::run(
            command = "Rscript",
            args    = c(
                file.path(pipeline_path, "scripts", "plot_regionplot.R"),
                "NULL",          # REGIONS_FILE (not needed for custom region)
                gff_topr_p,
                assoc_tables,
                "0",             # TOP_REGIONS
                genes_highlight,
                paste0(tmp_base, "/"),
                custom_region,
                trait,           # CUSTOM_TRAITS
                "NULL"           # CUSTOM_METHODS (all)
            ),
            wd      = pipeline_path,
            timeout = 180
        ),
        error = function(e) list(status = 1L, stderr = conditionMessage(e))
    )

    if (res$status != 0L) {
        stderr_tail <- paste(utils::tail(strsplit(res$stderr %||% "", "\n")[[1]], 5), collapse = "\n")
        return(list(error = paste0("Regionplot failed:\n", stderr_tail)))
    }

    # Find generated PNG in tmp_base
    pngs <- Sys.glob(file.path(tmp_base, "*.png"))
    if (length(pngs) == 0) {
        return(list(error = "No plot was generated. Check that the region contains SNPs in the p-value files."))
    }
    list(path = pngs[1])
}

# ── Async subprocess launchers ────────────────────────────────────────────────

#' Launch GO enrichment as a non-blocking background process.
#'
#' Returns a handle list immediately; the caller polls \code{handle$process$is_alive()}.
#' On completion, pass the handle to \code{.collect_enrichment_result()}.
#'
#' @return list(process, log_file, type, tmp_tables, tmp_plots, per_trait_id, trait)
#'   or list(error=) if pre-flight validation fails.
#' @noRd
launch_enrichment_subprocess <- function(genes_dt, region_id, trait, project_data) {
    pd <- project_data

    go_field   <- config_get(pd$config, "association", "go_field",  default = NULL)
    gff_feat   <- config_get(pd$config, "gff", "feature",           default = "mRNA")
    top_terms  <- config_get(pd$config, "enrichment", "top_terms",  default = 10L)
    plot_w     <- config_get(pd$config, "enrichment", "plot_width",  default = 10L)
    plot_h     <- config_get(pd$config, "enrichment", "plot_height", default = 8L)
    cnet_label <- config_get(pd$config, "enrichment", "cnet_label", default = "gene_id")

    if (is.null(go_field) || go_field == "" || go_field == "NULL")
        return(list(error = "GO field not configured (association.go_field)"))

    gff <- .resolve_gff_path(pd)
    if (!file_ok(gff))
        return(list(error = paste0("GFF not found: ", gff)))

    if (is.null(genes_dt) || nrow(genes_dt) == 0)
        return(list(error = "No genes in region"))

    pipeline_path <- get_pipeline_path()
    tmp_base  <- file.path(tempdir(), "adaptogene_enrichment",
                           paste0(pd$name, "_", gsub("[^a-z0-9]", "", tolower(region_id)),
                                  "_", trait))
    tmp_genes  <- file.path(tmp_base, "genes.tsv")
    tmp_tables <- file.path(tmp_base, "tables")
    tmp_inter  <- file.path(tmp_base, "intermediate")
    tmp_plots  <- file.path(tmp_base, "plots")
    log_file   <- file.path(tmp_base, "enrichment.log")

    dir.create(tmp_tables, recursive = TRUE, showWarnings = FALSE)
    dir.create(tmp_inter,  recursive = TRUE, showWarnings = FALSE)
    dir.create(tmp_plots,  recursive = TRUE, showWarnings = FALSE)

    per_trait_id <- paste0(region_id, "_", trait)
    genes_out <- data.table::copy(genes_dt)
    genes_out[, region_id := per_trait_id]
    data.table::fwrite(genes_out, tmp_genes, sep = "\t")

    enrich_args <- paste(sapply(c(
        file.path(pipeline_path, "scripts", "run_enrichment.R"),
        tmp_genes, gff, go_field, gff_feat, tmp_tables, tmp_inter
    ), shQuote), collapse = " ")

    plot_args <- paste(sapply(c(
        file.path(pipeline_path, "scripts", "plot_enrichment.R"),
        tmp_inter, tmp_plots,
        as.character(as.integer(top_terms)),
        as.character(as.numeric(plot_w)),
        as.character(as.numeric(plot_h)),
        cnet_label, tmp_genes, "0"
    ), shQuote), collapse = " ")

    # Chain both steps; log file captures stdout+stderr of both
    cmd <- paste0(
        "Rscript ", enrich_args,
        " >> ", shQuote(log_file), " 2>&1",
        " && Rscript ", plot_args,
        " >> ", shQuote(log_file), " 2>&1"
    )

    proc <- tryCatch(
        processx::process$new("bash", c("-c", cmd), wd = pipeline_path),
        error = function(e) NULL
    )
    if (is.null(proc))
        return(list(error = "Failed to launch enrichment process"))

    list(
        process      = proc,
        log_file     = log_file,
        type         = "enrichment",
        tmp_tables   = tmp_tables,
        tmp_plots    = tmp_plots,
        per_trait_id = per_trait_id,
        trait        = trait
    )
}

#' Collect enrichment result after the background process has finished.
#' @noRd
.collect_enrichment_result <- function(handle) {
    tsv_path <- file.path(handle$tmp_tables, handle$trait,
                          paste0("Region_", handle$per_trait_id, "_enrichment.tsv"))
    if (!file_ok(tsv_path))
        return(list(error = "No enrichment results (no significant GO terms)"))

    enrich_tbl <- tryCatch(
        data.table::fread(tsv_path, sep = "\t", header = TRUE),
        error = function(e) NULL
    )
    if (is.null(enrich_tbl) || nrow(enrich_tbl) == 0)
        return(list(error = "No enrichment results (no significant GO terms)"))

    dotplot_path  <- file.path(handle$tmp_plots, handle$trait,
                               paste0("Region_", handle$per_trait_id, "_dotplot.png"))
    emapplot_path <- file.path(handle$tmp_plots, handle$trait,
                               paste0("Region_", handle$per_trait_id, "_emapplot.png"))

    list(
        table         = enrich_tbl,
        dotplot_path  = if (file_ok(dotplot_path))  dotplot_path  else NULL,
        emapplot_path = if (file_ok(emapplot_path)) emapplot_path else NULL
    )
}

#' Build trait:method1|method2,trait2:method3 arg from region sig SNPs.
#' @noRd
.build_trait_method_arg <- function(region_snps) {
    if (is.null(region_snps) || nrow(region_snps) == 0) return(NULL)
    pairs <- unique(region_snps[, .(trait, method)])
    if (nrow(pairs) == 0) return(NULL)
    # Group methods by trait
    by_trait <- split(pairs$method, pairs$trait)
    paste(mapply(function(trait, methods) {
        paste0(trait, ":", paste(sort(unique(methods)), collapse = "|"))
    }, names(by_trait), by_trait), collapse = ",")
}

#' Launch regionplot as a non-blocking background process.
#' @param region_snps data.table with columns trait, method — sig SNPs in the region
#' @return list(process, log_file, type, tmp_base) or list(error=)
#' @noRd
launch_regionplot_subprocess <- function(region_row, region_snps, project_data, module = MOD_ASSOC) {
    pd <- project_data

    gff_topr_p <- .ensure_gff_topr(pd)
    if (!file_ok(gff_topr_p))
        return(list(error = "Could not find or generate topr gene annotation. Check GFF config."))

    method_files <- find_method_sigsnps_files(pd$name, module)
    if (length(method_files) == 0)
        return(list(error = "No method result files found. Run the association mode first."))
    # Filter to k_best only — glob returns files for all K values, causing each method
    # to appear twice in the regionplot legend (once per K) with different colors
    k_best <- pd$k_best
    if (!is.na(k_best)) {
        k_pattern <- paste0("_K", k_best, "_")
        method_files <- method_files[grepl(k_pattern, basename(method_files), fixed = TRUE)]
    }

    assoc_tables <- .build_assoc_tables_arg(method_files)
    if (is.null(assoc_tables))
        return(list(error = "Could not build association table arguments."))

    # Build structured trait:method pairs from region sig SNPs
    custom_traits <- .build_trait_method_arg(region_snps)
    if (is.null(custom_traits))
        return(list(error = "No significant trait/method pairs found in this region."))

    custom_region <- paste0(region_row$chr[1], ":", region_row$start[1], "-", region_row$end[1])
    genes_highlight <- config_get(pd$config, "regionplot", "genes", default = "NULL") %||% "NULL"
    if (length(genes_highlight) > 1) genes_highlight <- paste(genes_highlight, collapse = "|")

    pipeline_path <- get_pipeline_path()
    tmp_base <- file.path(tempdir(), "adaptogene_regionplot",
                          paste0(pd$name, "_",
                                 gsub("[^a-z0-9]", "", tolower(custom_region)),
                                 "_multi"))
    dir.create(tmp_base, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(tmp_base, "regionplot.log")

    rplot_args <- paste(sapply(c(
        file.path(pipeline_path, "scripts", "plot_regionplot.R"),
        "NULL", gff_topr_p, assoc_tables, "0",
        genes_highlight, paste0(tmp_base, "/"),
        custom_region, custom_traits, "NULL"
    ), shQuote), collapse = " ")

    cmd <- paste0("Rscript ", rplot_args, " >> ", shQuote(log_file), " 2>&1")

    proc <- tryCatch(
        processx::process$new("bash", c("-c", cmd), wd = pipeline_path),
        error = function(e) NULL
    )
    if (is.null(proc))
        return(list(error = "Failed to launch regionplot process"))

    list(
        process  = proc,
        log_file = log_file,
        type     = "regionplot",
        tmp_base = tmp_base
    )
}

#' Collect regionplot result after the background process has finished.
#' @noRd
.collect_regionplot_result <- function(handle) {
    pngs <- Sys.glob(file.path(handle$tmp_base, "*.png"))
    if (length(pngs) == 0)
        return(list(error = "No plot was generated. Check that the region contains SNPs in the p-value files."))
    list(path = pngs[1])
}

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Empty regions skeleton
#' @noRd
.empty_regions <- function() {
    data.table::data.table(
        region_id  = character(),
        chr        = character(),
        start      = integer(),
        end        = integer(),
        length     = integer(),
        snp_count  = integer(),
        snp_ids    = character(),
        traits     = character(),
        methods    = character(),
        min_pvalue = numeric()
    )
}

#' Resolve full GFF path from project config
#' @noRd
.resolve_gff_path <- function(pd) {
    cfg     <- pd$config
    inp_dir <- config_get(cfg, "input", "dir",  default = "data")
    gff_rel <- config_get(cfg, "input", "gff",  default = "")
    if (is.null(gff_rel) || gff_rel == "") return("")
    file.path(get_pipeline_path(), inp_dir, gff_rel)
}

#' Resolve topr-converted GFF path (from intermediate directory)
#' @noRd
.resolve_gff_topr_path <- function(pd) {
    # Pipeline places gff2topr output in _intermediate/annotation/
    Sys.glob(mod_path(pd$name, MOD_INTER, "annotation", "*.tsv"))[1] %||% ""
}

#' Ensure topr-format GFF exists, converting from normalized GFF3 if needed.
#' Returns the path to the topr TSV, or "" on failure.
#' @noRd
.ensure_gff_topr <- function(pd) {
    # Return immediately if already exists (pipeline ran gff2topr rule)
    existing <- .resolve_gff_topr_path(pd)
    if (file_ok(existing)) return(existing)

    # Locate the normalized GFF3 written by the processing mode
    norm_gff <- mod_path(pd$name, MOD_INTER, "annotation", "normalized.gff3")
    if (!file_ok(norm_gff)) return("")

    # Read GFF config values
    cfg       <- pd$config
    feature   <- config_get(cfg, "gff", "feature",   default = "mRNA")
    gene_name <- config_get(cfg, "gff", "gene_name", default = "description")
    biotype   <- config_get(cfg, "gff", "biotype",   default = "biotype")

    out_path <- mod_path(pd$name, MOD_INTER, "annotation", "topr_gene_annotation.tsv")
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

    pipeline_path <- get_pipeline_path()
    script <- file.path(pipeline_path, "scripts", "gff2topr.py")
    if (!file.exists(script)) return("")

    res <- tryCatch(
        processx::run(
            command = "python3",
            args    = c(script, norm_gff, feature, gene_name, biotype, out_path),
            wd      = pipeline_path,
            timeout = 30
        ),
        error = function(e) list(status = 1L)
    )

    if (res$status == 0L && file_ok(out_path)) out_path else ""
}

#' Build colon-separated ASSOC_TABLES arg for plot_regionplot.R
#' Format: "method:adjust_str:pvalues_K{k}.tsv,..."
#' @noRd
.build_assoc_tables_arg <- function(method_files) {
    if (length(method_files) == 0) return(NULL)
    # Each file is like .../methods/{method}/{method}_pvalues_K{k}_sig_snps_{adjust}.tsv
    # We need the pvalues (not sig_snps) file: {method}_pvalues_K{k}.tsv
    entries <- lapply(method_files, function(f) {
        parts <- strsplit(basename(f), "_")[[1]]
        # method_pvalues_K{k}_sig_snps_{adjust}.tsv → extract method, k, adjust
        sig_idx <- which(parts == "sig" & c(FALSE, parts[-length(parts)] == "snps") == FALSE)
        # Simpler: use regex on filename
        m <- regmatches(basename(f),
            regexec("^(.+?)_pvalues_(K[0-9.]+)_sig_snps_(.+)\\.tsv$",
                    basename(f)))[[1]]
        if (length(m) < 4) return(NULL)
        method <- m[2]; k_str <- m[3]; adjust <- m[4]
        pval_file <- file.path(dirname(f),
                               paste0(method, "_pvalues_", k_str, ".tsv"))
        if (!file.exists(pval_file)) return(NULL)
        paste(method, adjust, pval_file, sep = ":")
    })
    entries <- Filter(Negate(is.null), entries)
    if (length(entries) == 0) return(NULL)
    paste(entries, collapse = ",")
}
