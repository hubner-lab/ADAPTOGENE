# On-the-fly combine strategies for interactive Shiny exploration.
# Ports combine logic from scripts/combine_selected_snps.R using data.table
# (avoids GenomicRanges dependency in Shiny).

#' Okabe-Ito colour palette (8 colours, colour-blind safe)
#' @noRd
OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' Map traits to Okabe-Ito colours (sorted order, consistent with build_manhattan_plotly)
#'
#' @param traits character vector of trait names
#' @return named character vector: trait → hex colour
#' @noRd
trait_color_map <- function(traits) {
    traits_sorted <- sort(unique(traits))
    n <- length(traits_sorted)
    colors <- rep(OKABE_ITO, ceiling(n / length(OKABE_ITO)))[seq_len(n)]
    stats::setNames(colors, traits_sorted)
}

#' Combine per-method sig SNPs using a specified strategy (on-the-fly in Shiny)
#'
#' @param sigsnps_list Named list of data.tables, one per method.
#'   Each must have columns: SNPID, chr, pos, pvalue, method, trait.
#' @param strategy One of "Sum", "Overlap", "PairOverlap", or a method name.
#' @param gap Integer. Max bp distance for spatial overlap (Overlap/PairOverlap).
#' @return Long-format data.table with SNPID, chr, pos, pvalue, method, trait.
#'   Empty data.table if nothing passes.
#' @noRd
combine_sigsnps <- function(sigsnps_list, strategy = "Sum", gap = 200000L) {
    if (length(sigsnps_list) == 0) return(.empty_sigsnps())

    # Drop empty methods
    sigsnps_list <- sigsnps_list[vapply(sigsnps_list, nrow, integer(1)) > 0]
    if (length(sigsnps_list) == 0) return(.empty_sigsnps())

    methods_vec <- names(sigsnps_list)
    gap <- as.integer(gap)

    result <- if (strategy %in% methods_vec) {
        # ── Single-method ───────────────────────────────────────────────────────
        data.table::copy(sigsnps_list[[strategy]])

    } else if (strategy == "Sum") {
        # ── Union: all sig SNPs from all methods ────────────────────────────────
        data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)

    } else if (strategy == "Overlap") {
        # ── Cross-method spatial overlap (trait-agnostic) ───────────────────────
        # A SNP passes if any partner exists in another method within gap bp.
        selected_ids <- character(0)
        for (n1 in methods_vec) {
            for (n2 in methods_vec) {
                if (n1 == n2) next
                ids <- .spatial_overlap_ids(sigsnps_list[[n1]], sigsnps_list[[n2]], gap)
                selected_ids <- c(selected_ids, ids)
            }
        }
        selected_ids <- unique(selected_ids)
        if (length(selected_ids) == 0) return(.empty_sigsnps())
        all_snps <- data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)
        all_snps[SNPID %in% selected_ids]

    } else if (strategy == "PairOverlap") {
        # ── Cross-method spatial overlap, same-trait only ───────────────────────
        all_snps <- data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)
        all_traits <- unique(all_snps$trait)
        selected_rows <- list()
        for (n1 in methods_vec) {
            for (n2 in methods_vec) {
                if (n1 == n2) next
                for (bio in all_traits) {
                    m1_bio <- sigsnps_list[[n1]][trait == bio]
                    m2_bio <- sigsnps_list[[n2]][trait == bio]
                    if (nrow(m1_bio) == 0 || nrow(m2_bio) == 0) next
                    ids1 <- .spatial_overlap_ids(m1_bio, m2_bio, gap)
                    ids2 <- .spatial_overlap_ids(m2_bio, m1_bio, gap)
                    if (length(ids1) > 0)
                        selected_rows <- c(selected_rows, list(m1_bio[SNPID %in% ids1]))
                    if (length(ids2) > 0)
                        selected_rows <- c(selected_rows, list(m2_bio[SNPID %in% ids2]))
                }
            }
        }
        if (length(selected_rows) == 0) return(.empty_sigsnps())
        unique(data.table::rbindlist(selected_rows, use.names = TRUE, fill = TRUE))

    } else {
        stop(paste0("Unknown combine strategy: ", strategy))
    }

    if (is.null(result) || nrow(result) == 0) return(.empty_sigsnps())

    # Deduplicate — same SNP can appear in multiple methods/traits
    result <- unique(result, by = c("SNPID", "method", "trait"))

    # For combined view, y-position uses min pvalue across methods per SNP
    # (so dots align with the combined background PNG's y-calibration)
    result[, pvalue := min(pvalue, na.rm = TRUE), by = "SNPID"]

    result
}

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Find SNP IDs from dt1 that have a partner in dt2 within gap bp on same chr
#' @noRd
.spatial_overlap_ids <- function(dt1, dt2, gap) {
    if (nrow(dt1) == 0 || nrow(dt2) == 0) return(character(0))

    # Expand dt1 to windows [pos-gap, pos+gap]; dt2 stays as points [pos, pos]
    win <- data.table::data.table(
        chr   = as.character(dt1$chr),
        start = pmax(1L, as.integer(dt1$pos) - gap),
        end   = as.integer(dt1$pos) + gap,
        SNPID = dt1$SNPID
    )
    pts <- data.table::data.table(
        chr   = as.character(dt2$chr),
        start = as.integer(dt2$pos),
        end   = as.integer(dt2$pos),
        SNPID = dt2$SNPID
    )
    data.table::setkey(win, chr, start, end)
    data.table::setkey(pts, chr, start, end)
    ov <- data.table::foverlaps(pts, win, type = "within", nomatch = NULL)
    if (nrow(ov) == 0) return(character(0))
    # Return IDs from dt1 (window side)
    unique(ov$SNPID)
}

#' Empty sig SNPs skeleton
#' @noRd
.empty_sigsnps <- function() {
    data.table::data.table(
        SNPID  = character(),
        chr    = character(),
        pos    = integer(),
        pvalue = numeric(),
        method = character(),
        trait  = character()
    )
}
