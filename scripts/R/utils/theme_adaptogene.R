# theme_adaptogene.R — shared ggplot2 theme + palette matched to the Shiny
# app's "Precision Genomics" bslib theme (scripts/adaptogene.app/R/app_theme.R).
# Usage: source("/pipeline/scripts/R/utils/theme_adaptogene.R")
# Requires: ggplot2 loaded by the sourcing script.
#
# Rolled out to the Processing module first (plot_qc_processing.R). Intent is
# to extend this to the rest of the pipeline's plot scripts over time — see
# CLAUDE.md "Code Structure Rules". Colors only for now: no base_family is set
# (Docker installs no fonts; setting one would silently fall back to default
# sans while implying a match that isn't there — same trap as the existing
# dead base_family='Helvetica' settings elsewhere).

# Core palette — mirrors app_theme.R's bs_theme() colors exactly.
ADAPT_COL <- list(
    primary   = "#1B7A6E",  # deep teal-green
    secondary = "#4A5568",  # slate
    success   = "#38A169",
    info      = "#3182CE",
    warning   = "#D69E2E",  # amber
    danger    = "#E53E3E",  # red
    fg        = "#2D3748",
    muted     = "#718096"
)

# Semantic aliases — use these in plot code, not raw hex, so the meaning of a
# color stays consistent across every processing plot (and future ones).
# Data-point/fill colors below are the last 3 of wesanderson::wes_palette
# ("Rushmore1") = c("#E1BD6D","#EABE94","#0B775E","#35274A","#F2300F") — user
# preference, chosen over the app-matched teal/slate/amber. Note the middle
# value is a dark plum/purple, not literally blue.
ADAPT_RETAINED  <- "#0B775E"  # green — sample/SNP passed a filter
ADAPT_REMOVED   <- "#F2300F"  # red   — sample/SNP dropped by a filter
ADAPT_THRESHOLD <- "#35274A"  # dark plum/purple — threshold/cutoff lines
ADAPT_NEUTRAL   <- ADAPT_RETAINED  # plain/no-flag data points — same green as
                                   # retained by request (no more grey dots)

# Okabe-Ito colorblind-safe categorical palette (also used by
# scripts/R/utils/manhattan_utils.R — kept identical here to converge the
# pipeline's two independent categorical palettes onto one going forward).
ADAPT_CATEGORICAL <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                        "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' Shared ggplot2 theme — "Publication Classic"
#'
#' theme_classic() base (axis lines, no gridlines) + Precision Genomics sizing.
#' No base_family: device-default sans until the Docker image ships fonts.
#' @noRd
theme_adaptogene <- function(base_size = 11) {
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
        plot.title      = ggplot2::element_text(size = base_size, face = "bold", color = ADAPT_COL$fg),
        plot.subtitle   = ggplot2::element_text(size = base_size - 2, color = ADAPT_COL$muted),
        axis.title      = ggplot2::element_text(size = base_size - 1, color = ADAPT_COL$fg),
        axis.text       = ggplot2::element_text(color = ADAPT_COL$secondary),
        axis.line       = ggplot2::element_line(color = ADAPT_COL$secondary, linewidth = 0.4),
        axis.ticks      = ggplot2::element_line(color = ADAPT_COL$secondary, linewidth = 0.4),
        legend.position = "bottom",
        legend.title    = ggplot2::element_text(size = base_size - 1, color = ADAPT_COL$fg),
        legend.text     = ggplot2::element_text(size = base_size - 2, color = ADAPT_COL$secondary)
    )
}

#' Categorical color scale using the Okabe-Ito palette above 8 levels of grouping.
#' @noRd
scale_color_adaptogene <- function(...) {
    ggplot2::scale_color_manual(values = ADAPT_CATEGORICAL, ...)
}

#' @noRd
scale_fill_adaptogene <- function(...) {
    ggplot2::scale_fill_manual(values = ADAPT_CATEGORICAL, ...)
}
