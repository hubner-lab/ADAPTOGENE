#' Create a "plot not available" placeholder card body
#' @param message primary message text
#' @param suggestion optional smaller suggestion text (e.g. which pipeline mode to run)
#' @noRd
plot_placeholder <- function(message = "Plot not available", suggestion = NULL) {
    htmltools::div(
        class = "plot-placeholder",
        bsicons::bs_icon("image", size = "3em"),
        htmltools::p(message),
        if (!is.null(suggestion))
            htmltools::p(class = "text-muted small mt-1", suggestion)
    )
}

#' Region info bar widget
#' @param region_id region ID string
#' @param n_snps number of SNPs in region
#' @param n_exon exonic SNP count
#' @param n_promoter promoter SNP count
#' @param n_genes number of genes overlapping region
#' @param source optional source label (for overlap tab)
#' @noRd
region_info_bar <- function(region_id, n_snps = NULL, n_exon = NULL,
                              n_promoter = NULL, n_genes = NULL, source = NULL) {
    disp <- format_region_id(region_id)
    parts <- character(0)
    if (!is.null(n_snps))    parts <- c(parts, paste0(n_snps, " sig. SNPs"))
    if (!is.null(n_exon) && n_exon > 0)
        parts <- c(parts, paste0(n_exon, " exonic"))
    if (!is.null(n_promoter) && n_promoter > 0)
        parts <- c(parts, paste0(n_promoter, " promoter"))
    if (!is.null(n_genes))   parts <- c(parts, paste0(n_genes, " gene", if (n_genes != 1) "s"))
    if (!is.null(source))    parts <- c(parts, paste0("source: ", source))

    htmltools::div(
        class = "region-info-bar",
        htmltools::strong(disp),
        if (length(parts) > 0)
            htmltools::span(paste0(" \u2014 ", paste(parts, collapse = " | ")))
    )
}

#' Mode status indicator row
#' @param name display name
#' @param ok logical
#' @noRd
mode_status_row <- function(name, ok) {
    icon_cls <- if (ok) "status-ok"      else "status-missing"
    icon     <- if (ok) "\u2713"         else "\u2014"
    htmltools::div(
        class = "mode-status",
        htmltools::span(class = paste("status-icon", icon_cls), icon),
        htmltools::span(name)
    )
}

#' Small badge showing a config parameter value
#' @param label short label (e.g. "K", "combine")
#' @param value the value to display
#' @param class Bootstrap badge background class (default "bg-secondary")
#' @noRd
config_badge <- function(label, value, class = "bg-secondary") {
    htmltools::tags$span(
        class = paste("badge me-1", class),
        style = "font-size:0.75rem; font-weight:500;",
        paste0(label, ": ", value)
    )
}

#' Container bar for config parameter badges
#' @param ... config_badge() calls
#' @noRd
config_badges_bar <- function(...) {
    htmltools::div(
        class = "d-flex flex-wrap align-items-center gap-1 mb-2",
        style = "opacity:0.85;",
        bsicons::bs_icon("gear-fill", size = "0.8em", class = "text-muted me-1"),
        ...
    )
}

#' Info-circle popover explaining WZA window size and how to change it
#'
#' @param context One of "gea", "gwas", or "overlap"
#' @noRd
wza_window_note <- function(context = "gea") {
    body <- switch(context,
        gea = htmltools::tagList(
            htmltools::p(
                "WZA windows are ", htmltools::strong("precomputed by the pipeline"),
                " — they are not recalculated live."
            ),
            htmltools::p(
                "To change the WZA window size: open ",
                htmltools::strong("Config → GEA → Advanced → ‘WZA window size’"),
                " and set a fixed bp value, or ",
                htmltools::code("auto_genome_wide"),
                " / ",
                htmltools::code("auto_per_chromosome"),
                " (LD-derived). Then re-run the ",
                htmltools::strong("GEA"),
                " module; the new windows appear here after the run finishes."
            )
        ),
        gwas = htmltools::tagList(
            htmltools::p(
                "WZA windows are ", htmltools::strong("precomputed by the pipeline"),
                " — they are not recalculated live."
            ),
            htmltools::p(
                "To change the WZA window size: open ",
                htmltools::strong("Config → GWAS → Advanced → ‘WZA window size’"),
                " and set a fixed bp value, or ",
                htmltools::code("auto_genome_wide"),
                " / ",
                htmltools::code("auto_per_chromosome"),
                " (LD-derived). Then re-run the ",
                htmltools::strong("GWAS"),
                " module; the new windows appear here after the run finishes."
            )
        ),
        overlap = htmltools::tagList(
            htmltools::p(
                "The WZA windows in this Miami plot are ",
                htmltools::strong("inherited from the GEA and GWAS modules"),
                " — there is no separate WZA window setting for GEAxGWAS."
            ),
            htmltools::p(
                "To use a different WZA distance: change ",
                htmltools::strong("‘WZA window size’"),
                " in the GEA and/or GWAS Advanced config, re-run those modules, then reload —",
                " the Miami plot updates automatically."
            )
        )
    )

    bslib::popover(
        trigger = bsicons::bs_icon("info-circle",
                                   title = "WZA window size",
                                   class = "text-muted ms-1",
                                   style = "cursor:pointer;"),
        title = "WZA window size",
        body
    )
}

#' Hover note announcing the SNP -> WZA-window collapse under the WZA regime.
#'
#' Appears when the user switches the regime toggle ON; disappears when OFF.
#' stats is the list returned by wza_collapse_stats(); pass NULL to suppress.
#' Converted from a full-width alert banner to a compact hover badge (same
#' wording, now in the tooltip) — ambient FYI, not a page-wide banner.
#' @noRd
wza_collapse_note <- function(stats) {
    if (is.null(stats)) return(NULL)
    snp_txt <- if (!is.na(stats$n_snps))
        paste0(format(stats$n_snps, big.mark = ","), " SNPs collapsed into ") else ""
    win_txt <- if (!is.na(stats$window_bp))
        paste0(" (≈ ", format(stats$window_bp, big.mark = ","), " bp/window)") else ""
    htmltools::div(
        class = "d-flex justify-content-end mb-2",
        filter_note(
            paste0(format(stats$n_windows, big.mark = ","), " windows"),
            htmltools::p(
                htmltools::strong("WZA regime active"), " — ", snp_txt,
                htmltools::strong(format(stats$n_windows, big.mark = ","), " windows"),
                win_txt, "."
            ),
            class = "bg-secondary"
        )
    )
}

#' Small badge that reveals a note in a hover tooltip.
#'
#' Generic building block for the "instructional text belongs in Shiny, not baked into
#' plots" pattern: a compact badge (icon + short label) sits next to a plot's card
#' title; hovering (bslib::tooltip, not a click-popover) reveals `body` — one short
#' sentence or a few, kept out of the plot image itself.
#'
#' @param label short text shown on the badge itself (e.g. a threshold value, a count)
#' @param body tooltip content — a string, or htmltools tags/tagList for multi-paragraph notes
#' @param class Bootstrap badge background class (default neutral "bg-secondary";
#'   use "bg-danger"/"bg-success" only for a genuine status signal, not a fixed threshold)
#' @param placement tooltip placement, passed to bslib::tooltip()
#' @noRd
filter_note <- function(label, body, class = "bg-secondary", placement = "right") {
    trigger <- htmltools::tags$span(
        class = paste("badge", class),
        style = "cursor:default; font-size:0.7rem; font-weight:500;",
        bsicons::bs_icon("info-circle-fill", size = "0.75em"), " ", label
    )
    bslib::tooltip(trigger, body, placement = placement)
}

#' Colored hover note for the relatedness (IBD pi-hat) histogram.
#'
#' Replaces the old static caption alert with a compact green/red badge next to the
#' plot title. Hovering (bslib::tooltip, not a click-popover) reveals the high-selfing
#' caveat, the pair count, and action-specific guidance (nothing baked into the plot —
#' data-derived numbers only, per the "instructional text belongs in Shiny" rule).
#'
#' Unlike filter_note()'s default neutral badge, this one is green/red because pair
#' count is a genuine status signal ("related pairs found — look"), not a fixed threshold.
#'
#' @param threshold Filter.relatedness pi-hat threshold, or NULL/NA when unset
#' @param action Filter.relatedness_action ("keep" or "remove")
#' @param n_pairs_above number of pairs with pi-hat above threshold (NA if unknown)
#' @param n_would_remove number of samples the greedy filter would drop (NA if unknown)
#' @noRd
relatedness_note <- function(threshold, action, n_pairs_above, n_would_remove) {
    if (is.null(threshold) || is.na(threshold)) return(NULL)
    action <- tolower(action %||% "keep")
    is_red <- !is.na(n_pairs_above) && n_pairs_above > 0
    badge_class <- if (is_red) "bg-danger" else "bg-success"
    badge_label <- if (is.na(n_pairs_above)) "—"
                   else paste0(n_pairs_above, " pair", if (n_pairs_above != 1) "s" else "")

    action_txt <- if (action == "remove") {
        if (is.na(n_would_remove)) "Related samples removed."
        else paste0(n_would_remove, " related sample", if (n_would_remove != 1) "s" else "", " removed.")
    } else {
        paste0(
            "Nothing removed (action = keep). ",
            if (!is.na(n_would_remove))
                paste0("Would remove ", n_would_remove, " sample",
                       if (n_would_remove != 1) "s" else "", " if switched to remove. ")
            else "",
            "Set ", "Filter.relatedness_action: remove",
            " and re-run Processing to drop the higher-missingness member of each pair."
        )
    }

    body <- htmltools::tagList(
        htmltools::p(
            htmltools::strong("High-selfing species"),
            " show naturally inflated pi-hat — pick a threshold from the histogram, ",
            "not the standard outbred 0.2 default."
        ),
        htmltools::p(
            if (is.na(n_pairs_above)) "No pair data available."
            else paste0(n_pairs_above, " pair", if (n_pairs_above != 1) "s" else "",
                       " above pi-hat > ", threshold, ".")
        ),
        htmltools::p(action_txt)
    )

    filter_note(badge_label, body, class = badge_class)
}

#' A card header with title + download popover
#' @noRd
card_header_with_download <- function(ns, title, dl_id_svg = NULL, dl_id_png = NULL) {
    dl_btn <- if (!is.null(dl_id_svg) || !is.null(dl_id_png)) {
        bslib::popover(
            trigger = bsicons::bs_icon("download", title = "Download"),
            title = "Download",
            if (!is.null(dl_id_svg))
                shiny::downloadButton(ns(dl_id_svg), "SVG",
                                      class = "btn-sm btn-outline-secondary"),
            if (!is.null(dl_id_png))
                shiny::downloadButton(ns(dl_id_png), "PNG",
                                      class = "btn-sm btn-outline-secondary")
        )
    } else {
        NULL
    }

    bslib::card_header(
        class = "d-flex justify-content-between align-items-center",
        shiny::textOutput(ns(title), inline = TRUE),
        if (!is.null(dl_btn)) htmltools::span(class = "d-flex gap-2", dl_btn)
    )
}
