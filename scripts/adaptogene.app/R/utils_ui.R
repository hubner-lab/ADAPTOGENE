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
#' @param source optional source label (for overlap tab)
#' @noRd
region_info_bar <- function(region_id, n_snps = NULL, n_exon = NULL,
                              n_promoter = NULL, source = NULL) {
    disp <- format_region_id(region_id)
    parts <- character(0)
    if (!is.null(n_snps))    parts <- c(parts, paste0(n_snps, " SNPs"))
    if (!is.null(n_exon) && n_exon > 0)
        parts <- c(parts, paste0(n_exon, " exonic"))
    if (!is.null(n_promoter) && n_promoter > 0)
        parts <- c(parts, paste0(n_promoter, " promoter"))
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
