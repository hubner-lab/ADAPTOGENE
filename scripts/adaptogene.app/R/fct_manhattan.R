#' Load coordinate mapping JSON from pipeline
#' @param coords_path path to _coords.json file
#' @return list with chr_offsets, x_range, y_range, bonferroni_y etc.
#' @noRd
load_coords <- function(coords_path) {
    if (!file.exists(coords_path)) return(NULL)
    tryCatch(
        jsonlite::read_json(coords_path, simplifyVector = TRUE),
        error = function(e) NULL
    )
}

#' Compute chromosome midpoints from coords for plotly tick labels
#' @noRd
chr_midpoints <- function(coords) {
    offsets <- unlist(coords$chr_offsets)
    lengths <- unlist(coords$chr_lengths)
    offsets + lengths / 2
}

#' Encode a PNG file as a base64 data URI for plotly layout.images
#' @noRd
encode_background_png <- function(bg_path) {
    if (!file.exists(bg_path)) return(NULL)
    base64enc::dataURI(file = bg_path, mime = "image/png")
}

#' Build region rectangle shapes for plotly layout
#' @param regions data.table with chr, start, end, region_id columns
#' @param coords coordinate mapping list
#' @param y_lo y minimum for rectangles
#' @param y_hi y maximum for rectangles
#' @param fill_color plotly fillcolor string (rgba)
#' @param line_color plotly line color string (rgba)
#' @noRd
build_region_shapes <- function(regions, coords, y_lo = NULL, y_hi = NULL,
                                 fill_color = "rgba(100, 149, 237, 0.12)",
                                 line_color = "rgba(100, 149, 237, 0.3)") {
    if (is.null(regions) || nrow(regions) == 0) return(list())
    if (is.null(y_lo)) y_lo <- coords$y_range[1]
    if (is.null(y_hi)) y_hi <- coords$y_range[2]

    # Compute cumulative positions for regions
    chr_offsets <- unlist(coords$chr_offsets)
    regions <- data.table::copy(regions)
    regions[, chr := as.character(chr)]
    regions[, offset := chr_offsets[chr]]
    regions[, start_cum := start + offset]
    regions[, end_cum   := end   + offset]

    # Ensure minimum visible width (0.3% of x span for clickability)
    x_span <- diff(coords$x_range)
    min_hw  <- x_span * 0.003
    regions[, center_cum   := (start_cum + end_cum) / 2]
    regions[, half_width   := pmax((end_cum - start_cum) / 2, min_hw)]
    regions[, start_cum    := center_cum - half_width]
    regions[, end_cum      := center_cum + half_width]

    # Build plotly shapes list
    lapply(seq_len(nrow(regions)), function(i) {
        list(
            type      = "rect",
            xref      = "x",
            yref      = "y",
            x0        = regions$start_cum[i],
            x1        = regions$end_cum[i],
            y0        = y_lo,
            y1        = y_hi,
            fillcolor = fill_color,
            line      = list(color = line_color, width = 1),
            customdata = regions$region_id[i]
        )
    })
}

#' Build dual-color region shapes: current region blue, others gray
#' @param all_regions data.table of all computed regions
#' @param current_region_id character: region_id of the selected region (or NULL)
#' @param coords coordinate mapping list
#' @param y_lo y minimum; y_hi y maximum
#' @noRd
build_dual_region_shapes <- function(all_regions, current_region_id,
                                      coords, y_lo = NULL, y_hi = NULL) {
    if (is.null(all_regions) || nrow(all_regions) == 0) return(list())
    if (is.null(y_lo)) y_lo <- coords$y_range[1]
    if (is.null(y_hi)) y_hi <- coords$y_range[2]

    if (!is.null(current_region_id) && nzchar(current_region_id)) {
        others  <- all_regions[region_id != current_region_id]
        current <- all_regions[region_id == current_region_id]
    } else {
        others  <- all_regions
        current <- all_regions[0]
    }

    shapes_gray <- build_region_shapes(others, coords, y_lo, y_hi,
        fill_color = "rgba(180, 180, 180, 0.10)",
        line_color = "rgba(180, 180, 180, 0.25)")
    shapes_blue <- build_region_shapes(current, coords, y_lo, y_hi)

    c(shapes_gray, shapes_blue)
}

#' Build the plotly Manhattan overlay plot
#' @param bg_uri base64 encoded background PNG (from encode_background_png)
#' @param coords coordinate mapping list (from load_coords)
#' @param sig_snps data.table of significant SNPs with cum_pos, log10p, trait, method, SNPID, region_id
#' @param regions data.table of regions (optional, for shapes)
#' @param current_region_id character: region_id of the selected region; when set,
#'   that region is drawn blue and all others gray (dual-color mode)
#' @param trait_colors named character vector of trait colors
#' @param method_shapes named integer vector of plotly symbol codes per method
#' @param is_miami logical, TRUE for Miami plots (y axis has negative values)
#' @noRd
build_manhattan_plotly <- function(bg_uri, coords, sig_snps = NULL,
                                    regions = NULL,
                                    current_region_id = NULL,
                                    trait_colors = NULL,
                                    method_shapes = NULL,
                                    is_miami = FALSE,
                                    source = "manhattan_overlay") {
    if (is.null(coords)) return(plotly::event_register(plotly::plot_ly(source = source), "plotly_click"))

    x_range <- coords$x_range
    y_range <- coords$y_range

    # Chromosome tick positions and labels
    chr_mids  <- chr_midpoints(coords)
    chr_names <- names(unlist(coords$chr_offsets))

    # Region shapes (dual-color if a current region is selected)
    shapes <- if (!is.null(regions) && nrow(regions) > 0) {
        if (!is.null(current_region_id) && nzchar(current_region_id %||% "")) {
            build_dual_region_shapes(regions, current_region_id, coords,
                                     y_lo = y_range[1], y_hi = y_range[2])
        } else {
            build_region_shapes(regions, coords,
                                y_lo = y_range[1], y_hi = y_range[2])
        }
    } else {
        list()
    }

    # Threshold lines
    threshold_lines <- list()
    if (!is.null(coords$bonferroni_y)) {
        threshold_lines <- c(threshold_lines, list(
            list(type = "line", xref = "paper", yref = "y",
                 x0 = 0, x1 = 1, y0 = coords$bonferroni_y, y1 = coords$bonferroni_y,
                 line = list(color = "red", dash = "dash", width = 1))
        ))
    }
    if (is_miami) {
        if (!is.null(coords$gwas_threshold_y)) {
            threshold_lines <- c(threshold_lines, list(
                list(type = "line", xref = "paper", yref = "y",
                     x0 = 0, x1 = 1,
                     y0 = coords$gwas_threshold_y, y1 = coords$gwas_threshold_y,
                     line = list(color = "red", dash = "dash", width = 1))
            ))
        }
        # Zero line for Miami
        threshold_lines <- c(threshold_lines, list(
            list(type = "line", xref = "paper", yref = "y",
                 x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                 line = list(color = "black", width = 1.5))
        ))
    }

    p <- plotly::plot_ly(source = source) |>
        plotly::layout(
            paper_bgcolor = "white",
            plot_bgcolor  = "white",
            images = if (!is.null(bg_uri)) {
                list(list(
                    source  = bg_uri,
                    xref    = "x", yref = "y",
                    x       = x_range[1],
                    y       = y_range[2],
                    sizex   = x_range[2] - x_range[1],
                    sizey   = y_range[2] - y_range[1],
                    sizing  = "stretch",
                    layer   = "below"
                ))
            } else {
                list()
            },
            xaxis = list(
                range      = x_range,
                tickmode   = "array",
                tickvals   = unname(chr_mids),
                ticktext   = unname(chr_names),
                tickangle  = -60,
                tickfont   = list(size = 10),
                title      = "",
                showgrid   = FALSE,
                zeroline   = FALSE,
                showline   = FALSE
            ),
            yaxis = if (is_miami) {
                # Show absolute values on both halves of Miami plot
                tick_vals <- pretty(y_range, n = 6)
                tick_vals <- tick_vals[tick_vals != 0]
                list(
                    range      = y_range,
                    title      = "-log\u2081\u2080(p)",
                    tickfont   = list(size = 10),
                    showgrid   = FALSE,
                    zeroline   = FALSE,
                    tickmode   = "array",
                    tickvals   = unname(tick_vals),
                    ticktext   = unname(as.character(abs(tick_vals))),
                    tickformat = ""
                )
            } else {
                list(
                    range      = y_range,
                    title      = "-log\u2081\u2080(p)",
                    tickfont   = list(size = 10),
                    showgrid   = FALSE,
                    zeroline   = FALSE,
                    tickformat = ""
                )
            },
            shapes  = c(shapes, threshold_lines),
            margin  = list(l = 60, r = 20, t = 30, b = 60),
            legend  = list(orientation = "v", x = 1, xanchor = "left", y = 1),
            hovermode = "closest"
        )

    # Add sig SNP markers
    if (!is.null(sig_snps) && nrow(sig_snps) > 0) {
        # Map trait colors to plotly marker colors
        if (is.null(trait_colors)) {
            okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                           "#0072B2", "#D55E00", "#CC79A7", "#999999")
            traits <- unique(sig_snps$trait)
            trait_colors <- setNames(okabe_ito[seq_along(traits)], traits)
        }

        # Default method symbols (plotly symbol names)
        if (is.null(method_shapes)) {
            sym_names <- c("circle", "triangle-up", "square", "diamond",
                           "cross", "x", "star", "triangle-down")
            methods <- unique(sig_snps$method)
            method_shapes <- setNames(sym_names[seq_along(methods)], methods)
        }

        sig_snps <- data.table::copy(sig_snps)
        sig_snps[, marker_color  := trait_colors[trait]]
        sig_snps[, marker_symbol := method_shapes[method]]

        # Build hover text
        sig_snps[, hover_text := paste0(
            "SNP: ", SNPID, "<br>",
            "Trait: ", trait, "<br>",
            "Method: ", method, "<br>",
            "Region: ", region_id, "<br>",
            "-log10(p): ", round(abs(log10p), 3)
        )]

        # Data traces — one per (trait, method) combo, no legend entries
        combos <- unique(sig_snps[, .(method, trait)])
        for (ci in seq_len(nrow(combos))) {
            meth   <- combos$method[ci]
            tr     <- combos$trait[ci]
            df_sub <- sig_snps[method == meth & trait == tr]
            p <- p |>
                plotly::add_markers(
                    data       = df_sub,
                    x          = ~cum_pos,
                    y          = ~log10p,
                    marker     = list(
                        color   = trait_colors[tr],
                        symbol  = method_shapes[meth],
                        size    = 10,
                        opacity = 0.9,
                        line    = list(width = 1.2, color = "black")
                    ),
                    text       = ~hover_text,
                    hoverinfo  = "text",
                    customdata = ~region_id,
                    showlegend = FALSE
                )
        }

        # Legend block 1: Trait → Color
        # Place dummy point far outside visible range so it's clipped but trace is valid
        x_leg <- x_range[1] - abs(diff(x_range)) * 10
        y_leg <- y_range[2] + abs(diff(y_range)) * 10
        for (tr in unique(sig_snps$trait)) {
            p <- p |>
                plotly::add_markers(
                    x                = x_leg,
                    y                = y_leg,
                    marker           = list(color = trait_colors[tr],
                                           symbol = "circle", size = 10),
                    name             = tr,
                    legendgroup      = "trait",
                    legendgrouptitle = list(text = "Trait"),
                    showlegend       = TRUE,
                    hoverinfo        = "skip"
                )
        }

        # Legend block 2: Method → Shape (neutral gray)
        for (meth in unique(sig_snps$method)) {
            p <- p |>
                plotly::add_markers(
                    x                = x_leg,
                    y                = y_leg,
                    marker           = list(color = "rgba(80,80,80,0.85)",
                                           symbol = method_shapes[meth], size = 10,
                                           line = list(width = 1, color = "black")),
                    name             = meth,
                    legendgroup      = "method",
                    legendgrouptitle = list(text = "Method"),
                    showlegend       = TRUE,
                    hoverinfo        = "skip"
                )
        }
    }

    plotly::event_register(p, "plotly_click")
}

#' Compute cumulative positions for sig SNPs using coords chr_offsets
#' @noRd
add_cum_pos <- function(sig_snps, coords) {
    chr_offsets <- unlist(coords$chr_offsets)
    sig_snps <- data.table::copy(sig_snps)
    sig_snps[, chr := as.character(chr)]
    sig_snps[, cum_pos := pos + chr_offsets[chr]]
    sig_snps[, log10p  := -log10(pvalue)]
    sig_snps
}
