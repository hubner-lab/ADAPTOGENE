#' Manhattan overlay module UI
#'
#' Displays a plotly Manhattan with a static background PNG (non-sig SNPs)
#' and interactive sig-SNP markers overlaid.
#'
#' @param id module namespace id
#' @param height plotly output height (default "500px")
#' @noRd
mod_manhattan_overlay_ui <- function(id, height = "500px") {
    ns <- shiny::NS(id)
    bslib::card(
        full_screen = TRUE,
        bslib::card_header(
            class = "d-flex justify-content-between align-items-center",
            shiny::textOutput(ns("title"), inline = TRUE),
            htmltools::span(
                class = "d-flex gap-2 align-items-center",
                shiny::uiOutput(ns("show_regions_toggle")),
                bslib::popover(
                    trigger = bsicons::bs_icon("download", title = "Download plot"),
                    title = "Download",
                    shiny::downloadButton(ns("dl_svg"), "SVG",
                                          class = "btn-sm btn-outline-secondary"),
                    shiny::downloadButton(ns("dl_png"), "PNG",
                                          class = "btn-sm btn-outline-secondary")
                )
            )
        ),
        bslib::card_body(
            class = "p-1",
            plotly::plotlyOutput(ns("overlay"), height = height)
        )
    )
}

#' Manhattan overlay module server
#'
#' @param id module namespace id
#' @param project_data reactive returning project data bundle from make_project_data()
#' @param module character: MOD_ASSOC, MOD_PHENO, or MOD_OVERLAP
#' @param method character reactive: method name (or NULL for combined)
#' @param trait character reactive: trait name (or NULL for combined)
#' @param combined logical: TRUE = combined multi-trait/method view
#' @param is_miami logical: TRUE = Miami plot mode
#' @param title_label reactive character for card title
#' @param regions reactive data.table of regions (optional, for overlay shapes)
#' @param show_regions_control logical: show toggle for region highlights
#' @return reactive of clicked region_id (character or NULL)
#' @noRd
mod_manhattan_overlay_server <- function(id, project_data,
                                          module        = MOD_ASSOC,
                                          method        = shiny::reactive(NULL),
                                          trait         = shiny::reactive(NULL),
                                          combined      = FALSE,
                                          is_miami      = FALSE,
                                          title_label   = shiny::reactive("Manhattan"),
                                          regions       = shiny::reactive(NULL),
                                          show_regions_control = TRUE) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        output$title <- shiny::renderText(title_label())

        # Toggle for showing region rectangles (optional)
        output$show_regions_toggle <- shiny::renderUI({
            if (!show_regions_control) return(NULL)
            bslib::input_switch(ns("show_regions"), "Show regions", value = TRUE)
        })

        show_regions <- shiny::reactive({
            if (!show_regions_control) return(FALSE)
            isTRUE(input$show_regions)
        })

        # ── Resolve file paths ─────────────────────────────────────────────────
        bg_path <- shiny::reactive({
            pd  <- project_data()
            k   <- pd$k_best
            if (is.na(k)) return(NULL)

            if (is_miami) {
                if (show_regions() && file_ok(miami_regions_bg_path(pd$name, k)))
                    miami_regions_bg_path(pd$name, k)
                else
                    miami_bg_path(pd$name, k)
            } else if (combined) {
                if (show_regions() &&
                    file_ok(combined_manhattan_regions_bg_path(pd$name, module, k)))
                    combined_manhattan_regions_bg_path(pd$name, module, k)
                else
                    combined_manhattan_bg_path(pd$name, module, k)
            } else {
                m <- method(); t <- trait()
                shiny::req(m, t)
                cfg <- pd$config
                adj <- resolve_adjust(cfg, m, if (module == MOD_PHENO) "phenotype_association" else "association")
                if (is.null(adj)) return(NULL)
                if (show_regions() &&
                    file_ok(manhattan_regions_bg_path(pd$name, module, m, t, k, adj)))
                    manhattan_regions_bg_path(pd$name, module, m, t, k, adj)
                else
                    manhattan_bg_path(pd$name, module, m, t, k, adj)
            }
        })

        coords_path <- shiny::reactive({
            pd <- project_data()
            k  <- pd$k_best
            if (is.na(k)) return(NULL)

            if (is_miami)       miami_coords_path(pd$name, k)
            else if (combined)  combined_manhattan_coords_path(pd$name, module, k)
            else {
                m <- method(); t <- trait()
                shiny::req(m, t)
                adj <- resolve_adjust(pd$config, m,
                    if (module == MOD_PHENO) "phenotype_association" else "association")
                if (is.null(adj)) return(NULL)
                manhattan_coords_path(pd$name, module, m, t, k, adj)
            }
        })

        sig_snps_path <- shiny::reactive({
            pd <- project_data()
            selected_snps_path(pd$name, module)
        })

        # ── Load data ──────────────────────────────────────────────────────────
        coords <- shiny::reactive({
            load_cached(paste0("coords_", coords_path()), function() {
                load_coords(coords_path())
            })
        })

        sig_snps <- shiny::reactive({
            path <- sig_snps_path()
            load_cached(paste0("sig_snps_", path), function() {
                load_selected_snps(project_data()$name, module)
            })
        }) |> shiny::bindCache(project_data()$name, module)

        # ── Render plotly ──────────────────────────────────────────────────────
        output$overlay <- plotly::renderPlotly({
            co <- shiny::req(coords())
            bg <- bg_path()

            bg_uri <- if (file_ok(bg)) encode_background_png(bg) else NULL

            snps <- sig_snps()
            # Filter to current method/trait if not combined
            if (!combined && !is_miami && nrow(snps) > 0) {
                m <- method(); t <- trait()
                if (!is.null(m) && "method" %in% names(snps))
                    snps <- snps[snps$method == m, ]
                if (!is.null(t) && "trait" %in% names(snps))
                    snps <- snps[snps$trait == t, ]
            }

            # Add cumulative positions for sig SNPs
            if (nrow(snps) > 0 && !"cum_pos" %in% names(snps)) {
                snps <- add_cum_pos(snps, co)
            }

            # Regions overlay data
            reg_data <- if (show_regions()) regions() else NULL

            build_manhattan_plotly(
                bg_uri      = bg_uri,
                coords      = co,
                sig_snps    = if (nrow(snps) > 0) snps else NULL,
                regions     = reg_data,
                is_miami    = is_miami,
                source      = ns("overlay")
            )
        })

        # ── Sig SNP click → region_id reactive ────────────────────────────────
        selected_region <- shiny::reactive({
            click <- plotly::event_data("plotly_click", source = ns("overlay"))
            if (is.null(click) || is.null(click$customdata)) return(NULL)
            click$customdata[[1]]
        })

        # ── Download handlers (static files) ──────────────────────────────────
        output$dl_png <- shiny::downloadHandler(
            filename = function() {
                paste0("manhattan_", project_data()$name, "_K",
                       project_data()$k_best, ".png")
            },
            content = function(file) {
                p <- bg_path()
                if (file_ok(p)) file.copy(p, file) else writeLines("", file)
            }
        )

        output$dl_svg <- shiny::downloadHandler(
            filename = function() {
                paste0("manhattan_", project_data()$name, "_K",
                       project_data()$k_best, ".svg")
            },
            content = function(file) {
                p <- bg_path()
                svg_path <- sub("_background\\.png$", ".svg", p)
                src <- if (file_ok(svg_path)) svg_path else p
                if (file_ok(src)) file.copy(src, file) else writeLines("", file)
            }
        )

        # Return selected region for parent module
        selected_region
    })
}
