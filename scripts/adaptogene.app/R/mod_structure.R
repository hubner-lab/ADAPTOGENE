#' Structure tab UI
#'
#' Piemaps (bio variable + metric + zoom), climate plots, population stats,
#' LD decay, and summary tables.
#'
#' @param id module namespace id
#' @noRd
mod_structure_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Piemap controls inline above piemap
        htmltools::div(
            class = "control-bar",
            bslib::layout_columns(
                col_widths = c(3, 3, 3, 3),
                shiny::uiOutput(ns("bio_selector")),
                shiny::uiOutput(ns("metric_selector")),
                shiny::uiOutput(ns("zoom_selector")),
                htmltools::div(
                    class = "d-flex align-items-center h-100",
                    bslib::input_switch(ns("points"), "Points", value = FALSE)
                )
            )
        ),

        # Piemap + PCA-structure companion (constrains piemap to 2/3 width)
        bslib::layout_columns(
            col_widths = c(8, 4),
            mod_piemap_viewer_ui(ns("piemap")),
            mod_image_card_ui(ns("pca_structure_k"))
        ),

        # Climate correlation/density plots moved to the Climate tab (display
        # only — producers stay in structure.smk, zero DAG/path change).
        htmltools::div(
            class = "d-flex align-items-center gap-2 text-muted small mb-2",
            bsicons::bs_icon("thermometer-half"),
            "Climate correlation, density, and the invariant-predictor warning moved to the ",
            htmltools::tags$strong("Climate"), " tab (Predictors panel)."
        ),

        # LD decay, pop stats, tables in accordion
        bslib::accordion(
            id       = ns("sk_sections"),
            open     = c("ld_decay"),
            multiple = TRUE,

            bslib::accordion_panel(
                "LD Decay",
                value = "ld_decay",
                icon  = bsicons::bs_icon("bar-chart-line"),
                shiny::uiOutput(ns("ld_decay_skipped_alert")),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("ld_decay")),
                    mod_image_card_ui(ns("ld_decay_chr"))
                ),
                htmltools::h6("LD Decay Summary", class = "mt-3"),
                DT::DTOutput(ns("ld_decay_table"))
            ),

            bslib::accordion_panel(
                "Population Statistics",
                value = "pop_stats",
                icon  = bsicons::bs_icon("graph-up"),
                shiny::uiOutput(ns("pop_stats_info")),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("mantel")),
                    mod_image_card_ui(ns("amova"))
                )
            ),

            bslib::accordion_panel(
                "Tables",
                value  = "tables",
                icon   = bsicons::bs_icon("table"),
                shiny::uiOutput(ns("tables_content"))
            )
        )
    )
}

#' Structure K tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_structure_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        bio_values <- shiny::reactive(find_bio_values(project_data()$name))

        # ── Sidebar selectors ──────────────────────────────────────────────────
        output$bio_selector <- shiny::renderUI({
            available <- bio_values()
            all_bios  <- 1:19
            labels <- ifelse(
                all_bios %in% available,
                paste0("bio", all_bios),
                paste0("bio", all_bios, " (unavailable)")
            )
            choices  <- setNames(all_bios, labels)
            selected <- if (length(available) > 0) available[1] else 1L
            shiny::selectInput(ns("bio"), "Bio variable", choices = choices, selected = selected)
        })

        output$metric_selector <- shiny::renderUI({
            shiny::selectInput(ns("metric"), "Scale by metric",
                choices = c(
                    "None"           = "none",
                    "Tajima's D"     = "tajima_d",
                    "Pi Diversity"   = "pi_diversity"
                ),
                selected = "none"
            )
        })

        output$zoom_selector <- shiny::renderUI({
            # Zoom directories are subdirectories of piemap/zoom/
            pd       <- project_data()
            zoom_dir <- mod_path(pd$name, MOD_STRUCT, "plots", "piemap", "zoom")
            if (!dir.exists(zoom_dir)) return(NULL)
            tags <- list.dirs(zoom_dir, full.names = FALSE, recursive = FALSE)
            tags <- tags[nzchar(tags)]
            if (length(tags) == 0) return(NULL)
            shiny::selectInput(ns("zoom"), "Zoom region",
                choices  = c("Global" = "none", setNames(tags, tags)),
                selected = "none"
            )
        })

        selected_bio <- shiny::reactive({
            b <- input$bio
            if (is.null(b)) {
                bios <- bio_values()
                if (length(bios) == 0) return(1L)
                bios[1]
            } else as.integer(b)
        })

        selected_metric <- shiny::reactive(input$metric %||% "none")
        selected_zoom   <- shiny::reactive(input$zoom   %||% "none")
        selected_points <- shiny::reactive(isTRUE(input$points))

        # ── Piemap viewer ──────────────────────────────────────────────────────
        mod_piemap_viewer_server("piemap",
            project_data = project_data,
            bio          = selected_bio,
            metric       = selected_metric,
            zoom         = selected_zoom,
            points       = selected_points,
            note         = shiny::reactive(help_note("structure_piemap"))
        )

        # ── PCA-structure companion (k_best) ───────────────────────────────────
        shiny::observe({
            pd <- project_data()
            k  <- pd$k_best
            if (is.na(k)) return()
            mod_image_card_server("pca_structure_k",
                path    = shiny::reactive(pca_structure_path(pd$name, k)),
                title   = shiny::reactive("PCA"),
                dl_name = shiny::reactive(paste0("pca_structure_K", k)),
                note    = shiny::reactive(help_note("pca_structure_k", label = paste0("K=", k)))
            )
        })

        # ── Population stats images ────────────────────────────────────────────
        shiny::observe({
            pd       <- project_data()
            stats_dir <- mod_path(pd$name, MOD_STRUCT, "plots", "pop_stats")

            mantel_p <- file.path(stats_dir, "mantel_test.png")
            amova_p  <- file.path(stats_dir, "amova.png")

            mod_image_card_server("mantel",
                path    = shiny::reactive(if (file_ok(mantel_p)) mantel_p else NULL),
                title   = shiny::reactive("Mantel Test"),
                dl_name = shiny::reactive("mantel_test"),
                note    = shiny::reactive(help_note("mantel"))
            )
            mod_image_card_server("amova",
                path    = shiny::reactive(if (file_ok(amova_p)) amova_p else NULL),
                title   = shiny::reactive("AMOVA"),
                dl_name = shiny::reactive("amova"),
                note    = shiny::reactive(help_note("amova"))
            )
            # C3: Population statistics sample counts
            output$pop_stats_info <- shiny::renderUI({
                k <- pd$k_best
                if (is.na(k)) return(NULL)
                # Q-matrix has no discrete cluster column \u2014 derive pop counts via
                # argmax over C1..CK (same helper PreStructure uses)
                summary <- cluster_pop_summary(load_clusters(pd$name, k), k)
                if (is.null(summary)) return(NULL)
                rng <- if (summary$min_n == summary$max_n) as.character(summary$min_n)
                       else paste0(summary$min_n, "\u2013", summary$max_n)
                htmltools::div(
                    class = "d-flex justify-content-end mb-2",
                    filter_note(
                        paste0(summary$n_pops, " pops"),
                        htmltools::p(
                            htmltools::tags$strong(summary$n_pops, "populations"),
                            " at K=", k, " \u2014 ",
                            rng, " samples per population."
                        ),
                        class = "bg-secondary"
                    )
                )
            })

            # LD decay skipped-groups alert
            output$ld_decay_skipped_alert <- shiny::renderUI({
                log_path <- file.path(get_pipeline_path(), paste0(pd$name, "_logs"),
                                      "structure", "ld_decay_prepare.log")
                if (!file.exists(log_path)) return(NULL)
                lines    <- readLines(log_path, warn = FALSE)

                # Single-group dataset: the lone group IS 'All', so there is no
                # threshold to adjust. Distinct message, distinct badge.
                single <- grep("^WARNING: Single group", lines, value = TRUE)
                if (length(single) > 0) {
                    grp <- gsub("WARNING: Single group '([^']+)'.*", "\\1", single[1])
                    return(htmltools::div(
                        class = "d-flex justify-content-end mb-2",
                        filter_note(
                            "single group",
                            htmltools::p(
                                htmltools::tags$strong("Per-group LD decay skipped"),
                                " — every sample belongs to the same group (",
                                htmltools::tags$code(grp),
                                "), which makes it identical to ",
                                htmltools::tags$code("All"),
                                ". The genome-wide curve is the whole dataset."
                            ),
                            class = "bg-info text-dark"
                        )
                    ))
                }

                skipped  <- grep("^WARNING: Skipping group", lines, value = TRUE)
                if (length(skipped) == 0) return(NULL)
                # Parse "WARNING: Skipping group 'X' (N samples < min_samples=M)"
                groups <- gsub("WARNING: Skipping group '([^']+)'.*", "\\1", skipped)
                reasons <- gsub("WARNING: Skipping group '[^']+' \\((.+)\\)", "\\1", skipped)
                items <- mapply(function(g, r) {
                    htmltools::tags$li(htmltools::tags$strong(g), paste0(" — ", r))
                }, groups, reasons, SIMPLIFY = FALSE)
                htmltools::div(
                    class = "d-flex justify-content-end mb-2",
                    filter_note(
                        paste0(length(skipped), " skipped"),
                        htmltools::p(
                            htmltools::tags$strong(
                                length(skipped),
                                if (length(skipped) == 1) "group was" else "groups were",
                                "excluded from per-group LD decay"
                            ),
                            " — below the ",
                            htmltools::tags$code("LDdecay.min_samples"),
                            " threshold. Adjust the config to include them.",
                            htmltools::tags$ul(class = "mb-0 mt-1", items)
                        ),
                        class = "bg-warning text-dark"
                    )
                )
            })

            ld_r2_note <- shiny::reactive({
                filter_note(
                    "r² = 0.2",
                    paste0(
                        "Dashed line marks the r² = 0.2 background-LD reference level. ",
                        "Half-decay and r²=0.2 distances per group are tabulated below — ",
                        "use them to pick a clumping distance in the GEA/GWAS filter bar."
                    )
                )
            })

            mod_image_card_server("ld_decay",
                path    = shiny::reactive(ld_decay_path(pd$name, per_chr = FALSE)),
                title   = shiny::reactive("LD Decay (Genome-wide)"),
                dl_name = shiny::reactive("ld_decay"),
                note    = ld_r2_note
            )
            mod_image_card_server("ld_decay_chr",
                path    = shiny::reactive(ld_decay_path(pd$name, per_chr = TRUE)),
                title   = shiny::reactive("LD Decay (Per Chromosome)"),
                dl_name = shiny::reactive("ld_decay_per_chr"),
                note    = ld_r2_note
            )
        })

        # ── LD decay summary table (half-decay + r2=0.2 distance per group/scope) ──
        output$ld_decay_table <- DT::renderDataTable({
            pd  <- project_data()
            tbl <- as.data.frame(load_ld_decay_table(pd$name))
            if (nrow(tbl) > 0) {
                tbl$half_decay_kb <- round(tbl$half_decay_bp / 1000, 1)
                tbl$r2_02_kb      <- round(tbl$r2_02_bp / 1000, 1)
                tbl <- tbl[, c("group", "scope", "half_decay_kb", "r2_02_kb",
                               "n_samples", "n_pairs", "r2_intercept", "method")]
                names(tbl) <- c("Group", "Scope", "Half-decay (kb)", "r²=0.2 distance (kb)",
                                 "N samples", "N pairs", "r² intercept", "Method")
            }
            safe_datatable(
                tbl,
                options  = list(dom = "t", scrollX = TRUE, pageLength = -1),
                rownames = FALSE
            )
        })

        # ── Tables ─────────────────────────────────────────────────────────────
        output$tables_content <- shiny::renderUI({
            pd        <- project_data()
            tables_dir <- mod_path(pd$name, MOD_STRUCT, "tables", "pop_stats")
            if (!dir.exists(tables_dir)) {
                return(htmltools::p("No population statistics tables available.",
                                    class = "text-muted small"))
            }
            files <- list.files(tables_dir, pattern = "\\.tsv$", full.names = TRUE)
            if (length(files) == 0) {
                return(htmltools::p("No tables available.", class = "text-muted small"))
            }

            # Build accordion panels for each table file
            panels <- lapply(files, function(f) {
                tbl_name <- tools::file_path_sans_ext(basename(f))
                tbl_id   <- gsub("[^A-Za-z0-9]", "_", tbl_name)
                bslib::accordion_panel(
                    gsub("_", " ", tbl_name),
                    DT::DTOutput(ns(paste0("tbl_", tbl_id)))
                )
            })

            # Render each table server-side
            lapply(files, function(f) {
                tbl_name <- tools::file_path_sans_ext(basename(f))
                tbl_id   <- gsub("[^A-Za-z0-9]", "_", tbl_name)
                output[[paste0("tbl_", tbl_id)]] <- DT::renderDataTable({
                    dt <- data.table::fread(f)
                    safe_datatable(
                        as.data.frame(dt),
                        extensions = "Buttons",
                        options    = list(
                            dom       = "Bfrtip",
                            buttons   = list("csv"),
                            scrollX   = TRUE,
                            pageLength = 15
                        ),
                        rownames = FALSE
                    )
                })
            })

            do.call(bslib::accordion, c(panels, list(open = FALSE, multiple = TRUE)))
        })
    })
}
