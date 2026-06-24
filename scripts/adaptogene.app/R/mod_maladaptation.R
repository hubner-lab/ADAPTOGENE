#' Maladaptation tab UI
#'
#' Gradient Forest results: importance plots, genetic offset piemaps,
#' zoom maps, and site-level offset table.
#'
#' @param id module namespace id
#' @noRd
mod_maladaptation_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # SNP set / GF run selector + badges + set actions
        htmltools::div(
            class = "control-bar",
            shiny::uiOutput(ns("suffix_selector")),
            shiny::uiOutput(ns("snp_count_badge")),
            shiny::uiOutput(ns("gf_params_badge")),
            shiny::uiOutput(ns("climate_scenario_badge")),
            shiny::uiOutput(ns("set_details_btn")),
            shiny::uiOutput(ns("set_delete_btn"))
        ),

        # Importance plots
        bslib::layout_column_wrap(
            width = 1 / 2,
            mod_image_card_ui(ns("overall_importance")),
            mod_image_card_ui(ns("cumulative_importance"))
        ),

        shiny::hr(),

        # Piemap Type + Zoom right before piemap
        htmltools::div(
            class = "control-bar",
            bslib::layout_columns(
                col_widths = c(6, 6),
                shiny::uiOutput(ns("piemap_variant_ui")),
                shiny::uiOutput(ns("zoom_selector"))
            )
        ),

        # Genetic offset piemap (constrained to 2/3 width; full-screen for detail)
        bslib::layout_columns(
            col_widths = c(8, 4),
            htmltools::div(
                class = "piemap-container",
                mod_image_card_ui(ns("offset_piemap"))
            ),
            bslib::card(
                bslib::card_body(
                    class = "text-muted small",
                    htmltools::p(bsicons::bs_icon("info-circle"), " Genetic offset estimates local population vulnerability to projected climate change."),
                    htmltools::p(bsicons::bs_icon("arrows-fullscreen"), " Use the expand icon for a full-resolution view of the map."),
                    htmltools::p(bsicons::bs_icon("table"), " Site-level offset values are in the table below.")
                )
            )
        ),

        # Site table in accordion
        bslib::accordion(
            id       = ns("malad_sections"),
            open     = FALSE,
            multiple = TRUE,

            bslib::accordion_panel(
                "Site Offset Table",
                value = "site_table",
                icon  = bsicons::bs_icon("table"),
                bslib::card_body(
                    shiny::uiOutput(ns("offset_summary")),
                    shiny::downloadButton(ns("dl_site"), "Download CSV",
                                          class = "btn-sm btn-outline-secondary mb-2"),
                    DT::DTOutput(ns("site_table"))
                )
            ),

            bslib::accordion_panel(
                "Cross-set comparison",
                value = "cross_set",
                icon  = bsicons::bs_icon("bar-chart-steps"),
                bslib::card_body(
                    htmltools::div(
                        class = "text-muted small fst-italic",
                        "Coming soon: side-by-side comparison of genetic offset across saved SNP sets."
                    )
                )
            )
        )
    )
}

#' Maladaptation tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_maladaptation_server <- function(id, project_data, snp_sets_trigger = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # gf_suffixes depends on both project_data AND snp_sets_trigger
        # (project_data refreshes on pipeline run; trigger refreshes on delete)
        gf_suffixes <- shiny::reactive({
            if (!is.null(snp_sets_trigger)) snp_sets_trigger()
            find_gf_suffixes(project_data()$name)
        })

        output$suffix_selector <- shiny::renderUI({
            suf <- gf_suffixes()
            if (length(suf) == 0)
                return(shiny::p(
                    "No GF runs found. Save a SNP set in the GEA tab and run Maladaptation.",
                    class = "text-muted small"))
            shiny::selectInput(ns("suffix"), "SNP set / GF run", choices = suf, selected = suf[1])
        })

        selected_suffix <- shiny::reactive({
            s <- input$suffix
            if (is.null(s)) {
                suf <- gf_suffixes()
                if (length(suf) == 0) return(NULL)
                suf[1]
            } else s
        })

        # ── Resolve suffix → set name via manifest (NOT string-stripping) ─────
        # Pipeline emits {name}_spatial | {name}_nospatial | {name} (without spatial).
        # Underscores + 'both' mode make suffix-stripping ambiguous.
        selected_set_name <- shiny::reactive({
            suf <- selected_suffix(); if (is.null(suf)) return(NULL)
            if (!is.null(snp_sets_trigger)) snp_sets_trigger()
            names <- list_snp_sets(project_data()$name)$name
            # Longest match: suf == name OR suf == name_spatial OR suf == name_nospatial
            hit <- names[suf == names |
                         suf == paste0(names, "_spatial") |
                         suf == paste0(names, "_nospatial")]
            if (length(hit) == 0) NULL else hit[which.max(nchar(hit))]
        })

        # ── Manifest entry for the selected set ───────────────────────────────
        selected_set_entry <- shiny::reactive({
            nm <- selected_set_name(); if (is.null(nm)) return(NULL)
            man <- read_snp_sets_manifest(project_data()$name)
            entries <- Filter(function(x) identical(x$name, nm), man)
            if (length(entries) == 0) NULL else entries[[1]]
        })

        # ── View-details popover ───────────────────────────────────────────────
        output$set_details_btn <- shiny::renderUI({
            e <- selected_set_entry(); if (is.null(e)) return(NULL)
            body <- htmltools::tags$ul(
                class = "small mb-0 ps-3",
                htmltools::tags$li(htmltools::strong("Threshold: "),
                    paste0(e$threshold_type %||% "?", " = ", e$threshold_value %||% "?")),
                htmltools::tags$li(htmltools::strong("Strategy: "),  e$strategy %||% "?"),
                htmltools::tags$li(htmltools::strong("Regime: "),    e$regime %||% "?"),
                htmltools::tags$li(htmltools::strong("Source: "),    e$source_module %||% "?"),
                htmltools::tags$li(htmltools::strong("N SNPs: "),    e$n_snps %||% "?"),
                htmltools::tags$li(htmltools::strong("Created: "),   e$created %||% "?")
            )
            bslib::popover(
                shiny::actionButton(ns("set_details"),
                    label = htmltools::tagList(bsicons::bs_icon("info-circle"), " View details"),
                    class = "btn btn-outline-secondary btn-sm ms-2"),
                title = paste("SNP set:", e$name %||% ""),
                body
            )
        })

        # ── Delete button + confirm modal ──────────────────────────────────────
        output$set_delete_btn <- shiny::renderUI({
            if (is.null(selected_set_name())) return(NULL)
            shiny::actionButton(ns("set_delete"),
                label = htmltools::tagList(bsicons::bs_icon("trash3"), " Delete"),
                class = "btn btn-outline-danger btn-sm ms-1")
        })

        shiny::observeEvent(input$set_delete, {
            nm <- selected_set_name(); if (is.null(nm)) return()
            shiny::showModal(shiny::modalDialog(
                title = "Delete SNP set",
                htmltools::p(
                    "Delete set ", htmltools::strong(nm),
                    " and its Gradient Forest results? This cannot be undone."
                ),
                footer = htmltools::tagList(
                    shiny::modalButton("Cancel"),
                    shiny::actionButton(ns("set_delete_confirm"), "Delete",
                                        class = "btn btn-danger")
                ),
                easyClose = TRUE
            ))
        })

        shiny::observeEvent(input$set_delete_confirm, {
            nm <- selected_set_name()
            pd <- project_data()
            if (is.null(nm) || is.null(pd)) { shiny::removeModal(); return() }
            delete_snp_set(pd$name, nm, remove_gf_results = TRUE)
            shiny::removeModal()
            shiny::showNotification(
                sprintf("Deleted SNP set '%s'.", nm),
                type = "message", duration = 4)
            if (!is.null(snp_sets_trigger))
                snp_sets_trigger(snp_sets_trigger() + 1L)
        })

        # ── G1: SNP count badge with denominator ──────────────────────────────
        output$snp_count_badge <- shiny::renderUI({
            suf <- selected_suffix()
            pd  <- project_data()
            if (is.null(suf) || is.null(pd)) return(NULL)
            # Look up via curated set name (manifest-resolved)
            set_nm    <- selected_set_name()
            snps_path <- if (!is.null(set_nm)) snp_set_path(pd$name, set_nm) else NULL
            if (is.null(snps_path) || !file_ok(snps_path)) return(NULL)
            n <- tryCatch({
                nrow(data.table::fread(snps_path, select = 1L))
            }, error = function(e) NULL)
            if (is.null(n)) return(NULL)
            # Get total filtered SNP count for context
            summ <- load_pipeline_summary(pd$name)
            n_total_row <- summ[summ$step == "processing" & summ$metric == "snps_after_filtering", ]
            n_total <- if (nrow(n_total_row) > 0) suppressWarnings(as.integer(n_total_row$value[1])) else NA_integer_
            label <- if (!is.na(n_total) && n_total > 0) {
                pct <- round(100 * n / n_total, 1)
                paste0("\u25c6 ", n, " / ", format(n_total, big.mark = ","), " adaptive SNPs (", pct, "%)")
            } else {
                paste0("\u25c6 ", n, " adaptive SNPs")
            }
            htmltools::tags$span(
                class = "badge bg-secondary ms-2 align-self-center",
                style = "font-size:0.85rem; vertical-align:middle;",
                label
            )
        })

        # ── G2: GF run parameters badge ───────────────────────────────────────
        output$gf_params_badge <- shiny::renderUI({
            suf <- selected_suffix()
            pd  <- project_data()
            if (is.null(suf) || is.null(pd)) return(NULL)
            cfg     <- pd$config
            ntree   <- config_get(cfg, "Maladaptation", "methods", "gradient_forest",
                                  "ntree", default = 500L)
            cor_thr <- config_get(cfg, "Maladaptation", "methods", "gradient_forest",
                                  "cor_threshold", default = 0.5)
            spatial <- config_get(cfg, "Maladaptation", "methods", "gradient_forest",
                                  "spatial_correction", default = "with")
            rand    <- config_get(cfg, "Maladaptation", "methods", "gradient_forest",
                                  "random_model", default = TRUE)
            htmltools::div(
                class = "ms-2 align-self-center",
                config_badges_bar(
                    config_badge("ntree", ntree),
                    config_badge("cor.thr", cor_thr),
                    config_badge("spatial", spatial),
                    if (isTRUE(rand)) config_badge("random model", "yes") else NULL
                )
            )
        })

        # ── G4: Future climate scenario badge ─────────────────────────────────
        output$climate_scenario_badge <- shiny::renderUI({
            suf <- selected_suffix()
            pd  <- project_data()
            if (is.null(suf) || is.null(pd)) return(NULL)
            cfg    <- pd$config
            ssp    <- config_get(cfg, "Future", "ssp",    default = NULL)
            year   <- config_get(cfg, "Future", "year",   default = NULL)
            models <- config_get(cfg, "Future", "models", default = NULL)
            if (is.null(ssp) && is.null(year)) return(NULL)
            n_models <- if (is.null(models)) 0L else length(models)
            label_parts <- c(
                if (!is.null(ssp))  paste0("SSP", ssp),
                if (!is.null(year)) as.character(year),
                if (n_models > 0)   paste0(n_models, " GCMs")
            )
            htmltools::div(
                class = "ms-1 align-self-center",
                htmltools::tags$span(
                    class = "badge bg-info text-dark align-self-center",
                    style = "font-size:0.75rem;",
                    paste0("\u2601 Future: ", paste(label_parts, collapse = " "))
                )
            )
        })

        # ── G3: Offset range summary ──────────────────────────────────────────
        output$offset_summary <- shiny::renderUI({
            dt <- site_data()
            if (is.null(dt) || nrow(dt) == 0 || !"genetic_offset" %in% names(dt))
                return(NULL)
            offs <- dt$genetic_offset
            offs <- offs[!is.na(offs)]
            if (length(offs) == 0) return(NULL)
            htmltools::div(
                class = "alert alert-info d-flex gap-2 align-items-start mb-2 py-2",
                bsicons::bs_icon("info-circle-fill", class = "flex-shrink-0 mt-1"),
                htmltools::div(
                    htmltools::tags$strong("Genetic offset range: "),
                    sprintf("%.4f", min(offs)), " \u2013 ", sprintf("%.4f", max(offs)),
                    " (mean: ", sprintf("%.4f", mean(offs)), ", n=", length(offs), " sites)"
                )
            )
        })

        # ── Piemap variant selector (only shows variants that exist) ───────────
        output$piemap_variant_ui <- shiny::renderUI({
            suf <- shiny::req(selected_suffix())
            pd  <- project_data()
            choices <- c("Genetic Offset" = "base")
            if (file_ok(gf_offset_piemap_path(pd$name, suf, "tajima_d")))
                choices <- c(choices, c("Offset \u00d7 Tajima's D" = "tajima_d"))
            if (file_ok(gf_offset_piemap_path(pd$name, suf, "pi_diversity")))
                choices <- c(choices, c("Offset \u00d7 Pi Diversity" = "pi_diversity"))
            if (length(choices) == 1) return(NULL)  # only base — no selector needed
            shiny::selectInput(ns("piemap_variant"), "Piemap Type",
                               choices = choices, selected = "base")
        })

        # ── Zoom selector (only shows if zoom maps exist) ──────────────────────
        output$zoom_selector <- shiny::renderUI({
            suf   <- shiny::req(selected_suffix())
            pd    <- project_data()
            zooms <- find_gf_zooms(pd$name, suf)
            if (length(zooms) == 0) return(NULL)
            shiny::selectInput(ns("zoom"), "Zoom Region",
                choices  = c("Full view" = "none", setNames(zooms, zooms)),
                selected = "none"
            )
        })

        selected_variant <- shiny::reactive(input$piemap_variant %||% "base")
        selected_zoom    <- shiny::reactive(input$zoom %||% "none")

        # ── Importance images ──────────────────────────────────────────────────
        shiny::observe({
            suf <- shiny::req(selected_suffix())
            pd  <- project_data()

            mod_image_card_server("overall_importance",
                path    = shiny::reactive(gf_importance_path(pd$name, suf, "overall")),
                title   = shiny::reactive("Overall Variable Importance"),
                dl_name = shiny::reactive(paste0("overall_importance_", suf))
            )
            mod_image_card_server("cumulative_importance",
                path    = shiny::reactive(gf_importance_path(pd$name, suf, "cumulative")),
                title   = shiny::reactive("Cumulative Importance"),
                dl_name = shiny::reactive(paste0("cumulative_importance_", suf))
            )
        })

        # ── Single selector-driven offset piemap ───────────────────────────────
        piemap_path <- shiny::reactive({
            suf     <- shiny::req(selected_suffix())
            pd      <- project_data()
            variant <- selected_variant()
            zoom    <- selected_zoom()
            if (zoom != "none") {
                mod_path(pd$name, MOD_MALAD, "plots", suf, "zoom",
                         paste0(zoom, ".png"))
            } else {
                gf_offset_piemap_path(pd$name, suf, variant)
            }
        })

        piemap_title <- shiny::reactive({
            variant <- selected_variant()
            zoom    <- selected_zoom()
            label   <- c(
                base         = "Genetic Offset",
                tajima_d     = "Offset \u00d7 Tajima's D",
                pi_diversity = "Offset \u00d7 Pi Diversity"
            )[variant]
            if (zoom != "none") paste0(label, " \u2014 Zoom: ", zoom) else label
        })

        mod_image_card_server("offset_piemap",
            path        = piemap_path,
            title       = piemap_title,
            dl_name     = shiny::reactive(paste0("offset_piemap_", selected_variant(),
                                                  "_", selected_suffix() %||% "gf")),
            placeholder = shiny::reactive("Run mode=maladaptation to generate Gradient Forest results")
        )

        # ── Site table ─────────────────────────────────────────────────────────
        site_data <- shiny::reactive({
            suf <- selected_suffix()
            if (is.null(suf)) return(data.table::data.table())
            p  <- gf_site_table_path(project_data()$name, suf)
            # mtime fingerprint: re-running GF for the same suffix otherwise serves stale data
            fp <- if (file.exists(p)) as.character(file.info(p)$mtime) else "missing"
            load_cached(paste0("gf_site_", project_data()$name, "_", suf),
                        function() {
                p <- gf_site_table_path(project_data()$name, suf)
                if (!file_ok(p)) return(data.table::data.table())
                dt <- data.table::fread(p, colClasses = c("site" = "character",
                                                            "sample" = "character"))
                # Join pop stats (TajimaD, PI) if available
                tajima_p <- mod_path(project_data()$name, MOD_STRUCT,
                                     "tables", "pop_stats", "tajima_d_by_pop.tsv")
                pi_p     <- mod_path(project_data()$name, MOD_STRUCT,
                                     "tables", "pop_stats", "pi_diversity_by_pop.tsv")
                if (file_ok(tajima_p)) {
                    taj <- data.table::fread(tajima_p, colClasses = c("site" = "character"))
                    dt  <- taj[dt, on = "site"]
                }
                if (file_ok(pi_p)) {
                    pi_dt <- data.table::fread(pi_p, colClasses = c("site" = "character"))
                    dt    <- pi_dt[dt, on = "site"]
                }
                dt
            }, fingerprint = fp)
        })

        output$site_table <- DT::renderDataTable({
            dt <- site_data()
            safe_datatable(
                as.data.frame(dt),
                extensions = "Buttons",
                options    = list(
                    dom       = "Bfrtip",
                    buttons   = list("csv"),
                    scrollX   = TRUE,
                    pageLength = 20
                ),
                rownames = FALSE
            )
        })

        output$dl_site <- shiny::downloadHandler(
            filename = function() {
                paste0("genetic_offset_site_", project_data()$name,
                       "_", selected_suffix(), ".csv")
            },
            content = function(file) {
                dt <- site_data()
                utils::write.csv(as.data.frame(dt), file, row.names = FALSE)
            }
        )
    })
}
