#' PreGEA tab UI
#'
#' Hyperparameter-exploration mode (LD-pruned only): sweeps LFMM-K,
#' EMMAX-#PC, and RDA Condition()-PC ladders + a shared spatial
#' variance-partitioning panel, before committing to the expensive
#' full-SNP GEA run. Every plot here is grid-level (histogram grid, QQ
#' grid, lambda-vs-rung, hits-vs-rung) — plot_pregea_ladder.R never emits
#' one file per K/#PC rung, so unlike PreStructure/Structure there is no
#' per-rung selector to build; the ladder tables carry the per-rung numbers.
#'
#' @param id module namespace id
#' @noRd
mod_pregea_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Recommendations — headline handoff artifact, cross-cuts every ladder
        # below. Read-only here; pre-fill/Apply into the GEA config is a
        # separate follow-up (see project notes).
        bslib::card(
            full_screen = TRUE,
            bslib::card_header("Recommended Hyperparameters"),
            htmltools::p(class = "text-muted small",
                "One row per (method, hyperparameter): the rule that produced the ",
                "recommendation and the evidence behind it. Nothing here writes to ",
                "the config automatically — copy values you agree with into the ",
                "GEA tab's method config."),
            DT::DTOutput(ns("recommendations_table"))
        ),

        bslib::navset_card_tab(
            id = ns("pregea_tabs"),

            bslib::nav_panel("Climate & Predictors",
                shiny::uiOutput(ns("climate_invariant_warning")),
                mod_image_card_ui(ns("climate_heatmap")),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("climate_density")),
                    mod_image_card_ui(ns("phenotype_density"))
                )
            ),

            bslib::nav_panel("LFMM K ladder",
                mod_image_card_ui(ns("lfmm_screeplot")),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("lfmm_hist")),
                    mod_image_card_ui(ns("lfmm_qq")),
                    mod_image_card_ui(ns("lfmm_lambda")),
                    mod_image_card_ui(ns("lfmm_hits"))
                ),
                htmltools::h6("Ladder Table", class = "mt-3"),
                DT::DTOutput(ns("lfmm_table"))
            ),

            bslib::nav_panel("EMMAX #PC ladder",
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("emmax_hist")),
                    mod_image_card_ui(ns("emmax_qq")),
                    mod_image_card_ui(ns("emmax_lambda")),
                    mod_image_card_ui(ns("emmax_hits"))
                ),
                htmltools::h6("Ladder Table", class = "mt-3"),
                DT::DTOutput(ns("emmax_table"))
            ),

            bslib::nav_panel("RDA setup",
                bslib::layout_column_wrap(
                    width = 1 / 3,
                    mod_image_card_ui(ns("rda_collinearity")),
                    mod_image_card_ui(ns("rda_screeplot")),
                    mod_image_card_ui(ns("rda_condition_ladder"))
                ),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("rda_ordir2step")),
                    mod_image_card_ui(ns("rda_biplot"))
                ),
                htmltools::h6("Tables", class = "mt-3"),
                bslib::accordion(
                    open = FALSE, multiple = TRUE,
                    bslib::accordion_panel("Condition-PC ladder",
                        DT::DTOutput(ns("rda_ladder_table"))),
                    bslib::accordion_panel("Predictor collinearity",
                        DT::DTOutput(ns("rda_collinearity_table"))),
                    bslib::accordion_panel("ordiR2step selection path",
                        DT::DTOutput(ns("rda_ordir2step_table")))
                )
            ),

            bslib::nav_panel("Varpart & spatial (dbMEM)",
                bslib::layout_column_wrap(
                    width = 1 / 3,
                    mod_image_card_ui(ns("dbmem_screeplot")),
                    mod_image_card_ui(ns("dbmem_selection_path")),
                    mod_image_card_ui(ns("varpart_venn"))
                ),
                mod_image_card_ui(ns("px_barplot")),
                htmltools::h6("Tables", class = "mt-3"),
                bslib::accordion(
                    open = FALSE, multiple = TRUE,
                    bslib::accordion_panel("Variance-partition fractions",
                        DT::DTOutput(ns("varpart_fractions_table"))),
                    bslib::accordion_panel("Testable-fraction significance",
                        DT::DTOutput(ns("varpart_anova_table"))),
                    bslib::accordion_panel("Px per variable",
                        DT::DTOutput(ns("px_table"))),
                    bslib::accordion_panel("dbMEM diagnostics",
                        DT::DTOutput(ns("dbmem_diagnostics_table")))
                )
            )
        ),

        # Transfer guard — opt-in (PreGEA.TransferGuard.enabled), cross-cuts
        # LFMM + EMMAX, collapsed by default since most runs never enable it.
        bslib::accordion(
            open = FALSE,
            bslib::accordion_panel(
                "Transfer Guard (full-set validation)",
                icon = bsicons::bs_icon("shield-check"),
                mod_image_card_ui(ns("transfer_guard_plot")),
                htmltools::h6("Table", class = "mt-3"),
                DT::DTOutput(ns("transfer_guard_table"))
            )
        )
    )
}

#' PreGEA tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_pregea_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Small local helper: safe_datatable() over a TSV that may not exist yet.
        load_tsv_dt <- function(path) {
            if (!file_ok(path)) return(data.table::data.table())
            tryCatch(data.table::fread(path, sep = "\t", header = TRUE),
                     error = function(e) data.table::data.table())
        }
        render_tsv_table <- function(path) {
            safe_datatable(as.data.frame(load_tsv_dt(path)),
                           options = list(dom = "t", scrollX = TRUE, pageLength = -1),
                           rownames = FALSE)
        }

        # ── Recommendations ─────────────────────────────────────────────────
        output$recommendations_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "", "pregea_recommendations"))
        })

        # ── Climate & Predictors (display moved from Structure tab; producers
        #    stay in structure.smk — zero DAG/path change) ─────────────────────
        shiny::observe({
            pd <- project_data()
            mod_image_card_server("climate_heatmap",
                path    = shiny::reactive(climate_heatmap_path(pd$name)),
                title   = shiny::reactive("Climate Correlation Heatmap"),
                dl_name = shiny::reactive("climate_heatmap"),
                note    = shiny::reactive(help_note("climate_heatmap"))
            )
            mod_image_card_server("climate_density",
                path    = shiny::reactive(climate_density_path(pd$name)),
                title   = shiny::reactive("Climate Density (Present)"),
                dl_name = shiny::reactive("climate_density_present"),
                note    = shiny::reactive(help_note("climate_density"))
            )
            mod_image_card_server("phenotype_density",
                path    = shiny::reactive(phenotype_density_path(pd$name)),
                title   = shiny::reactive("Phenotype Density"),
                dl_name = shiny::reactive("phenotype_density"),
                note    = shiny::reactive(help_note("phenotype_density"))
            )

            # I3: Climate invariant predictors warning (moved verbatim from
            # mod_structure.R — was nested inside its pop-stats observe(), not
            # the climate one, despite belonging to the climate section).
            output$climate_invariant_warning <- shiny::renderUI({
                p <- climate_invariant_path(pd$name)
                if (!file.exists(p)) return(NULL)
                dt <- tryCatch(
                    data.table::fread(p, sep = "\t", header = TRUE),
                    error = function(e) data.table::data.table()
                )
                if (nrow(dt) == 0 || !"predictor" %in% names(dt)) return(NULL)
                preds <- dt$predictor
                items <- lapply(preds, function(pr) {
                    reason <- if ("reason" %in% names(dt)) {
                        r <- dt[dt$predictor == pr, ]$reason
                        if (length(r) > 0) paste0(" — ", r[1]) else ""
                    } else ""
                    htmltools::tags$li(htmltools::tags$code(pr), reason)
                })
                htmltools::div(
                    class = "d-flex justify-content-end mb-2",
                    filter_note(
                        paste0(length(preds), " excluded"),
                        htmltools::p(
                            htmltools::tags$strong(length(preds),
                                if (length(preds) == 1) "climate predictor" else "climate predictors",
                                "excluded due to zero variance at sample sites"),
                            " — not used in GEA or Gradient Forest analyses. ",
                            "Consider removing from ", htmltools::tags$code("Climate.predictors"), ".",
                            htmltools::tags$ul(class = "mb-0 mt-1", items)
                        ),
                        class = "bg-warning text-dark"
                    )
                )
            })
        })

        # ── LFMM K ladder ────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Enable PreGEA.LFMM.enabled and re-run PreGEA.")

            mod_image_card_server("lfmm_screeplot",
                path       = shiny::reactive(pregea_plot_path(pd$name, "structure", "pruned_pca_screeplot")),
                title      = shiny::reactive("LD-pruned PCA Screeplot"),
                dl_name    = shiny::reactive("pregea_pruned_pca_screeplot"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_screeplot"))
            )
            mod_image_card_server("lfmm_hist",
                path       = shiny::reactive(pregea_plot_path(pd$name, "lfmm", "lfmm_pvalue_histogram_grid")),
                title      = shiny::reactive("P-value Histogram (per K)"),
                dl_name    = shiny::reactive("lfmm_pvalue_histogram_grid"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_lfmm_hist"))
            )
            mod_image_card_server("lfmm_qq",
                path       = shiny::reactive(pregea_plot_path(pd$name, "lfmm", "lfmm_qq_grid")),
                title      = shiny::reactive("QQ Plot (per K)"),
                dl_name    = shiny::reactive("lfmm_qq_grid"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_lfmm_qq"))
            )
            mod_image_card_server("lfmm_lambda",
                path       = shiny::reactive(pregea_plot_path(pd$name, "lfmm", "lfmm_lambda_vs_K")),
                title      = shiny::reactive("Lambda vs K"),
                dl_name    = shiny::reactive("lfmm_lambda_vs_K"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_lfmm_lambda"))
            )
            mod_image_card_server("lfmm_hits",
                path       = shiny::reactive(pregea_plot_path(pd$name, "lfmm", "lfmm_hits_vs_K")),
                title      = shiny::reactive("Hit Count vs K"),
                dl_name    = shiny::reactive("lfmm_hits_vs_K"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_lfmm_hits"))
            )
        })
        output$lfmm_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "lfmm", "lfmm_ladder"))
        })

        # ── EMMAX #PC ladder ─────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Enable PreGEA.EMMAX.enabled and re-run PreGEA.")

            mod_image_card_server("emmax_hist",
                path       = shiny::reactive(pregea_plot_path(pd$name, "emmax", "emmax_pvalue_histogram_grid")),
                title      = shiny::reactive("P-value Histogram (per #PCs)"),
                dl_name    = shiny::reactive("emmax_pvalue_histogram_grid"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_emmax_hist"))
            )
            mod_image_card_server("emmax_qq",
                path       = shiny::reactive(pregea_plot_path(pd$name, "emmax", "emmax_qq_grid")),
                title      = shiny::reactive("QQ Plot (per #PCs)"),
                dl_name    = shiny::reactive("emmax_qq_grid"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_emmax_qq"))
            )
            mod_image_card_server("emmax_lambda",
                path       = shiny::reactive(pregea_plot_path(pd$name, "emmax", "emmax_lambda_vs_n_pcs")),
                title      = shiny::reactive("Lambda vs #PCs"),
                dl_name    = shiny::reactive("emmax_lambda_vs_n_pcs"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_emmax_lambda"))
            )
            mod_image_card_server("emmax_hits",
                path       = shiny::reactive(pregea_plot_path(pd$name, "emmax", "emmax_hits_vs_n_pcs")),
                title      = shiny::reactive("Hit Count vs #PCs"),
                dl_name    = shiny::reactive("emmax_hits_vs_n_pcs"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_emmax_hits"))
            )
        })
        output$emmax_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "emmax", "emmax_ladder"))
        })

        # ── RDA setup ────────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Enable PreGEA.RDA.enabled and re-run PreGEA.")

            mod_image_card_server("rda_collinearity",
                path       = shiny::reactive(pregea_plot_path(pd$name, "rda", "rda_predictor_collinearity")),
                title      = shiny::reactive("Predictor Collinearity"),
                dl_name    = shiny::reactive("rda_predictor_collinearity"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_rda_collinearity"))
            )
            mod_image_card_server("rda_screeplot",
                path       = shiny::reactive(pregea_plot_path(pd$name, "rda", "rda_axis_screeplot")),
                title      = shiny::reactive("Constrained-axis Screeplot"),
                dl_name    = shiny::reactive("rda_axis_screeplot"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_rda_screeplot"))
            )
            mod_image_card_server("rda_condition_ladder",
                path       = shiny::reactive(pregea_plot_path(pd$name, "rda", "rda_condition_ladder")),
                title      = shiny::reactive("Condition-PC Ladder"),
                dl_name    = shiny::reactive("rda_condition_ladder"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_rda_condition_ladder"))
            )
            mod_image_card_server("rda_ordir2step",
                path       = shiny::reactive(pregea_plot_path(pd$name, "rda", "rda_ordir2step_path")),
                title      = shiny::reactive("ordiR2step Selection Path"),
                dl_name    = shiny::reactive("rda_ordir2step_path"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_rda_ordir2step"))
            )
            mod_image_card_server("rda_biplot",
                path       = shiny::reactive(pregea_plot_path(pd$name, "rda", "rda_biplot_best")),
                title      = shiny::reactive("RDA Biplot (recommended rung)"),
                dl_name    = shiny::reactive("rda_biplot_best"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_rda_biplot"))
            )
        })
        output$rda_ladder_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "rda", "rda_condition_ladder"))
        })
        output$rda_collinearity_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "rda", "rda_predictor_collinearity"))
        })
        output$rda_ordir2step_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "rda", "rda_ordir2step_path"))
        })

        # ── Varpart & spatial (dbMEM) ────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Enable PreGEA.Varpart.enabled and re-run PreGEA.")

            mod_image_card_server("dbmem_screeplot",
                path       = shiny::reactive(pregea_plot_path(pd$name, "spatial", "dbmem_screeplot")),
                title      = shiny::reactive("dbMEM Screeplot"),
                dl_name    = shiny::reactive("dbmem_screeplot"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_dbmem_screeplot"))
            )
            mod_image_card_server("dbmem_selection_path",
                path       = shiny::reactive(pregea_plot_path(pd$name, "varpart", "dbmem_selection_path")),
                title      = shiny::reactive("dbMEM Forward-selection Path"),
                dl_name    = shiny::reactive("dbmem_selection_path"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_dbmem_selection_path"))
            )
            mod_image_card_server("varpart_venn",
                path       = shiny::reactive(pregea_plot_path(pd$name, "varpart", "varpart_venn")),
                title      = shiny::reactive("Variance-partition Venn"),
                dl_name    = shiny::reactive("varpart_venn"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_varpart_venn"))
            )
            mod_image_card_server("px_barplot",
                path       = shiny::reactive(pregea_plot_path(pd$name, "varpart", "px_barplot")),
                title      = shiny::reactive("Px per Predictor (Lasky et al. 2012)"),
                dl_name    = shiny::reactive("px_barplot"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("pregea_px_barplot"))
            )
        })
        output$varpart_fractions_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "varpart", "varpart_fractions"))
        })
        output$varpart_anova_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "varpart", "varpart_anova"))
        })
        output$px_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "varpart", "px_per_variable"))
        })
        output$dbmem_diagnostics_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "spatial", "dbmem_diagnostics"))
        })

        # ── Transfer guard (opt-in) ──────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            mod_image_card_server("transfer_guard_plot",
                path       = shiny::reactive(pregea_plot_path(pd$name, "transfer", "transfer_guard_lambda")),
                title      = shiny::reactive("Transfer Guard: Lambda Pruned vs Full"),
                dl_name    = shiny::reactive("pregea_transfer_guard_lambda"),
                suggestion = shiny::reactive("Enable PreGEA.TransferGuard.enabled and re-run PreGEA."),
                note       = shiny::reactive(help_note("pregea_transfer_guard"))
            )
        })
        output$transfer_guard_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "", "pregea_transfer_guard"))
        })
    })
}
