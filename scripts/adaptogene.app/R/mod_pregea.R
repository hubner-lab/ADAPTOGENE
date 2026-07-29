#' PreGEA tab UI
#'
#' Focused hyperparameter-choice module: "how many K / how many PCs" for
#' LFMM/EMMAX/GAPIT, plus RDA's Condition()-PC setup — sweeps LFMM-K,
#' EMMAX-#PC, and RDA Condition()-PC ladders (one shared PC range for the
#' latter two) before committing to the expensive full-SNP GEA run.
#' Predictor characterization (correlations, densities, varpart, dbMEM) lives
#' on the Climate tab (mod_climate.R) — see that file's header for the
#' module-split rationale.
#'
#' LFMM/EMMAX plots are grid-level only (one histogram grid / QQ grid /
#' lambda-vs-K-or-#PCs / hits-vs-K-or-#PCs PNG per engine, never one file per
#' rung) — the ladder tables carry the per-value numbers. RDA IS per-model:
#' one biplot/screeplot/axis-anova set per swept Condition()-PC value, picked
#' with the selector below the model-comparison plot.
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

            bslib::nav_panel("LFMM K ladder",
                mod_image_card_ui(ns("lfmm_screeplot")),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("lfmm_hist")),
                    mod_image_card_ui(ns("lfmm_qq")),
                    mod_image_card_ui(ns("lfmm_lambda")),
                    mod_image_card_ui(ns("lfmm_hits"))
                ),
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mt-3",
                    htmltools::h6("Ladder Table", class = "mb-0"),
                    shiny::uiOutput(ns("lfmm_table_note"), inline = TRUE)
                ),
                DT::DTOutput(ns("lfmm_table"))
            ),

            bslib::nav_panel("EMMAX / GAPIT #PC ladder",
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("emmax_hist")),
                    mod_image_card_ui(ns("emmax_qq")),
                    mod_image_card_ui(ns("emmax_lambda")),
                    mod_image_card_ui(ns("emmax_hits"))
                ),
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mt-3",
                    htmltools::h6("Ladder Table", class = "mb-0"),
                    shiny::uiOutput(ns("emmax_table_note"), inline = TRUE)
                ),
                DT::DTOutput(ns("emmax_table"))
            ),

            bslib::nav_panel("RDA setup",
                mod_image_card_ui(ns("rda_comparison")),
                htmltools::div(
                    class = "control-bar mt-3",
                    shiny::uiOutput(ns("rda_model_selector"))
                ),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("rda_biplot")),
                    mod_image_card_ui(ns("rda_screeplot"))
                ),
                htmltools::h6("Tables", class = "mt-3"),
                bslib::accordion(
                    open = FALSE, multiple = TRUE,
                    bslib::accordion_panel("Selected model — axis anova",
                        DT::DTOutput(ns("rda_model_axis_table"))),
                    bslib::accordion_panel("Condition-PC comparison (all models)",
                        DT::DTOutput(ns("rda_ladder_table"))),
                    bslib::accordion_panel("Predictor collinearity",
                        DT::DTOutput(ns("rda_collinearity_table"))),
                    bslib::accordion_panel("ordiR2step selection path",
                        DT::DTOutput(ns("rda_ordir2step_table")))
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
        render_tsv_table <- function(path, columns = NULL, colnames = NULL) {
            dt <- load_tsv_dt(path)
            if (!is.null(columns) && nrow(dt) > 0) {
                keep <- intersect(columns, names(dt))
                dt <- dt[, ..keep]
                if (!is.null(colnames) && length(colnames) == length(columns)) {
                    # Rename by matching against the ORIGINAL columns list
                    # (not `keep`) so a column missing from this particular
                    # TSV doesn't shift every later header off by one.
                    names(dt) <- colnames[match(names(dt), columns)]
                }
            }
            safe_datatable(as.data.frame(dt),
                           options = list(dom = "t", scrollX = TRUE, pageLength = -1),
                           rownames = FALSE)
        }

        # ── Recommendations ─────────────────────────────────────────────────
        output$recommendations_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "", "pregea_recommendations"))
        })

        # ── LFMM K ladder ────────────────────────────────────────────────────
        # Ladder table trimmed to the 6 columns that actually decide K — see
        # help_note("pregea_ladder_table") for how to read them together.
        # Hidden (still on disk in the full TSV): engine, rung_param (constant
        # per table), frac_p_gt_half (redundant w/ flatness), hist_spike0
        # (already encoded in hist_shape), hits_bonf (never plotted), n_tests,
        # rung_value_num (dcast artifact).
        LADDER_COLS <- c("trait", "rung_value", "hist_shape", "hist_flatness_ks", "hits_qval", "lambda_gc")
        LADDER_COLNAMES <- c("Trait", "Value", "Histogram shape", "Flatness (KS)", "Hits (FDR 0.1)", "lambda_GC")

        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Run PreGEA.")

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
        output$lfmm_table_note <- shiny::renderUI(help_note("pregea_ladder_table"))
        output$lfmm_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "lfmm", "lfmm_ladder"),
                             columns = LADDER_COLS, colnames = LADDER_COLNAMES)
        })

        # ── EMMAX / GAPIT #PC ladder ─────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Run PreGEA.")

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
        output$emmax_table_note <- shiny::renderUI(help_note("pregea_ladder_table"))
        output$emmax_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(pregea_table_path(pd$name, "emmax", "emmax_ladder"),
                             columns = LADDER_COLS, colnames = LADDER_COLNAMES)
        })

        # ── RDA setup ────────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            mod_image_card_server("rda_comparison",
                path       = shiny::reactive(pregea_plot_path(pd$name, "rda", "rda_model_comparison")),
                title      = shiny::reactive("Model Comparison (all Condition()-PC values)"),
                dl_name    = shiny::reactive("rda_model_comparison"),
                suggestion = shiny::reactive("Run PreGEA."),
                note       = shiny::reactive(help_note("pregea_rda_model_comparison"))
            )
        })

        # Model selector — one biplot/screeplot/axis-anova set per swept
        # Condition()-PC value (pregea_rda_setup.R writes them all in one
        # rule's internal loop, no Snakemake fan-out — see that rule's header).
        # Defaults to the recommended condition_pcs when available.
        rda_selected_pc <- shiny::reactiveVal(NULL)
        output$rda_model_selector <- shiny::renderUI({
            pd <- project_data()
            models <- find_rda_models(pd$name)
            if (length(models) == 0) {
                return(htmltools::p(class = "text-muted small mb-0",
                    "No RDA models yet — run PreGEA."))
            }
            recs <- pregea_recommendations_lookup(pd$name)
            recommended <- suppressWarnings(as.integer(recs$RDA$condition_pcs$value))
            default_sel <- if (!is.null(rda_selected_pc()) && rda_selected_pc() %in% models) {
                rda_selected_pc()
            } else if (!is.na(recommended) && recommended %in% models) {
                recommended
            } else {
                models[1]
            }
            choices <- setNames(models, paste0("Condition PCs = ", models,
                                               ifelse(models == recommended, " (recommended)", "")))
            shiny::selectInput(ns("rda_model_pc"), "RDA model", choices = choices, selected = default_sel, width = "300px")
        })
        shiny::observeEvent(input$rda_model_pc, rda_selected_pc(input$rda_model_pc), ignoreInit = TRUE)

        shiny::observe({
            pd <- project_data()
            pc <- input$rda_model_pc
            shiny::req(pc)
            mod_image_card_server("rda_biplot",
                path       = shiny::reactive(rda_model_biplot_path(pd$name, pc)),
                title      = shiny::reactive(paste0("Site-scores Biplot (Condition PCs = ", pc, ")")),
                dl_name    = shiny::reactive(paste0("rda_biplot_pc", pc)),
                suggestion = shiny::reactive("Run PreGEA."),
                note       = shiny::reactive(help_note("pregea_rda_model_biplot"))
            )
            mod_image_card_server("rda_screeplot",
                path       = shiny::reactive(rda_model_screeplot_path(pd$name, pc)),
                title      = shiny::reactive(paste0("Axis Screeplot (Condition PCs = ", pc, ")")),
                dl_name    = shiny::reactive(paste0("rda_screeplot_pc", pc)),
                suggestion = shiny::reactive("Run PreGEA."),
                note       = shiny::reactive(help_note("pregea_rda_model_screeplot"))
            )
            output$rda_model_axis_table <- DT::renderDataTable({
                render_tsv_table(rda_model_axis_anova_path(pd$name, pc))
            })
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
