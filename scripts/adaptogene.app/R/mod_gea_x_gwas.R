#' GEAxGWAS tab UI
#'
#' Two independent filter bars (GEA above Miami, GWAS below) drive independent
#' region sets on each half of the Miami plot. Overlap regions — where a GEA region
#' and a GWAS region intersect — are highlighted in orange across both halves and
#' feed the region explorer (genes / SNPs / GO / regionplot / haplotype).
#'
#' @param id module namespace id
#' @noRd
mod_gea_x_gwas_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Notification about filter bars and fill buttons
        htmltools::div(
            class = "alert alert-info py-2 mb-2 small",
            bsicons::bs_icon("info-circle"), " ",
            "The GEA and GWAS filter bars below drive the Miami plot, region overlap explorer, ",
            "and Pairwise Trait Overlap table. ",
            "Use the ", htmltools::strong("Fill from …"), " buttons to copy your settings from the GEA or GWAS tabs."
        ),

        # GEA filter bar
        bslib::card(
            class = "mb-0",
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                htmltools::div(
                    htmltools::span("GEA Filter",
                        class = "fw-semibold text-primary me-2"),
                    htmltools::span("(above zero line)",
                        class = "text-muted small")
                ),
                shiny::actionButton(ns("fill_gea"), "Fill from GEA tab",
                                    class = "btn-sm btn-outline-primary")
            ),
            bslib::card_body(
            class = "p-2",
            shiny::uiOutput(ns("gea_threshold_bar")),
            shiny::uiOutput(ns("gea_filter_bar"))
        )
        ),

        # Miami Manhattan — GEA filter above, GWAS filter below
        mod_manhattan_overlay_ui(
            ns("miami"),
            height    = "450px",
            filter_ui = NULL
        ),

        # GWAS filter bar
        bslib::card(
            class = "mb-0 mt-1",
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                htmltools::div(
                    htmltools::span("GWAS Filter",
                        class = "fw-semibold text-success me-2"),
                    htmltools::span("(below zero line)",
                        class = "text-muted small")
                ),
                shiny::actionButton(ns("fill_gwas"), "Fill from GWAS tab",
                                    class = "btn-sm btn-outline-success")
            ),
            bslib::card_body(
            class = "p-2",
            shiny::uiOutput(ns("gwas_threshold_bar")),
            shiny::uiOutput(ns("gwas_filter_bar"))
        )
        ),

        # Overlap bounds selector (persisted per project)
        bslib::card(
            class = "mt-2 mb-3",
            bslib::card_body(
                class = "py-2 px-3",
                htmltools::div(
                    class = "d-flex align-items-center gap-4",
                    htmltools::span("Overlap region bounds:", class = "fw-semibold me-1"),
                    bslib::tooltip(
                        bsicons::bs_icon("info-circle", class = "text-muted", size = "0.95em"),
                        htmltools::tags$div(
                            htmltools::tags$strong("Union"), ": merged extent of GEA \u222a GWAS region \u2014 covers the full area of either source.",
                            htmltools::tags$br(),
                            htmltools::tags$strong("Intersection"), ": only the overlapping portion between GEA and GWAS regions.",
                            htmltools::tags$br(),
                            htmltools::tags$strong("GEA only"), ": GEA region boundaries, ignoring GWAS extent.",
                            htmltools::tags$br(),
                            htmltools::tags$strong("GWAS only"), ": GWAS region boundaries, ignoring GEA extent."
                        ),
                        placement = "right"
                    ),
                    shiny::radioButtons(
                        ns("overlap_bounds"), label = NULL,
                        choices = c(
                            "Union (GEA \u222a GWAS)"         = "union",
                            "Intersection (GEA \u2229 GWAS)"  = "intersection",
                            "GEA region only"                  = "gea_only",
                            "GWAS region only"                 = "gwas_only"
                        ),
                        selected = "union",
                        inline   = TRUE
                    )
                )
            )
        ),

        # Overlap region explorer
        mod_region_explorer_ui(ns("region_explorer")),

        # Pairwise Trait Overlap — collapsed by default (independent analysis)
        bslib::accordion(
            id       = ns("pairwise_accordion"),
            open     = FALSE,
            multiple = FALSE,
            class    = "mt-4",
            bslib::accordion_panel(
                "Pairwise Trait Overlap",
                value = "pairwise",
                icon  = bsicons::bs_icon("intersect"),
                mod_pairwise_overlap_ui(ns("pairwise"))
            )
        )
    )
}

#' Overlapping Regions tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_gea_x_gwas_server <- function(id, project_data, run_trigger = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        module <- MOD_GEAXGWAS

        # ── GEA side data ──────────────────────────────────────────────────────
        gea_methods <- shiny::reactive(find_assoc_methods(project_data()$name, MOD_GEA))

        gea_method_pvalues <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_pvalues_cached(pd$name, MOD_GEA, pd$k_best)
        })
        gea_method_wza_pvalues <- shiny::reactive({
            if (!gea_regime_wza()) return(list())
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_wza_pvalues_cached(pd$name, MOD_GEA, pd$k_best)
        })
        gea_effective_pvalues <- shiny::reactive({
            if (gea_regime_wza()) gea_method_wza_pvalues() else gea_method_pvalues()
        })

        gea_all_trait_names <- shiny::reactive({
            pv <- gea_method_pvalues()
            if (length(pv) == 0) return(character(0))
            fixed <- c("SNPID","chr","pos","n_snps","mean_maf")
            sort(unique(unlist(lapply(pv, function(dt) setdiff(names(dt), fixed)))))
        })

        gea_effective_sigsnps <- shiny::reactive({
            pd <- project_data()
            compute_method_sigsnps_cached(
                pvalues_list = gea_effective_pvalues(),
                type         = gea_threshold_type(),
                value        = gea_threshold_value(),
                k            = pd$k_best,
                regime       = if (gea_regime_wza()) "wza" else "snp",
                project      = pd$name,
                module       = MOD_GEA,
                cutoffs      = gea_combo_thresholds()
            )
        })

        gea_traits <- shiny::reactive({
            sig <- gea_effective_sigsnps()
            if (length(sig) == 0) return(character(0))
            all_dt <- data.table::rbindlist(sig, use.names = TRUE, fill = TRUE)
            if (nrow(all_dt) == 0) return(character(0))
            sort(unique(all_dt$trait))
        })

        gea_trait_colors  <- shiny::reactive(trait_color_map(gea_all_trait_names()))
        gea_method_shapes <- shiny::reactive(method_shape_map(gea_methods()))

        gea_combo_counts <- shiny::reactive({
            snps <- gea_effective_sigsnps()
            counts <- list()
            for (m in names(snps)) {
                dt <- snps[[m]]
                if (nrow(dt) > 0) {
                    tc <- dt[, .(n = .N), by = "trait"]
                    for (i in seq_len(nrow(tc)))
                        counts[[paste0(tc$trait[i], "::", m)]] <- tc$n[i]
                }
            }
            counts
        })

        # ── GWAS side data ─────────────────────────────────────────────────────
        gwas_methods <- shiny::reactive(find_assoc_methods(project_data()$name, MOD_GWAS))

        gwas_method_pvalues <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_pvalues_cached(pd$name, MOD_GWAS, pd$k_best)
        })
        gwas_method_wza_pvalues <- shiny::reactive({
            if (!gwas_regime_wza()) return(list())
            if (!is.null(run_trigger)) run_trigger()  # invalidate when pipeline completes
            pd <- project_data()
            load_all_method_wza_pvalues_cached(pd$name, MOD_GWAS, pd$k_best)
        })
        gwas_effective_pvalues <- shiny::reactive({
            if (gwas_regime_wza()) gwas_method_wza_pvalues() else gwas_method_pvalues()
        })

        gwas_all_trait_names <- shiny::reactive({
            pv <- gwas_method_pvalues()
            if (length(pv) == 0) return(character(0))
            fixed <- c("SNPID","chr","pos","n_snps","mean_maf")
            sort(unique(unlist(lapply(pv, function(dt) setdiff(names(dt), fixed)))))
        })

        gwas_effective_sigsnps <- shiny::reactive({
            pd <- project_data()
            compute_method_sigsnps_cached(
                pvalues_list = gwas_effective_pvalues(),
                type         = gwas_threshold_type(),
                value        = gwas_threshold_value(),
                k            = pd$k_best,
                regime       = if (gwas_regime_wza()) "wza" else "snp",
                project      = pd$name,
                module       = MOD_GWAS,
                cutoffs      = gwas_combo_thresholds()
            )
        })

        gwas_traits <- shiny::reactive({
            sig <- gwas_effective_sigsnps()
            if (length(sig) == 0) return(character(0))
            all_dt <- data.table::rbindlist(sig, use.names = TRUE, fill = TRUE)
            if (nrow(all_dt) == 0) return(character(0))
            sort(unique(all_dt$trait))
        })

        gwas_trait_colors  <- shiny::reactive(trait_color_map(gwas_all_trait_names()))
        gwas_method_shapes <- shiny::reactive(method_shape_map(gwas_methods()))

        gwas_combo_counts <- shiny::reactive({
            snps <- gwas_effective_sigsnps()
            counts <- list()
            for (m in names(snps)) {
                dt <- snps[[m]]
                if (nrow(dt) > 0) {
                    tc <- dt[, .(n = .N), by = "trait"]
                    for (i in seq_len(nrow(tc)))
                        counts[[paste0(tc$trait[i], "::", m)]] <- tc$n[i]
                }
            }
            counts
        })

        # File fingerprint reactives — change when the pipeline regenerates pvalue TSVs.
        gea_pvalue_fingerprint <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()
            pd <- project_data()
            pvalues_file_fingerprint(pd$name, MOD_GEA, pd$k_best,
                                     regime = if (gea_regime_wza()) "wza" else "snp")
        })
        gwas_pvalue_fingerprint <- shiny::reactive({
            if (!is.null(run_trigger)) run_trigger()
            pd <- project_data()
            pvalues_file_fingerprint(pd$name, MOD_GWAS, pd$k_best,
                                     regime = if (gwas_regime_wza()) "wza" else "snp")
        })

        # Per-cell significance thresholds for filter bar cells (bindCached per side)
        gea_combo_thresholds <- shiny::reactive({
            compute_method_thresholds(
                pvalues_list = gea_effective_pvalues(),
                type         = gea_threshold_type(),
                value        = gea_threshold_value()
            )
        }) |> shiny::bindCache(gea_threshold_type(), gea_threshold_value(),
                               gea_regime_wza(), gea_pvalue_fingerprint())

        gwas_combo_thresholds <- shiny::reactive({
            compute_method_thresholds(
                pvalues_list = gwas_effective_pvalues(),
                type         = gwas_threshold_type(),
                value        = gwas_threshold_value()
            )
        }) |> shiny::bindCache(gwas_threshold_type(), gwas_threshold_value(),
                               gwas_regime_wza(), gwas_pvalue_fingerprint())

        # ── Initialize per-side filter params from config/region_params ─────────
        .init_dist_param <- function(input_id, rp_key, config_default_fn) {
            shiny::observe({
                pd  <- project_data()
                rp  <- read_region_params(pd$name)
                val <- get_global_param(rp, module, rp_key)
                d <- if (!is.null(val)) as.integer(val)
                     else config_default_fn(pd$config, pd$name)
                shiny::updateNumericInput(session, input_id, value = d)
            })
            shiny::observeEvent(input[[input_id]], {
                v  <- input[[input_id]]
                pd <- project_data()
                if (is.null(v) || is.na(v) || is.null(pd)) return()
                rp <- read_region_params(pd$name)
                rp <- set_global_param(rp, module, rp_key, as.integer(v))
                save_region_params(pd$name, rp)
            }, ignoreInit = TRUE)
        }

        .init_dist_param("gea_snp_clumping_distance", "gea_snp_clumping_distance",
            function(cfg, proj) resolve_ui_snp_clumping_distance(cfg, "GEA", proj))
        .init_dist_param("gwas_snp_clumping_distance", "gwas_snp_clumping_distance",
            function(cfg, proj) resolve_ui_snp_clumping_distance(cfg, "GWAS", proj))

        gea_snp_clumping_distance <- shiny::reactive({
            v <- input$gea_snp_clumping_distance
            if (is.null(v) || is.na(v) || v < 1000L) 1000000L else as.integer(v)
        })
        gwas_snp_clumping_distance <- shiny::reactive({
            v <- input$gwas_snp_clumping_distance
            if (is.null(v) || is.na(v) || v < 1000L) 1000000L else as.integer(v)
        })

        # ── Regime (independent per side) ─────────────────────────────────────
        .init_regime_param <- function(input_id, rp_key) {
            shiny::observe({
                pd    <- project_data()
                rp    <- read_region_params(pd$name)
                saved <- get_global_param(rp, module, rp_key)
                if (!is.null(saved))
                    bslib::update_switch(input_id, value = isTRUE(saved), session = session)
            })
            shiny::observeEvent(input[[input_id]], {
                pd <- project_data()
                if (is.null(pd)) return()
                rp <- read_region_params(pd$name)
                rp <- set_global_param(rp, module, rp_key, isTRUE(input[[input_id]]))
                save_region_params(pd$name, rp)
            }, ignoreInit = TRUE)
        }

        .init_regime_param("gea_regime", "gea_regime")
        .init_regime_param("gwas_regime", "gwas_regime")

        gea_regime_wza  <- shiny::reactive(isTRUE(input$gea_regime))
        gwas_regime_wza <- shiny::reactive(isTRUE(input$gwas_regime))

        # ── Threshold params per side ──────────────────────────────────────────
        .init_threshold_param <- function(type_id, value_id, rp_type_key, rp_value_key,
                                          config_src_module) {
            shiny::observe({
                pd  <- project_data(); cfg <- pd$config; rp <- read_region_params(pd$name)
                saved_t <- get_global_param(rp, module, rp_type_key)
                saved_v <- get_global_param(rp, module, rp_value_key)
                if (!is.null(saved_t)) {
                    shiny::updateSelectInput(session, type_id, selected = saved_t)
                } else {
                    def <- default_threshold(cfg, config_src_module)
                    shiny::updateSelectInput(session, type_id, selected = def$type)
                }
                if (!is.null(saved_v)) {
                    shiny::updateNumericInput(session, value_id, value = as.numeric(saved_v))
                } else if (is.null(saved_t)) {
                    def <- default_threshold(cfg, config_src_module)
                    shiny::updateNumericInput(session, value_id, value = def$value)
                }
            })
            shiny::observeEvent(input[[type_id]], {
                pd <- project_data(); if (is.null(pd)) return()
                rp <- read_region_params(pd$name)
                rp <- set_global_param(rp, module, rp_type_key, input[[type_id]])
                save_region_params(pd$name, rp)
            }, ignoreInit = TRUE)
            shiny::observeEvent(input[[value_id]], {
                v <- input[[value_id]]; pd <- project_data()
                if (is.null(pd) || is.null(v) || is.na(v)) return()
                rp <- read_region_params(pd$name)
                rp <- set_global_param(rp, module, rp_value_key, as.numeric(v))
                save_region_params(pd$name, rp)
            }, ignoreInit = TRUE)
        }

        .init_threshold_param("gea_threshold_type", "gea_threshold_value",
            "gea_threshold_type", "gea_threshold_value", MOD_GEA)
        .init_threshold_param("gwas_threshold_type", "gwas_threshold_value",
            "gwas_threshold_type", "gwas_threshold_value", MOD_GWAS)

        # ── Strategy persistence per side ─────────────────────────────────────
        # Mirrors .init_regime_param; strategy is the radioButton input not a switch.
        .init_strategy_param <- function(input_id, rp_key) {
            shiny::observe({
                pd    <- project_data()
                rp    <- read_region_params(pd$name)
                saved <- get_global_param(rp, module, rp_key)
                if (!is.null(saved))
                    shiny::updateRadioButtons(session, input_id,
                                              selected = .normalize_strategy(saved))
            })
            shiny::observeEvent(input[[input_id]], {
                pd <- project_data(); if (is.null(pd)) return()
                rp <- read_region_params(pd$name)
                rp <- set_global_param(rp, module, rp_key, input[[input_id]])
                save_region_params(pd$name, rp)
            }, ignoreInit = TRUE)
        }

        .init_strategy_param("gea_combine_strategy", "gea_strategy")
        .init_strategy_param("gwas_combine_strategy", "gwas_strategy")

        # ── Fill from GEA / GWAS tab buttons ─────────────────────────────────
        # Reads the currently persisted params from the source module's region_params.json
        # keys and pushes them into the prefixed GEAxGWAS filter-bar inputs.
        # The existing .init_* observeEvent handlers auto-persist the copied values.
        .fill_from <- function(src_module, prefix) {
            pd  <- project_data()
            rp  <- read_region_params(pd$name)
            cfg <- pd$config
            src_str <- if (src_module == MOD_GEA) "GEA" else "GWAS"

            tt <- get_global_param(rp, src_module, "threshold_type") %||%
                      default_threshold(cfg, src_module)$type
            tv <- as.numeric(get_global_param(rp, src_module, "threshold_value") %||%
                      default_threshold(cfg, src_module)$value)
            cd <- as.integer(get_global_param(rp, src_module, "snp_clumping_distance") %||%
                      resolve_ui_snp_clumping_distance(cfg, src_str, pd$name))
            rg <- isTRUE(get_global_param(rp, src_module, "regime"))
            st <- get_global_param(rp, src_module, "combine_strategy")

            shiny::updateSelectInput(session,  paste0(prefix, "threshold_type"),  selected = tt)
            shiny::updateNumericInput(session, paste0(prefix, "threshold_value"), value = tv)
            shiny::updateNumericInput(session, paste0(prefix, "snp_clumping_distance"), value = cd)
            bslib::update_switch(paste0(prefix, "regime"), value = rg, session = session)
            if (!is.null(st))
                shiny::updateRadioButtons(session, paste0(prefix, "combine_strategy"),
                                          selected = .normalize_strategy(st))
            shiny::showNotification(
                sprintf("Copied %s tab settings into the GEAxGWAS %s filter bar.", src_module, toupper(prefix)),
                type = "message", duration = 3
            )
        }

        shiny::observeEvent(input$fill_gea,  .fill_from(MOD_GEA,  "gea_"))
        shiny::observeEvent(input$fill_gwas, .fill_from(MOD_GWAS, "gwas_"))

        gea_threshold_type_raw  <- shiny::reactive(input$gea_threshold_type  %||% "bonf")
        gwas_threshold_type_raw <- shiny::reactive(input$gwas_threshold_type %||% "bonf")
        # Debounce type so switching to "qval" doesn't fire the full 2M-row scan instantly
        gea_threshold_type  <- shiny::debounce(gea_threshold_type_raw,  500)
        gwas_threshold_type <- shiny::debounce(gwas_threshold_type_raw, 500)

        .make_threshold_value <- function(input_id, src_module) {
            raw <- shiny::reactive({
                v <- input[[input_id]]
                if (is.null(v) || is.na(v) || v <= 0)
                    default_threshold(project_data()$config, src_module)$value
                else as.numeric(v)
            })
            shiny::debounce(raw, 500)
        }
        gea_threshold_value  <- .make_threshold_value("gea_threshold_value",  MOD_GEA)
        gwas_threshold_value <- .make_threshold_value("gwas_threshold_value", MOD_GWAS)

        output$gea_threshold_hint <- shiny::renderUI({
            t <- gea_threshold_type()
            hint <- switch(t, bonf="α for Bonferroni", qval="FDR q-value (0–1)",
                           top="Top N SNPs", custom="Raw p-value (e.g. 1e-5)", "")
            htmltools::span(class="text-muted small mt-1", hint)
        })
        output$gwas_threshold_hint <- shiny::renderUI({
            t <- gwas_threshold_type()
            hint <- switch(t, bonf="α for Bonferroni", qval="FDR q-value (0–1)",
                           top="Top N SNPs", custom="Raw p-value (e.g. 1e-5)", "")
            htmltools::span(class="text-muted small mt-1", hint)
        })

        # Coerce threshold values when type changes
        shiny::observeEvent(input$gea_threshold_type, {
            v <- input$gea_threshold_value
            if (!threshold_value_valid_for_type(input$gea_threshold_type, v))
                shiny::updateNumericInput(session, "gea_threshold_value",
                    value = threshold_value_default_for_type(input$gea_threshold_type))
        }, ignoreInit = TRUE)
        shiny::observeEvent(input$gwas_threshold_type, {
            v <- input$gwas_threshold_value
            if (!threshold_value_valid_for_type(input$gwas_threshold_type, v))
                shiny::updateNumericInput(session, "gwas_threshold_value",
                    value = threshold_value_default_for_type(input$gwas_threshold_type))
        }, ignoreInit = TRUE)

        # ── Always-present threshold bars ─────────────────────────────────────
        output$gea_threshold_bar <- shiny::renderUI({
            pd <- project_data()
            shiny::isolate({
                build_threshold_bar_ui(
                    ns                    = ns,
                    input_prefix          = "gea_",
                    regime_value          = isTRUE(input$gea_regime),
                    threshold_type_value  = input$gea_threshold_type  %||% "bonf",
                    threshold_value_value = input$gea_threshold_value %||%
                        default_threshold(pd$config, MOD_GEA)$value,
                    regime_context        = "overlap"
                )
            })
        })
        output$gwas_threshold_bar <- shiny::renderUI({
            pd <- project_data()
            shiny::isolate({
                build_threshold_bar_ui(
                    ns                    = ns,
                    input_prefix          = "gwas_",
                    regime_value          = isTRUE(input$gwas_regime),
                    threshold_type_value  = input$gwas_threshold_type  %||% "bonf",
                    threshold_value_value = input$gwas_threshold_value %||%
                        default_threshold(pd$config, MOD_GWAS)$value,
                    regime_context        = "overlap"
                )
            })
        })

        # ── Overlap bounds (persisted) ─────────────────────────────────────────
        shiny::observe({
            pd  <- project_data()
            rp  <- read_region_params(pd$name)
            val <- get_global_param(rp, module, "overlap_bounds")
            if (!is.null(val)) shiny::updateRadioButtons(session, "overlap_bounds", selected = val)
        })
        shiny::observeEvent(input$overlap_bounds, {
            pd <- project_data()
            if (is.null(pd)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, module, "overlap_bounds", input$overlap_bounds)
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        overlap_bounds <- shiny::reactive({
            input$overlap_bounds %||% "union"
        })

        # ── GEA filter bar UI (matrix + strategy + clumping) ──────────────────
        output$gea_filter_bar <- shiny::renderUI({
            pd <- project_data()
            shiny::isolate({
                rp           <- read_region_params(pd$name)
                saved_strat  <- get_global_param(rp, module, "gea_strategy")
                strat_val    <- if (!is.null(saved_strat)) .normalize_strategy(saved_strat) else "Union"
            })
            build_filter_bar_ui(
                ns                          = ns,
                traits                      = gea_all_trait_names(),
                methods                     = gea_methods(),
                trait_colors                = gea_trait_colors(),
                combo_counts                = gea_combo_counts(),
                combo_thresholds            = gea_combo_thresholds(),
                default_strategy_value      = strat_val,
                snp_clumping_distance_value = gea_snp_clumping_distance(),
                input_prefix                = "gea_"
            )
        })

        # ── GWAS filter bar UI (matrix + strategy + clumping) ─────────────────
        output$gwas_filter_bar <- shiny::renderUI({
            pd <- project_data()
            shiny::isolate({
                rp           <- read_region_params(pd$name)
                saved_strat  <- get_global_param(rp, module, "gwas_strategy")
                strat_val    <- if (!is.null(saved_strat)) .normalize_strategy(saved_strat) else "Union"
            })
            build_filter_bar_ui(
                ns                          = ns,
                traits                      = gwas_all_trait_names(),
                methods                     = gwas_methods(),
                trait_colors                = gwas_trait_colors(),
                combo_counts                = gwas_combo_counts(),
                combo_thresholds            = gwas_combo_thresholds(),
                default_strategy_value      = strat_val,
                snp_clumping_distance_value = gwas_snp_clumping_distance(),
                input_prefix                = "gwas_"
            )
        })

        # ── Interactive sig SNPs (independent per side) ────────────────────────
        gea_interactive_sigsnps <- shiny::reactive({
            compute_interactive_sigsnps(
                all_method_sigsnps  = gea_effective_sigsnps(),
                tm_selection_json   = input$gea_tm_selection,
                combo_counts        = gea_combo_counts(),
                known_traits        = gea_all_trait_names(),
                strategy            = input$gea_combine_strategy %||% "Union",
                clumping_distance   = gea_snp_clumping_distance(),
                project_name        = project_data()$name,
                module              = MOD_GEA
            )
        })

        gwas_interactive_sigsnps <- shiny::reactive({
            compute_interactive_sigsnps(
                all_method_sigsnps  = gwas_effective_sigsnps(),
                tm_selection_json   = input$gwas_tm_selection,
                combo_counts        = gwas_combo_counts(),
                known_traits        = gwas_all_trait_names(),
                strategy            = input$gwas_combine_strategy %||% "Union",
                clumping_distance   = gwas_snp_clumping_distance(),
                project_name        = project_data()$name,
                module              = MOD_GWAS
            )
        })

        # ── Independent regions (one set per source) ───────────────────────────
        gea_regions <- shiny::reactive({
            snps <- gea_interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(.empty_regions())
            compute_all_regions(snps, gea_snp_clumping_distance())
        })

        gwas_regions <- shiny::reactive({
            snps <- gwas_interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(.empty_regions())
            compute_all_regions(snps, gwas_snp_clumping_distance())
        })

        # ── Overlap detection ──────────────────────────────────────────────────
        overlap_pairs <- shiny::reactive({
            compute_region_overlaps(gea_regions(), gwas_regions())
        })

        # ── Overlap regions for the explorer (with bounds strategy applied) ────
        overlap_regions <- shiny::reactive({
            pairs <- overlap_pairs()
            if (is.null(pairs) || nrow(pairs) == 0) return(.empty_overlap_regions())
            compute_all_overlap_regions(
                overlap_pairs  = pairs,
                gea_regions    = gea_regions(),
                gwas_regions   = gwas_regions(),
                strategy       = overlap_bounds(),
                gea_sig_snps   = gea_interactive_sigsnps(),
                gwas_sig_snps  = gwas_interactive_sigsnps()
            )
        })

        # ── Unified sig SNPs for the detail panel (region_id = overlap region_id) ──
        unified_sigsnps <- shiny::reactive({
            gea_snps  <- gea_interactive_sigsnps()
            gwas_snps <- gwas_interactive_sigsnps()

            parts <- list()
            if (!is.null(gea_snps) && nrow(gea_snps) > 0) {
                gea_copy <- data.table::copy(gea_snps)
                gea_copy[, source := "gea"]
                parts[["gea"]] <- gea_copy
            }
            if (!is.null(gwas_snps) && nrow(gwas_snps) > 0) {
                gwas_copy <- data.table::copy(gwas_snps)
                gwas_copy[, source := "gwas"]
                parts[["gwas"]] <- gwas_copy
            }
            if (length(parts) == 0) return(NULL)

            unified <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)

            # Assign overlap region_ids so Miami plot click → correct region
            ov_regs <- overlap_regions()
            if (!is.null(ov_regs) && nrow(ov_regs) > 0) {
                unified <- assign_overlap_region_ids(unified, ov_regs)
            }
            unified
        })

        # ── Region explorer ────────────────────────────────────────────────────
        explorer <- mod_region_explorer_server("region_explorer",
            project_data        = project_data,
            module              = module,
            interactive_sigsnps = unified_sigsnps,
            region_distance     = shiny::reactive(1000000L),  # not used when override_regions supplied
            override_regions    = overlap_regions,
            regime              = shiny::reactive(gea_regime_wza() && gwas_regime_wza())
        )

        # ── Miami plot ─────────────────────────────────────────────────────────
        # Combine trait/method maps for both sources (unified legend, stable colors)
        all_traits  <- shiny::reactive(sort(union(gea_all_trait_names(), gwas_all_trait_names())))
        all_methods <- shiny::reactive(union(gea_methods(), gwas_methods()))
        all_trait_colors  <- shiny::reactive(trait_color_map(all_traits()))
        all_method_shapes <- shiny::reactive(method_shape_map(all_methods()))

        # Miami threshold y: min threshold across both sides (most stringent for display)
        # Derived from bindCached combo_thresholds — no redundant full-vector scan.
        miami_threshold_y <- shiny::reactive({
            all_ct <- c(unlist(gea_combo_thresholds()), unlist(gwas_combo_thresholds()))
            all_ct <- all_ct[!is.na(all_ct) & all_ct > 0]
            if (length(all_ct) == 0) return(NULL)
            -log10(min(all_ct))
        })

        miami_wza_bg <- shiny::reactive({
            if (!(gea_regime_wza() && gwas_regime_wza())) return(NULL)
            pd <- project_data(); k <- pd$k_best
            if (is.na(k)) return(NULL)
            miami_wza_bg_path(pd$name, k)
        })
        miami_wza_co <- shiny::reactive({
            if (!(gea_regime_wza() && gwas_regime_wza())) return(NULL)
            pd <- project_data(); k <- pd$k_best
            if (is.na(k)) return(NULL)
            miami_wza_coords_path(pd$name, k)
        })

        miami_click <- mod_manhattan_overlay_server("miami",
            project_data         = project_data,
            module               = module,
            is_miami             = TRUE,
            combined             = TRUE,
            title_label          = shiny::reactive({
                if (gea_regime_wza() && gwas_regime_wza())
                    "Miami Plot WZA (GEA \u2191 | GWAS \u2193)"
                else
                    "Miami Plot (GEA \u2191 | GWAS \u2193)"
            }),
            note                 = shiny::reactive(help_note("miami_plot")),
            regions              = gea_regions,
            current_region_id    = explorer$selected_region_id,
            show_regions_control = FALSE,
            sig_snps_override    = unified_sigsnps,
            trait_colors         = all_trait_colors,
            method_shapes        = all_method_shapes,
            gwas_regions         = gwas_regions,
            overlap_regions      = overlap_regions,
            bg_path_override     = miami_wza_bg,
            coords_path_override = miami_wza_co,
            threshold_y          = miami_threshold_y
        )

        # Miami SNP click → select the enclosing overlap region
        shiny::observeEvent(miami_click(), {
            rid <- miami_click()
            if (!is.null(rid) && nzchar(rid)) {
                # Only select if the clicked region_id is an overlap region
                ov_regs <- overlap_regions()
                if (!is.null(ov_regs) && rid %in% ov_regs$region_id) {
                    explorer$selected_region_id(rid)
                }
                # Clicks on GEA/GWAS-only region IDs are intentionally ignored
                # (they're contextual highlights, not selectable overlap regions)
            }
        }, ignoreNULL = TRUE)

        # ── Pairwise Trait Overlap (on-demand async, driven by filter-bar params) ─
        miami_coords <- shiny::reactive({
            pd   <- project_data()
            k    <- pd$k_best
            if (is.na(k)) return(NULL)
            cp   <- miami_coords_path(pd$name, k)
            fp   <- if (file.exists(cp)) as.character(file.info(cp)$mtime) else "missing"
            load_cached(paste0("coords_miami_", pd$name, "_", k), function() {
                load_coords(cp)
            }, fingerprint = fp)
        })

        # Parameter bundle for pairwise hash — includes a file fingerprint so
        # a pipeline re-run (new pvalue files) invalidates the cache automatically.
        pairwise_params <- shiny::reactive({
            pd <- project_data()
            k  <- pd$k_best
            fp <- pairwise_file_fingerprint(pd$name)
            list(
                k              = k,
                gea_thr_type   = gea_threshold_type(),
                gea_thr_value  = gea_threshold_value(),
                gea_regime     = isTRUE(input$gea_regime),
                gwas_thr_type  = gwas_threshold_type(),
                gwas_thr_value = gwas_threshold_value(),
                gwas_regime    = isTRUE(input$gwas_regime),
                gea_window     = gea_snp_clumping_distance(),
                gwas_window    = gwas_snp_clumping_distance(),
                min_snps       = as.integer(
                    tryCatch(pd$config$GEAxGWAS$pairwise$min_snps, error = function(e) 2L) %||% 2L
                ),
                fp             = fp
            )
        })

        mod_pairwise_overlap_server("pairwise",
            project_data    = project_data,
            coords          = miami_coords,
            gea_sigsnps     = gea_effective_sigsnps,
            gwas_sigsnps    = gwas_effective_sigsnps,
            pairwise_params = pairwise_params
        )
    })
}
