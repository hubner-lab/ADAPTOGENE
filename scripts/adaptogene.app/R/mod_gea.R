#' GEA tab UI
#'
#' GEA results: combined Manhattan with interactive filter/strategy bar,
#' per-method accordion, and region-centric detail panel.
#'
#' @param id module namespace id
#' @noRd
mod_gea_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Config parameter badges
        shiny::uiOutput(ns("config_badges")),

        # Combined Manhattan — threshold bar (always visible) + filter bar (matrix/strategy)
        mod_manhattan_overlay_ui(
            ns("combined_manhattan"),
            filter_ui = htmltools::tagList(
                shiny::uiOutput(ns("threshold_bar")),
                shiny::uiOutput(ns("filter_bar"))
            )
        ),

        # Warning: traits with no significant SNPs
        shiny::uiOutput(ns("no_sig_snps_warning")),

        # Per-method accordion (collapsed by default)
        bslib::accordion(
            id       = ns("per_method_accordion"),
            open     = FALSE,
            multiple = FALSE,

            bslib::accordion_panel(
                "Per-Method Details",
                value = "per_method",
                icon  = bsicons::bs_icon("layers"),
                shiny::uiOutput(ns("method_tabs_ui")),
                shiny::uiOutput(ns("per_method_trait_ui")),
                bslib::layout_columns(
                    col_widths = c(9, 3),
                    mod_manhattan_overlay_ui(ns("method_manhattan"), height = "400px"),
                    mod_image_card_ui(ns("qq_plot"))
                )
            )
        ),

        # Interactive region explorer (replaces static region dropdown)
        mod_region_explorer_ui(ns("region_explorer"))
    )
}

#' Association tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param module character: MOD_GEA or MOD_GWAS
#' @noRd
mod_gea_server <- function(id, project_data, module = MOD_GEA) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ── Data loading ───────────────────────────────────────────────────────
        methods <- shiny::reactive(find_assoc_methods(project_data()$name, module))

        # Full p-value tables (all SNPs, all traits, wide format). Loaded once per
        # project switch; kept in reactive memo (not cachem — can exceed 200MB at WGS scale).
        all_method_pvalues <- shiny::reactive({
            pd <- project_data()
            load_all_method_pvalues(pd$name, module, pd$k_best)
        })

        all_method_wza_pvalues <- shiny::reactive({
            if (!regime_wza()) return(list())
            pd <- project_data()
            load_all_method_wza_pvalues(pd$name, module, pd$k_best)
        })

        effective_method_pvalues <- shiny::reactive({
            if (regime_wza()) all_method_wza_pvalues() else all_method_pvalues()
        })

        # All trait names (threshold-independent; for color stability + no_sig warning)
        all_trait_names_rv <- shiny::reactive({
            pv <- all_method_pvalues()
            if (length(pv) == 0) return(character(0))
            fixed <- c("SNPID", "chr", "pos", "n_snps", "mean_maf")
            sort(unique(unlist(lapply(pv, function(dt) setdiff(names(dt), fixed)))))
        })

        # Significant traits at current threshold (derived from effective_method_sigsnps)
        traits <- shiny::reactive({
            sig <- effective_method_sigsnps()
            if (length(sig) == 0) return(character(0))
            all_dt <- data.table::rbindlist(sig, use.names = TRUE, fill = TRUE)
            if (nrow(all_dt) == 0) return(character(0))
            sort(unique(all_dt$trait))
        })

        # ── Filter bar ─────────────────────────────────────────────────────────

        # Trait chip colours (Okabe-Ito, keyed on ALL traits for color stability)
        trait_colors <- shiny::reactive({
            tr <- all_trait_names_rv()
            if (length(tr) == 0) return(character(0))
            trait_color_map(tr)
        })

        # Method shapes (stable across filtering — keyed on full method list)
        method_shapes <- shiny::reactive({
            ms <- methods()
            if (length(ms) == 0) return(character(0))
            method_shape_map(ms)
        })

        # ── Trait x Method combo counts (how many sig SNPs per combo) ─────────
        combo_counts <- shiny::reactive({
            all_snps <- effective_method_sigsnps()
            counts <- list()
            for (m in names(all_snps)) {
                dt <- all_snps[[m]]
                if (nrow(dt) > 0) {
                    tc <- dt[, .(n = .N), by = "trait"]
                    for (i in seq_len(nrow(tc)))
                        counts[[paste0(tc$trait[i], "::", m)]] <- tc$n[i]
                }
            }
            counts
        })

        # Strategy always defaults to "All" (pipeline pre-computes All; user explores interactively)
        default_strategy <- shiny::reactive("All")

        active_strategy <- shiny::reactive(input$combine_strategy %||% default_strategy())

        # ── SNP clumping distance (single param for both merging and overlap) ───
        shiny::observe({
            pd      <- project_data()
            rp      <- read_region_params(pd$name)
            saved_d <- get_global_param(rp, MOD_GEA, "snp_clumping_distance")
            d <- if (!is.null(saved_d)) as.integer(saved_d)
                 else resolve_ui_snp_clumping_distance(pd$config, "GEA", pd$name)
            shiny::updateNumericInput(session, "snp_clumping_distance", value = d)
        })

        shiny::observeEvent(input$snp_clumping_distance, {
            v  <- input$snp_clumping_distance
            pd <- project_data()
            if (is.null(v) || is.na(v) || v < 1000L || is.null(pd)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, MOD_GEA, "snp_clumping_distance", as.integer(v))
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        snp_clumping_distance <- shiny::reactive({
            v <- input$snp_clumping_distance
            if (is.null(v) || is.na(v) || v < 1000L) {
                pd <- project_data()
                resolve_ui_snp_clumping_distance(pd$config, "GEA", pd$name)
            } else {
                as.integer(v)
            }
        })

        # ── Regime (per-SNP vs WZA) ────────────────────────────────────────────
        shiny::observe({
            pd    <- project_data()
            rp    <- read_region_params(pd$name)
            saved <- get_global_param(rp, MOD_GEA, "regime")
            if (!is.null(saved))
                bslib::update_switch("regime", value = isTRUE(saved), session = session)
        })

        shiny::observeEvent(input$regime, {
            pd <- project_data()
            if (is.null(pd)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, MOD_GEA, "regime", isTRUE(input$regime))
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        regime_wza <- shiny::reactive(isTRUE(input$regime))

        # ── Threshold type + value (interactive significance level) ───────────
        shiny::observe({
            pd      <- project_data()
            cfg     <- pd$config
            rp      <- read_region_params(pd$name)
            saved_t <- get_global_param(rp, module, "threshold_type")
            saved_v <- get_global_param(rp, module, "threshold_value")
            if (!is.null(saved_t)) {
                shiny::updateSelectInput(session,  "threshold_type",  selected = saved_t)
            } else {
                def <- default_threshold(cfg, module)
                shiny::updateSelectInput(session,  "threshold_type",  selected = def$type)
            }
            if (!is.null(saved_v)) {
                shiny::updateNumericInput(session, "threshold_value", value = as.numeric(saved_v))
            } else if (is.null(saved_t)) {
                def <- default_threshold(cfg, module)
                shiny::updateNumericInput(session, "threshold_value", value = def$value)
            }
        })

        shiny::observeEvent(input$threshold_type, {
            pd <- project_data(); if (is.null(pd)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, module, "threshold_type", input$threshold_type)
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$threshold_value, {
            v  <- input$threshold_value
            pd <- project_data(); if (is.null(pd) || is.null(v) || is.na(v)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, module, "threshold_value", as.numeric(v))
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        threshold_type  <- shiny::reactive(input$threshold_type  %||% "bonf")
        threshold_value_raw <- shiny::reactive({
            v <- input$threshold_value
            if (is.null(v) || is.na(v) || v <= 0) {
                default_threshold(project_data()$config, module)$value
            } else {
                as.numeric(v)
            }
        })
        # Debounce threshold value to avoid firing qval computation on every keystroke
        threshold_value <- shiny::debounce(threshold_value_raw, 500)

        # Threshold hint text (context-aware help under the value input)
        output$threshold_hint <- shiny::renderUI({
            t <- threshold_type()
            hint <- switch(t,
                bonf   = "α for Bonferroni correction",
                qval   = "FDR q-value target (0–1)",
                top    = "Number of top SNPs per trait",
                custom = "Raw p-value cutoff (e.g. 1e-5)",
                ""
            )
            htmltools::span(class = "text-muted small mt-1", hint)
        })

        # Coerce threshold value when type changes to avoid invalid combos
        shiny::observeEvent(input$threshold_type, {
            v <- input$threshold_value
            if (!threshold_value_valid_for_type(input$threshold_type, v)) {
                shiny::updateNumericInput(
                    session, "threshold_value",
                    value = threshold_value_default_for_type(input$threshold_type)
                )
            }
        }, ignoreInit = TRUE)

        # ── Always-present threshold bar (regime + threshold type/value) ─────
        output$threshold_bar <- shiny::renderUI({
            pd <- project_data()
            shiny::isolate({
                build_threshold_bar_ui(
                    ns                    = ns,
                    regime_value          = isTRUE(input$regime),
                    threshold_type_value  = input$threshold_type  %||% "bonf",
                    threshold_value_value = input$threshold_value %||%
                        default_threshold(pd$config, module)$value
                )
            })
        })

        # ── Interactive sig SNPs: computed from full pvalue tables + threshold ─
        effective_method_sigsnps <- shiny::reactive({
            pd <- project_data()
            compute_method_sigsnps_cached(
                pvalues_list = effective_method_pvalues(),
                type         = threshold_type(),
                value        = threshold_value(),
                k            = pd$k_best,
                regime       = if (regime_wza()) "wza" else "snp",
                project      = pd$name,
                module       = module
            )
        })

        # ── Threshold y-value for live threshold line on Manhattan ────────────
        # combined: min raw p threshold across all methods+traits → highest -log10 line
        combined_threshold_y <- shiny::reactive({
            pv <- effective_method_pvalues()
            if (length(pv) == 0) return(NULL)
            type  <- threshold_type()
            value <- threshold_value()
            all_thrs <- unlist(lapply(names(pv), function(m) {
                dt <- pv[[m]]
                if (nrow(dt) == 0) return(numeric(0))
                tcols <- setdiff(names(dt), c("SNPID", "chr", "pos", "n_snps", "mean_maf"))
                sapply(tcols, function(tr) {
                    r <- tryCatch(compute_pval_threshold(dt[[tr]], type, value),
                                  error = function(e) list(status = "error"))
                    if (r$status == "ok") r$threshold else NA_real_
                })
            }))
            all_thrs <- all_thrs[!is.na(all_thrs) & all_thrs > 0]
            if (length(all_thrs) == 0) return(NULL)
            -log10(min(all_thrs))  # min raw p = most stringent = highest y line
        })

        # per-method panel: threshold for selected (method, trait)
        per_method_threshold_y <- shiny::reactive({
            m <- active_method(); t <- per_method_trait()
            if (is.null(m) || is.null(t)) return(NULL)
            pv <- effective_method_pvalues()
            dt <- pv[[m]]
            if (is.null(dt) || nrow(dt) == 0 || !(t %in% names(dt))) return(NULL)
            r <- tryCatch(
                compute_pval_threshold(dt[[t]], threshold_type(), threshold_value()),
                error = function(e) list(status = "error")
            )
            if (r$status != "ok") return(NULL)
            -log10(r$threshold)
        })

        # ── I4: Config parameter badges ────────────────────────────────────────
        output$config_badges <- shiny::renderUI({
            pd <- project_data()
            k  <- pd$k_best
            if (is.na(k) && length(methods()) == 0) return(NULL)
            cfg  <- pd$config
            cmb  <- config_get(cfg, "GEA", "combine_method", default = "Union")
            cdist <- config_get(cfg, "GEA", "snp_clumping_distance",
                                default = "auto_per_chromosome")
            cdist_str <- if (grepl("^auto", cdist)) paste0(cdist, " (LD-derived)")
                         else paste0(format(as.integer(cdist), big.mark = ","), " bp")
            config_badges_bar(
                if (!is.na(k)) config_badge("K", k, "bg-primary"),
                config_badge("combine", cmb),
                config_badge("clumping", cdist_str),
                if (regime_wza()) config_badge("regime", "WZA", "bg-info")
            )
        })

        # ── D5: Traits with zero significant SNPs warning ──────────────────────
        output$no_sig_snps_warning <- shiny::renderUI({
            all_traits  <- all_trait_names_rv()
            sig_traits  <- traits()
            missing     <- setdiff(all_traits, sig_traits)
            if (length(missing) == 0) return(NULL)
            items <- lapply(missing, function(tr)
                htmltools::tags$li(htmltools::tags$code(tr)))
            htmltools::div(
                class = "alert alert-warning d-flex gap-2 align-items-start mb-2",
                bsicons::bs_icon("exclamation-triangle-fill", class = "flex-shrink-0 mt-1"),
                htmltools::div(
                    htmltools::tags$strong(
                        length(missing),
                        if (length(missing) == 1) "trait" else "traits",
                        "yielded no significant SNPs across all methods"
                    ),
                    " \u2014 not shown in filter bar. ",
                    "Adjust the Significance threshold above the Manhattan plot.",
                    htmltools::tags$ul(class = "mb-0 mt-1", items)
                )
            )
        })

        # ── Filter bar: Trait x Method matrix + Strategy + Clumping ──────────
        output$filter_bar <- shiny::renderUI({
            build_filter_bar_ui(
                ns                          = ns,
                traits                      = traits(),
                methods                     = methods(),
                trait_colors                = trait_colors(),
                combo_counts                = combo_counts(),
                default_strategy_value      = default_strategy(),
                snp_clumping_distance_value = snp_clumping_distance()
            )
        })

        # ── Interactive sig SNPs (combined + filtered by matrix selection) ─────
        # effective_method_sigsnps() is already threshold-filtered by compute_method_sigsnps_cached.
        # compute_interactive_sigsnps further applies matrix filter, strategy, and clumping.
        interactive_sigsnps <- shiny::reactive({
            compute_interactive_sigsnps(
                all_method_sigsnps  = effective_method_sigsnps(),
                tm_selection_json   = input$tm_selection,
                combo_counts        = combo_counts(),
                known_traits        = all_trait_names_rv(),
                strategy            = active_strategy(),
                clumping_distance   = snp_clumping_distance(),
                project_name        = project_data()$name,
                module              = module
            )
        })

        # ── WZA path overrides for Manhattan ──────────────────────────────────
        combined_wza_bg <- shiny::reactive({
            if (!regime_wza()) return(NULL)
            pd <- project_data(); k <- pd$k_best
            if (is.na(k)) return(NULL)
            combined_manhattan_wza_bg_path(pd$name, module, k)
        })
        combined_wza_coords <- shiny::reactive({
            if (!regime_wza()) return(NULL)
            pd <- project_data(); k <- pd$k_best
            if (is.na(k)) return(NULL)
            combined_manhattan_wza_coords_path(pd$name, module, k)
        })

        # ── Interactive region explorer ────────────────────────────────────────
        explorer <- mod_region_explorer_server("region_explorer",
            project_data        = project_data,
            module              = module,
            interactive_sigsnps = interactive_sigsnps,
            region_distance     = snp_clumping_distance,
            regime              = regime_wza
        )

        # ── Combined Manhattan ─────────────────────────────────────────────────
        manhattan_click <- mod_manhattan_overlay_server("combined_manhattan",
            project_data         = project_data,
            module               = module,
            combined             = TRUE,
            title_label          = shiny::reactive({
                if (regime_wza()) "Combined Manhattan (WZA)" else "Combined Manhattan"
            }),
            regions              = explorer$computed_regions,
            current_region_id    = explorer$selected_region_id,
            show_regions_control = FALSE,
            sig_snps_override    = interactive_sigsnps,
            trait_colors         = trait_colors,
            method_shapes        = method_shapes,
            bg_path_override     = combined_wza_bg,
            coords_path_override = combined_wza_coords,
            threshold_y          = combined_threshold_y
        )

        # Manhattan SNP click → select the enclosing region in the explorer
        shiny::observeEvent(manhattan_click(), {
            rid <- manhattan_click()
            if (!is.null(rid) && nzchar(rid)) {
                explorer$selected_region_id(rid)
            }
        }, ignoreNULL = TRUE)

        # ── Per-method tabs UI ─────────────────────────────────────────────────
        output$method_tabs_ui <- shiny::renderUI({
            ms <- methods()
            if (length(ms) == 0) return(NULL)
            panels <- lapply(ms, function(m) {
                bslib::nav_panel(m, value = m)
            })
            do.call(bslib::navset_card_underline, c(list(id = ns("method_tab")), panels))
        })

        # Per-method trait selector UI — uses all trait names so user can view any
        # trait regardless of current threshold (useful for QC / comparison)
        output$per_method_trait_ui <- shiny::renderUI({
            tr <- all_trait_names_rv()
            if (length(tr) == 0) return(NULL)
            shiny::selectInput(ns("per_method_trait"), "Trait",
                               choices = tr, selected = tr[1], width = "200px")
        })

        per_method_trait <- shiny::reactive({
            tr <- all_trait_names_rv()
            if (length(tr) == 0) return(NULL)
            input$per_method_trait %||% tr[1]
        })

        active_method <- shiny::reactive(input$method_tab %||% methods()[1])

        # Per-method Manhattan overlay — use interactive sig SNPs (avoids pipeline file read)
        per_method_sigsnps_override <- shiny::reactive({
            m <- active_method(); t <- per_method_trait()
            all_ms <- effective_method_sigsnps()
            if (is.null(m) || is.null(t) || length(all_ms) == 0) return(data.table::data.table())
            dt <- all_ms[[m]]
            if (is.null(dt) || nrow(dt) == 0) return(data.table::data.table())
            dt[trait == t]
        })
        mod_manhattan_overlay_server("method_manhattan",
            project_data         = project_data,
            module               = module,
            method               = active_method,
            trait                = per_method_trait,
            combined             = FALSE,
            sig_snps_override    = per_method_sigsnps_override,
            title_label          = shiny::reactive({
                m <- active_method()
                t <- per_method_trait()
                if (!is.null(m) && !is.null(t)) paste0(m, " — ", t)
                else "Method Manhattan"
            }),
            show_regions_control = FALSE,
            threshold_y          = per_method_threshold_y
        )

        # QQ plot for current method + trait
        qq_path <- shiny::reactive({
            m   <- active_method()
            t   <- per_method_trait()
            pd  <- project_data()
            k   <- pd$k_best
            if (is.null(m) || is.null(t) || is.na(k)) return(NULL)
            adj <- resolve_adjust(pd$config, m,
                if (module == MOD_GWAS) MOD_GWAS else MOD_GEA)
            if (is.null(adj)) return(NULL)
            qq_plot_path(pd$name, module, m, t, k, adj)
        })

        mod_image_card_server("qq_plot",
            path    = qq_path,
            title   = shiny::reactive("QQ Plot"),
            dl_name = shiny::reactive({
                paste0("qq_", active_method() %||% "method",
                       "_", per_method_trait() %||% "trait")
            })
        )

    })
}
