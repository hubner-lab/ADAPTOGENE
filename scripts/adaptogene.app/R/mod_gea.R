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

        # Combined Manhattan — filter bar injected inside the card
        mod_manhattan_overlay_ui(
            ns("combined_manhattan"),
            filter_ui = shiny::uiOutput(ns("filter_bar"))
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
        traits  <- shiny::reactive(find_assoc_traits(project_data()$name, module))

        # Load ALL per-method sig SNPs once per project (for interactive combine)
        all_method_sigsnps <- shiny::reactive({
            pd <- project_data()
            load_all_method_sigsnps(pd$name, module)
        })

        # ── Filter bar ─────────────────────────────────────────────────────────

        # Trait chip colours (Okabe-Ito, sorted order matching fct_manhattan.R)
        trait_colors <- shiny::reactive({
            tr <- traits()
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
            all_snps <- all_method_sigsnps()
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

        # ── WZA sig windows (loaded when regime=wza) ──────────────────────────
        all_method_wza_sigsnps <- shiny::reactive({
            if (!regime_wza()) return(list())
            pd <- project_data()
            load_all_method_wza_sigwindows(pd$name, module)
        })

        # Regime-aware source of sig SNPs (snp or wza)
        effective_method_sigsnps <- shiny::reactive({
            if (regime_wza()) all_method_wza_sigsnps() else all_method_sigsnps()
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
            pd          <- project_data()
            all_traits  <- load_all_trait_names(pd$name, module)
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
                    "Try relaxing thresholds in ",
                    htmltools::tags$code("GEA.configs"), ".",
                    htmltools::tags$ul(class = "mb-0 mt-1", items)
                )
            )
        })

        # ── Filter bar: Regime + Trait x Method matrix + Strategy + Distance ──
        output$filter_bar <- shiny::renderUI({
            build_filter_bar_ui(
                ns                          = ns,
                traits                      = traits(),
                methods                     = methods(),
                trait_colors                = trait_colors(),
                combo_counts                = combo_counts(),
                default_strategy_value      = default_strategy(),
                snp_clumping_distance_value = snp_clumping_distance(),
                regime_value                = regime_wza()
            )
        })

        # ── Interactive sig SNPs (combined + filtered by matrix selection) ─────
        interactive_sigsnps <- shiny::reactive({
            compute_interactive_sigsnps(
                all_method_sigsnps  = effective_method_sigsnps(),
                tm_selection_json   = input$tm_selection,
                combo_counts        = combo_counts(),
                known_traits        = traits(),
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
            coords_path_override = combined_wza_coords
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

        # Per-method trait selector UI
        output$per_method_trait_ui <- shiny::renderUI({
            tr <- traits()
            if (length(tr) == 0) return(NULL)
            shiny::selectInput(ns("per_method_trait"), "Trait",
                               choices = tr, selected = tr[1], width = "200px")
        })

        per_method_trait <- shiny::reactive({
            tr <- traits()
            if (length(tr) == 0) return(NULL)
            input$per_method_trait %||% tr[1]
        })

        active_method <- shiny::reactive(input$method_tab %||% methods()[1])

        # Per-method Manhattan overlay (no region shapes — shows method-specific data only)
        mod_manhattan_overlay_server("method_manhattan",
            project_data         = project_data,
            module               = module,
            method               = active_method,
            trait                = per_method_trait,
            combined             = FALSE,
            title_label          = shiny::reactive({
                m <- active_method()
                t <- per_method_trait()
                if (!is.null(m) && !is.null(t)) paste0(m, " — ", t)
                else "Method Manhattan"
            }),
            show_regions_control = FALSE
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
