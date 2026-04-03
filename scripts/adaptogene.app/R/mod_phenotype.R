#' Phenotype Association tab UI
#'
#' GWAS results: phenotype piemap + combined Manhattan with interactive
#' filter/strategy bar + per-method accordion + region explorer.
#'
#' @param id module namespace id
#' @noRd
mod_phenotype_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Trait selector inline above phenomap
        htmltools::div(
            class = "control-bar",
            shiny::uiOutput(ns("trait_selector"))
        ),

        # Phenomap piemap (constrained to 2/3 width)
        bslib::layout_columns(
            col_widths = c(8, 4),
            bslib::card(
                full_screen = TRUE,
                bslib::card_header(
                    bsicons::bs_icon("geo-fill"),
                    " Phenotype Map"
                ),
                bslib::card_body(
                    class = "piemap-container",
                    shiny::uiOutput(ns("phenomap_content"))
                )
            ),
            bslib::card(
                bslib::card_body(
                    class = "text-muted small",
                    htmltools::p(bsicons::bs_icon("info-circle"), " Click a significant SNP in the Manhattan plot below to explore GO enrichment and genes for that region."),
                    htmltools::p(bsicons::bs_icon("arrows-fullscreen"), " Use the expand icon on any card for full-screen detail view.")
                )
            )
        ),

        # Combined Manhattan — filter bar injected inside the card
        mod_manhattan_overlay_ui(
            ns("combined_manhattan"),
            filter_ui = shiny::uiOutput(ns("filter_bar"))
        ),

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

#' Phenotype Association tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_phenotype_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        module <- MOD_PHENO

        # ── Data loading ───────────────────────────────────────────────────────
        methods <- shiny::reactive(find_assoc_methods(project_data()$name, module))
        traits  <- shiny::reactive(find_assoc_traits(project_data()$name, module))

        # Load ALL per-method sig SNPs once per project (for interactive combine)
        all_method_sigsnps <- shiny::reactive({
            pd <- project_data()
            load_all_method_sigsnps(pd$name, module)
        })

        # ── Filter bar ─────────────────────────────────────────────────────────

        trait_colors <- shiny::reactive({
            tr <- traits()
            if (length(tr) == 0) return(character(0))
            trait_color_map(tr)
        })

        method_shapes <- shiny::reactive({
            ms <- methods()
            if (length(ms) == 0) return(character(0))
            method_shape_map(ms)
        })

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

        # Strategy: try phenotype_association section, fall back to association
        default_strategy <- shiny::reactive({
            pd <- project_data()
            ds <- .normalize_strategy(config_get(pd$config, "phenotype_association", "combine_method",
                             default = config_get(pd$config, "association", "combine_method",
                                                  default = "All")))
            if (!ds %in% c("All", "Overlap", "MethodOverlap")) "All" else ds
        })

        active_strategy <- shiny::reactive(input$combine_strategy %||% default_strategy())

        # ── Region distance (owned here, passed to explorer) ───────────────────
        shiny::observe({
            pd      <- project_data()
            rp      <- read_region_params(pd$name)
            saved_d <- get_global_param(rp, MOD_PHENO, "region_distance")
            d <- if (!is.null(saved_d)) as.integer(saved_d)
                 else as.integer(config_get(pd$config, "phenotype_association", "region_distance",
                                  default = config_get(pd$config, "association", "region_distance",
                                                       default = 2000000L)))
            shiny::updateNumericInput(session, "region_distance", value = d)
        })

        shiny::observeEvent(input$region_distance, {
            v  <- input$region_distance
            pd <- project_data()
            if (is.null(v) || is.na(v) || v < 1000L || is.null(pd)) return()
            rp <- read_region_params(pd$name)
            rp <- set_global_param(rp, MOD_PHENO, "region_distance", as.integer(v))
            save_region_params(pd$name, rp)
        }, ignoreInit = TRUE)

        region_distance <- shiny::reactive({
            v <- input$region_distance
            if (is.null(v) || is.na(v) || v < 1000L) {
                pd <- project_data()
                as.integer(config_get(
                    pd$config, "phenotype_association", "region_distance",
                    default = config_get(pd$config, "association", "region_distance",
                                        default = 2000000L)))
            } else {
                as.integer(v)
            }
        })

        # ── Filter bar UI ──────────────────────────────────────────────────────
        output$filter_bar <- shiny::renderUI({
            build_filter_bar_ui(
                ns                     = ns,
                traits                 = traits(),
                methods                = methods(),
                trait_colors           = trait_colors(),
                combo_counts           = combo_counts(),
                default_strategy_value = default_strategy(),
                region_distance_value  = region_distance()
            )
        })

        # ── Interactive sig SNPs ───────────────────────────────────────────────
        interactive_sigsnps <- shiny::reactive({
            pd  <- project_data()
            gap <- config_get(pd$config, "phenotype_association", "combine_gap",
                              default = config_get(pd$config, "association", "combine_gap",
                                                   default = 200000L))
            compute_interactive_sigsnps(
                all_method_sigsnps = all_method_sigsnps(),
                tm_selection_json  = input$tm_selection,
                combo_counts       = combo_counts(),
                known_traits       = traits(),
                strategy           = active_strategy(),
                gap                = gap,
                project_name       = pd$name,
                module             = module
            )
        })

        # ── Phenotype trait selector (for phenomap) ────────────────────────────
        output$trait_selector <- shiny::renderUI({
            tr <- traits()
            if (length(tr) == 0)
                return(shiny::p("No traits found.", class = "text-muted small"))
            shiny::selectInput(ns("pheno_trait"), "Trait (phenomap)",
                               choices = tr, selected = tr[1])
        })

        selected_pheno_trait <- shiny::reactive(input$pheno_trait %||% traits()[1])

        # ── Phenomap ───────────────────────────────────────────────────────────
        output$phenomap_content <- shiny::renderUI({
            tr <- selected_pheno_trait()
            if (is.null(tr)) return(plot_placeholder("Select a trait"))
            pd   <- project_data()
            path <- pheno_piemap_path(pd$name, tr)
            if (file_ok(path)) {
                shiny::imageOutput(ns("phenomap_img"), height = "auto", width = "100%")
            } else {
                plot_placeholder("Phenotype map not available",
                    "Run mode=association_phenotypes to generate phenotype piemaps")
            }
        })

        output$phenomap_img <- shiny::renderImage({
            tr <- shiny::req(selected_pheno_trait())
            pd <- project_data()
            p  <- pheno_piemap_path(pd$name, tr)
            shiny::validate(shiny::need(file_ok(p), "Phenomap not found"))
            list(src = p, contentType = "image/png", width = "100%", alt = paste("Phenomap", tr))
        }, deleteFile = FALSE)

        # ── Interactive region explorer ────────────────────────────────────────
        explorer <- mod_region_explorer_server("region_explorer",
            project_data        = project_data,
            module              = module,
            interactive_sigsnps = interactive_sigsnps,
            region_distance     = region_distance
        )

        # ── Combined Manhattan ─────────────────────────────────────────────────
        manhattan_click <- mod_manhattan_overlay_server("combined_manhattan",
            project_data         = project_data,
            module               = module,
            combined             = TRUE,
            title_label          = shiny::reactive("Combined GWAS Manhattan"),
            regions              = explorer$computed_regions,
            current_region_id    = explorer$selected_region_id,
            show_regions_control = FALSE,
            sig_snps_override    = interactive_sigsnps,
            trait_colors         = trait_colors,
            method_shapes        = method_shapes
        )

        # Manhattan SNP click → select the enclosing region in the explorer
        shiny::observeEvent(manhattan_click(), {
            rid <- manhattan_click()
            if (!is.null(rid) && nzchar(rid)) {
                explorer$selected_region_id(rid)
            }
        }, ignoreNULL = TRUE)

        # ── Per-method tabs ────────────────────────────────────────────────────
        output$method_tabs_ui <- shiny::renderUI({
            ms <- methods()
            if (length(ms) == 0) return(NULL)
            panels <- lapply(ms, function(m) bslib::nav_panel(m, value = m))
            do.call(bslib::navset_card_underline, c(list(id = ns("method_tab")), panels))
        })

        active_method <- shiny::reactive(input$method_tab %||% methods()[1])

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

        mod_manhattan_overlay_server("method_manhattan",
            project_data         = project_data,
            module               = module,
            method               = active_method,
            trait                = per_method_trait,
            combined             = FALSE,
            title_label          = shiny::reactive({
                m <- active_method()
                t <- per_method_trait()
                if (!is.null(m) && !is.null(t)) paste0(m, " — ", t) else "Method Manhattan"
            }),
            show_regions_control = FALSE
        )

        qq_path <- shiny::reactive({
            m   <- active_method()
            t   <- per_method_trait()
            pd  <- project_data()
            k   <- pd$k_best
            if (is.null(m) || is.null(t) || is.na(k)) return(NULL)
            adj <- resolve_adjust(pd$config, m, "phenotype_association")
            if (is.null(adj)) return(NULL)
            qq_plot_path(pd$name, module, m, t, k, adj)
        })

        mod_image_card_server("qq_plot",
            path    = qq_path,
            title   = shiny::reactive("QQ Plot"),
            dl_name = shiny::reactive(paste0("qq_", active_method() %||% "method",
                                              "_", per_method_trait() %||% "trait"))
        )

    })
}
