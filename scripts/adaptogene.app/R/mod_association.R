#' Association tab UI
#'
#' GEA results: combined Manhattan with interactive filter/strategy bar,
#' per-method accordion, and region-centric detail panel.
#'
#' @param id module namespace id
#' @noRd
mod_association_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
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

#' Association tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param module character: MOD_ASSOC or MOD_PHENO
#' @noRd
mod_association_server <- function(id, project_data, module = MOD_ASSOC) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ── Data loading ───────────────────────────────────────────────────────
        methods <- shiny::reactive(find_assoc_methods(project_data()$name, module))
        traits  <- shiny::reactive(find_assoc_traits(project_data()$name, module))

        # Load ALL per-method sig SNPs once per project (for interactive combine)
        all_method_sigsnps <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("all_method_sigsnps_", pd$name, "_", module), function() {
                load_all_method_sigsnps(pd$name, module)
            })
        })

        # ── Filter bar ─────────────────────────────────────────────────────────

        # Trait chip colours (Okabe-Ito, sorted order matching fct_manhattan.R)
        trait_colors <- shiny::reactive({
            tr <- traits()
            if (length(tr) == 0) return(character(0))
            trait_color_map(tr)
        })

        output$filter_bar <- shiny::renderUI({
            tr   <- traits()
            ms   <- methods()
            cols <- trait_colors()
            pd   <- project_data()

            if (length(tr) == 0) return(NULL)

            # Default strategy from config
            default_strategy <- config_get(pd$config, "association", "combine_method",
                                           default = "Sum")
            # Normalise: only Sum/Overlap/PairOverlap are valid interactive strategies
            valid_strategies <- c("Sum", "Overlap", "PairOverlap")
            if (!default_strategy %in% valid_strategies) default_strategy <- "Sum"

            # Trait checkboxes with colored dot HTML labels
            trait_choices <- stats::setNames(
                lapply(tr, function(t) {
                    htmltools::HTML(paste0(
                        '<span class="filter-trait-dot" style="background:', cols[t], '"></span>',
                        htmltools::htmlEscape(t)
                    ))
                }),
                tr
            )

            htmltools::div(
                class = "manhattan-filter-bar",

                # Row 1: trait toggles
                htmltools::div(
                    class = "filter-row",
                    htmltools::span("Traits", class = "filter-label"),
                    shiny::checkboxGroupInput(
                        ns("filter_traits"), label = NULL,
                        choices  = trait_choices,
                        selected = tr,
                        inline   = TRUE
                    )
                ),

                # Row 2: method toggles + strategy
                htmltools::div(
                    class = "filter-row",
                    htmltools::span("Methods", class = "filter-label"),
                    shiny::checkboxGroupInput(
                        ns("filter_methods"), label = NULL,
                        choices  = ms,
                        selected = ms,
                        inline   = TRUE
                    ),
                    htmltools::span("Strategy", class = "filter-label ms-3"),
                    shiny::radioButtons(
                        ns("combine_strategy"), label = NULL,
                        choices  = c("Sum", "Overlap", "PairOverlap"),
                        selected = default_strategy,
                        inline   = TRUE
                    )
                )
            )
        })

        # Reactive filter state
        active_traits   <- shiny::reactive(input$filter_traits   %||% traits())
        active_methods  <- shiny::reactive(input$filter_methods  %||% methods())
        active_strategy <- shiny::reactive(input$combine_strategy %||% "Sum")

        # ── Interactive sig SNPs (combined + filtered) ─────────────────────────
        interactive_sigsnps <- shiny::reactive({
            pd       <- project_data()
            all_snps <- all_method_sigsnps()

            # Fall back to precomputed file when no per-method files exist
            if (length(all_snps) == 0) return(NULL)

            tr_sel   <- active_traits()
            ms_sel   <- active_methods()
            strategy <- active_strategy()
            gap      <- config_get(pd$config, "association", "combine_gap",
                                   default = 200000L)

            # Apply trait × method filter to each method's data
            filtered <- lapply(names(all_snps), function(m) {
                if (!m %in% ms_sel) return(NULL)
                dt <- all_snps[[m]]
                if (nrow(dt) == 0) return(NULL)
                dt[trait %in% tr_sel]
            })
            names(filtered) <- names(all_snps)
            filtered <- Filter(function(x) !is.null(x) && nrow(x) > 0, filtered)

            if (length(filtered) == 0) return(.empty_sigsnps_assoc())

            # Combine using selected strategy
            combined_dt <- combine_sigsnps(filtered, strategy = strategy,
                                           gap = as.integer(gap))

            if (nrow(combined_dt) == 0) return(.empty_sigsnps_assoc())

            # Assign region IDs from precomputed regions_combined.tsv
            combined_dt <- assign_region_ids(combined_dt, pd$name, module)
            combined_dt
        })

        # ── Interactive region explorer ────────────────────────────────────────
        explorer <- mod_region_explorer_server("region_explorer",
            project_data        = project_data,
            module              = module,
            interactive_sigsnps = interactive_sigsnps
        )

        # ── Combined Manhattan ─────────────────────────────────────────────────
        manhattan_click <- mod_manhattan_overlay_server("combined_manhattan",
            project_data         = project_data,
            module               = module,
            combined             = TRUE,
            title_label          = shiny::reactive("Combined Manhattan"),
            regions              = explorer$computed_regions,
            current_region_id    = explorer$selected_region_id,
            show_regions_control = FALSE,
            sig_snps_override    = interactive_sigsnps
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
                if (module == MOD_PHENO) "phenotype_association" else "association")
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

# Helper: empty data.table matching the expected sig_snps schema
.empty_sigsnps_assoc <- function() {
    data.table::data.table(
        SNPID     = character(),
        chr       = character(),
        pos       = integer(),
        pvalue    = numeric(),
        method    = character(),
        trait     = character(),
        region_id = character()
    )
}
