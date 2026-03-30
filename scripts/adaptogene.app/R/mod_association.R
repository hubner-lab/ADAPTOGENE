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

        # Strategy from config (default for radio button)
        default_strategy <- shiny::reactive({
            pd <- project_data()
            ds <- config_get(pd$config, "association", "combine_method", default = "Sum")
            if (!ds %in% c("Sum", "Overlap", "PairOverlap")) "Sum" else ds
        })

        active_strategy <- shiny::reactive(input$combine_strategy %||% default_strategy())

        # ── Region distance (owned here, passed to explorer) ───────────────────
        shiny::observe({
            pd <- project_data()
            d  <- config_get(pd$config, "association", "region_distance", default = 2000000L)
            shiny::updateNumericInput(session, "region_distance", value = as.integer(d))
        })

        region_distance <- shiny::reactive({
            v <- input$region_distance
            if (is.null(v) || is.na(v) || v < 1000L) 2000000L else as.integer(v)
        })

        # ── Filter bar: Trait x Method matrix + Strategy + Distance ───────────
        output$filter_bar <- shiny::renderUI({
            tr   <- traits()
            ms   <- methods()
            cols <- trait_colors()
            cnts <- combo_counts()

            if (length(tr) == 0) return(NULL)

            # Build table header row: blank + method column headers
            header_cells <- lapply(ms, function(m) {
                htmltools::tags$th(
                    class = "tm-col-header",
                    `data-method` = m,
                    m
                )
            })
            header_row <- htmltools::tags$tr(
                htmltools::tags$th(),  # blank corner
                header_cells
            )

            # Build body rows: one per trait
            body_rows <- lapply(tr, function(t) {
                dot_html <- paste0(
                    '<span class="filter-trait-dot" style="background:', cols[t],
                    ';display:inline-block;width:8px;height:8px;border-radius:50%;',
                    'margin-right:4px;"></span>'
                )
                row_header <- htmltools::tags$th(
                    class = "tm-row-header",
                    `data-trait` = t,
                    htmltools::HTML(paste0(dot_html, htmltools::htmlEscape(t)))
                )
                cells <- lapply(ms, function(m) {
                    key <- paste0(t, "::", m)
                    n   <- cnts[[key]] %||% 0L
                    if (n > 0) {
                        htmltools::tags$td(
                            htmltools::tags$button(
                                class = "tm-cell tm-active",
                                `data-trait` = t,
                                `data-method` = m,
                                as.character(n)
                            )
                        )
                    } else {
                        htmltools::tags$td(
                            htmltools::tags$button(
                                class = "tm-cell tm-empty",
                                `data-trait` = t,
                                `data-method` = m,
                                "—"
                            )
                        )
                    }
                })
                htmltools::tags$tr(row_header, cells)
            })

            matrix_table <- htmltools::tags$table(
                class = "tm-matrix",
                htmltools::tags$thead(header_row),
                htmltools::tags$tbody(body_rows)
            )

            # Hidden input bridge updated by JS
            hidden_input <- shiny::textInput(
                ns("tm_selection"), label = NULL, value = ""
            )

            # JS: delegated click on matrix container
            container_id <- ns("tm_container")
            input_id     <- ns("tm_selection")
            js_code <- sprintf('
(function() {
    var container = document.getElementById("%s");
    if (!container) return;
    function syncSelection() {
        var active = container.querySelectorAll(".tm-cell.tm-active");
        var pairs = Array.from(active).map(function(el) {
            return el.dataset.trait + "::" + el.dataset.method;
        });
        Shiny.setInputValue("%s", JSON.stringify(pairs), {priority: "event"});
    }
    container.addEventListener("click", function(e) {
        var cell = e.target.closest(".tm-cell");
        if (cell && !cell.classList.contains("tm-empty")) {
            cell.classList.toggle("tm-active");
            syncSelection();
            return;
        }
        var rh = e.target.closest(".tm-row-header");
        if (rh) {
            var trait = rh.dataset.trait;
            var cells = container.querySelectorAll(".tm-cell[data-trait=\\"" + trait + "\\"]:not(.tm-empty)");
            var allActive = Array.from(cells).every(function(c) { return c.classList.contains("tm-active"); });
            cells.forEach(function(c) { c.classList.toggle("tm-active", !allActive); });
            syncSelection();
            return;
        }
        var ch = e.target.closest(".tm-col-header");
        if (ch) {
            var method = ch.dataset.method;
            var cells = container.querySelectorAll(".tm-cell[data-method=\\"" + method + "\\"]:not(.tm-empty)");
            var allActive = Array.from(cells).every(function(c) { return c.classList.contains("tm-active"); });
            cells.forEach(function(c) { c.classList.toggle("tm-active", !allActive); });
            syncSelection();
        }
    });
    syncSelection();
})();
', container_id, input_id)

            htmltools::div(
                class = "manhattan-filter-bar",
                htmltools::div(
                    class = "filter-row align-items-start",
                    # Matrix
                    htmltools::div(
                        id    = container_id,
                        class = "tm-container me-4",
                        matrix_table,
                        hidden_input
                    ),
                    # Strategy
                    htmltools::div(
                        class = "d-flex flex-column me-4",
                        htmltools::span("Strategy", class = "filter-label mb-1"),
                        shiny::radioButtons(
                            ns("combine_strategy"), label = NULL,
                            choices  = c("Sum", "Overlap", "PairOverlap"),
                            selected = default_strategy(),
                            inline   = FALSE
                        )
                    ),
                    # Distance
                    htmltools::div(
                        class = "d-flex flex-column",
                        htmltools::span("Distance (bp)", class = "filter-label mb-1"),
                        shiny::numericInput(
                            ns("region_distance"), label = NULL,
                            value = 2000000L, min = 1000L, step = 100000L,
                            width = "130px"
                        )
                    )
                ),
                htmltools::tags$script(htmltools::HTML(js_code))
            )
        })

        # ── Interactive sig SNPs (combined + filtered by matrix selection) ─────
        interactive_sigsnps <- shiny::reactive({
            pd       <- project_data()
            all_snps <- all_method_sigsnps()

            # Fall back when no per-method files exist
            if (length(all_snps) == 0) return(NULL)

            strategy <- active_strategy()
            gap      <- config_get(pd$config, "association", "combine_gap",
                                   default = 200000L)

            # Parse selected trait::method pairs from matrix
            # Restrict to traits shown in the matrix (climate/GEA traits only —
            # per-method files may also contain phenotypic traits from a combined run)
            known_traits <- traits()
            all_valid_pairs <- names(combo_counts())[
                vapply(combo_counts(), function(n) n > 0, logical(1))
            ]
            all_valid_pairs <- all_valid_pairs[
                sub("::.*", "", all_valid_pairs) %in% known_traits
            ]

            sel_json <- input$tm_selection
            if (is.null(sel_json) || !nzchar(sel_json)) {
                selected_pairs <- all_valid_pairs
            } else {
                selected_pairs <- tryCatch(
                    jsonlite::fromJSON(sel_json),
                    error = function(e) all_valid_pairs
                )
            }

            # Filter each method to selected traits for that method
            filtered <- lapply(names(all_snps), function(m) {
                paired_traits <- sub("::.*", "",
                    selected_pairs[grepl(paste0("::", m, "$"), selected_pairs, fixed = FALSE)])
                if (length(paired_traits) == 0) return(NULL)
                dt <- all_snps[[m]]
                if (nrow(dt) == 0) return(NULL)
                dt[trait %in% paired_traits]
            })
            names(filtered) <- names(all_snps)
            filtered <- Filter(function(x) !is.null(x) && nrow(x) > 0, filtered)

            if (length(filtered) == 0) return(.empty_sigsnps_assoc())

            combined_dt <- combine_sigsnps(filtered, strategy = strategy,
                                           gap = as.integer(gap))

            if (nrow(combined_dt) == 0) return(.empty_sigsnps_assoc())

            combined_dt <- assign_region_ids(combined_dt, pd$name, module)
            combined_dt
        })

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
            title_label          = shiny::reactive("Combined Manhattan"),
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
