#' Association tab UI
#'
#' GEA results: combined Manhattan overlay, per-method accordion,
#' and region-centric detail panel.
#'
#' @param id module namespace id
#' @noRd
mod_association_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Combined Manhattan (always visible)
        mod_manhattan_overlay_ui(ns("combined_manhattan")),

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

        # Region selector inline above detail panel
        htmltools::div(
            class = "control-bar",
            bslib::layout_columns(
                col_widths = c(9, 3),
                shiny::uiOutput(ns("region_selector")),
                shiny::actionButton(ns("clear_region"), "Clear selection",
                                    class = "btn-sm btn-outline-secondary mt-4 w-100")
            )
        ),

        # Region detail panel (conditional)
        shiny::uiOutput(ns("region_detail_ui"))
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
        regions_data <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("regions_combined_", pd$name, "_", module), function() {
                load_regions_combined(pd$name, module)
            })
        })

        genes_data <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("genes_", pd$name, "_", module), function() {
                load_genes(pd$name, module)
            })
        })

        methods <- shiny::reactive(find_assoc_methods(project_data()$name, module))
        traits  <- shiny::reactive(find_assoc_traits(project_data()$name, module))

        # ── Region selection state ─────────────────────────────────────────────
        selected_region <- shiny::reactiveVal(NULL)

        shiny::observeEvent(input$clear_region, {
            selected_region(NULL)
            shiny::updateSelectInput(session, "region_id", selected = "")
        })

        shiny::observeEvent(input$region_id, {
            rid <- input$region_id
            selected_region(if (nzchar(rid %||% "")) rid else NULL)
        }, ignoreInit = TRUE)

        # ── Combined Manhattan ─────────────────────────────────────────────────
        manhattan_click <- mod_manhattan_overlay_server("combined_manhattan",
            project_data         = project_data,
            module               = module,
            combined             = TRUE,
            title_label          = shiny::reactive("Combined Manhattan"),
            regions              = regions_data,
            show_regions_control = TRUE
        )

        # Click on a sig SNP → update selected region
        shiny::observeEvent(manhattan_click(), {
            rid <- manhattan_click()
            if (!is.null(rid) && nzchar(rid)) {
                selected_region(rid)
                shiny::updateSelectInput(session, "region_id", selected = rid)
            }
        }, ignoreNULL = TRUE)

        # ── Region dropdown ────────────────────────────────────────────────────
        output$region_selector <- shiny::renderUI({
            rg <- regions_data()
            if (nrow(rg) == 0) {
                return(htmltools::p("No regions found.", class = "text-muted small"))
            }
            genes <- genes_data()
            labels <- build_region_labels(rg, genes)
            current <- selected_region()
            shiny::selectInput(ns("region_id"), "Select region",
                choices  = c(setNames("", ""), labels),
                selected = current %||% ""
            )
        })

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

        # Per-method Manhattan overlay
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
            regions              = regions_data,
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

        # ── Region detail panel ────────────────────────────────────────────────
        output$region_detail_ui <- shiny::renderUI({
            rid <- selected_region()
            if (is.null(rid)) return(NULL)
            mod_region_detail_ui(ns("region_detail"))
        })

        # Current trait for the selected region (first trait in region, or NULL)
        region_trait <- shiny::reactive({
            rid <- selected_region()
            rg  <- regions_data()
            if (is.null(rid) || nrow(rg) == 0) return(NULL)
            row <- rg[rg$region_id == rid, ]
            if (nrow(row) == 0) return(NULL)
            # regions_combined uses 'traits' column (comma-separated); fall back to 'trait'
            trait_col <- intersect(c("trait", "traits"), names(row))[1]
            if (is.na(trait_col)) return(NULL)
            trimws(strsplit(as.character(row[[trait_col]][1]), ",")[[1]])[1]
        })

        mod_region_detail_server("region_detail",
            project_data = project_data,
            region_id    = selected_region,
            module       = module,
            trait        = region_trait,
            genes_data   = genes_data
        )
    })
}
