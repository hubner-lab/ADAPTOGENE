#' Overlapping Regions tab UI
#'
#' Miami Manhattan with interactive filter/strategy bar, overlap pairs table,
#' region explorer, and pairwise trait overlap section.
#'
#' @param id module namespace id
#' @noRd
mod_overlapping_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Miami Manhattan — filter bar injected inside the card
        mod_manhattan_overlay_ui(
            ns("miami"),
            height     = "420px",
            filter_ui  = shiny::uiOutput(ns("filter_bar"))
        ),

        # Overlap pairs table
        bslib::card(
            full_screen = TRUE,
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                htmltools::span(
                    bsicons::bs_icon("table"),
                    " Overlap Summary"
                ),
                shiny::downloadButton(ns("dl_overlap"), "CSV",
                                      class = "btn-sm btn-outline-secondary")
            ),
            bslib::card_body(
                DT::DTOutput(ns("overlap_table"))
            )
        ),

        # Interactive region explorer (replaces static region dropdown)
        mod_region_explorer_ui(ns("region_explorer")),

        # Pairwise Trait Overlap section
        htmltools::hr(class = "my-4"),
        mod_pairwise_overlap_ui(ns("pairwise"))
    )
}

#' Overlapping Regions tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_overlapping_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        module <- MOD_OVERLAP

        # ── Data loading ───────────────────────────────────────────────────────
        # Methods and traits come from both GEA and GWAS modules
        methods <- shiny::reactive({
            nm <- project_data()$name
            union(find_assoc_methods(nm, MOD_ASSOC), find_assoc_methods(nm, MOD_PHENO))
        })

        traits <- shiny::reactive({
            nm <- project_data()$name
            union(find_assoc_traits(nm, MOD_ASSOC), find_assoc_traits(nm, MOD_PHENO))
        })

        # Load per-method sig SNPs from both GEA + GWAS modules
        all_method_sigsnps <- shiny::reactive({
            pd <- project_data()
            load_all_overlap_method_sigsnps(pd$name)
        })

        overlap_data <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("overlap_summary_", pd$name), function() {
                load_overlap_summary(pd$name)
            })
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

        default_strategy <- shiny::reactive({
            pd <- project_data()
            ds <- config_get(pd$config, "overlap", "combine_method",
                             default = "Sum")
            if (!ds %in% c("Sum", "Overlap", "PairOverlap")) "Sum" else ds
        })

        active_strategy <- shiny::reactive(input$combine_strategy %||% default_strategy())

        # ── Region distance (owned here, passed to explorer) ───────────────────
        shiny::observe({
            pd <- project_data()
            d  <- config_get(pd$config, "overlap", "region_distance",
                             default = config_get(pd$config, "association", "region_distance",
                                                  default = 2000000L))
            shiny::updateNumericInput(session, "region_distance", value = as.integer(d))
        })

        region_distance <- shiny::reactive({
            v <- input$region_distance
            if (is.null(v) || is.na(v) || v < 1000L) 2000000L else as.integer(v)
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
            gap <- config_get(pd$config, "overlap", "combine_gap", default = 200000L)
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

        # ── Interactive region explorer ────────────────────────────────────────
        explorer <- mod_region_explorer_server("region_explorer",
            project_data        = project_data,
            module              = module,
            interactive_sigsnps = interactive_sigsnps,
            region_distance     = region_distance
        )

        # ── Miami Manhattan ────────────────────────────────────────────────────
        miami_click <- mod_manhattan_overlay_server("miami",
            project_data         = project_data,
            module               = module,
            is_miami             = TRUE,
            combined             = TRUE,
            title_label          = shiny::reactive("Miami Plot (GEA \u2191 | GWAS \u2193)"),
            regions              = explorer$computed_regions,
            current_region_id    = explorer$selected_region_id,
            show_regions_control = FALSE,
            sig_snps_override    = interactive_sigsnps,
            trait_colors         = trait_colors,
            method_shapes        = method_shapes
        )

        # Miami SNP click → select the enclosing region in the explorer
        shiny::observeEvent(miami_click(), {
            rid <- miami_click()
            if (!is.null(rid) && nzchar(rid)) {
                explorer$selected_region_id(rid)
            }
        }, ignoreNULL = TRUE)

        # ── Overlap table ──────────────────────────────────────────────────────
        output$overlap_table <- DT::renderDataTable({
            ov <- overlap_data()
            # Filter out SUMMARY rows for main table display
            if (nrow(ov) > 0 && "region_id" %in% names(ov)) {
                ov <- ov[ov$region_id != "SUMMARY", ]
            }
            safe_datatable(
                as.data.frame(ov),
                extensions = "Buttons",
                options    = list(
                    dom       = "Bfrtip",
                    buttons   = list("csv"),
                    scrollX   = TRUE,
                    pageLength = 15
                ),
                rownames = FALSE
            )
        })

        output$dl_overlap <- shiny::downloadHandler(
            filename = function() paste0("overlap_summary_", project_data()$name, ".csv"),
            content  = function(file) {
                utils::write.csv(as.data.frame(overlap_data()), file, row.names = FALSE)
            }
        )

        # ── Signal Comparison ──────────────────────────────────────────────────
        miami_coords <- shiny::reactive({
            pd <- project_data()
            k  <- pd$k_best
            if (is.na(k)) return(NULL)
            load_cached(paste0("coords_miami_", pd$name, "_", k), function() {
                load_coords(miami_coords_path(pd$name, k))
            })
        })

        mod_pairwise_overlap_server("pairwise",
            project_data = project_data,
            coords       = miami_coords
        )
    })
}
