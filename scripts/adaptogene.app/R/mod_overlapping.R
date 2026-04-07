#' Overlapping Regions tab UI
#'
#' Two independent filter bars (GEA above Miami, GWAS below) drive independent
#' region sets on each half of the Miami plot. Overlap regions — where a GEA region
#' and a GWAS region intersect — are highlighted in orange across both halves and
#' feed the region explorer (genes / SNPs / GO / regionplot / haplotype).
#'
#' @param id module namespace id
#' @noRd
mod_overlapping_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # GEA filter bar
        bslib::card(
            class = "mb-0",
            bslib::card_header(
                htmltools::span("GEA Filter",
                    class = "fw-semibold text-primary me-2"),
                htmltools::span("(above zero line)",
                    class = "text-muted small")
            ),
            bslib::card_body(class = "p-2", shiny::uiOutput(ns("gea_filter_bar")))
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
                htmltools::span("GWAS Filter",
                    class = "fw-semibold text-success me-2"),
                htmltools::span("(below zero line)",
                    class = "text-muted small")
            ),
            bslib::card_body(class = "p-2", shiny::uiOutput(ns("gwas_filter_bar")))
        ),

        # Overlap bounds selector (persisted per project)
        bslib::card(
            class = "mt-2 mb-3",
            bslib::card_body(
                class = "py-2 px-3 d-flex align-items-center gap-4",
                htmltools::span("Overlap region bounds:", class = "fw-semibold me-1"),
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
mod_overlapping_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        module <- MOD_OVERLAP

        # ── GEA side data ──────────────────────────────────────────────────────
        gea_methods <- shiny::reactive(find_assoc_methods(project_data()$name, MOD_ASSOC))
        gea_traits  <- shiny::reactive(find_assoc_traits(project_data()$name,  MOD_ASSOC))

        gea_method_sigsnps <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("all_sigsnps_", pd$name, "_", MOD_ASSOC),
                        function() load_all_method_sigsnps(pd$name, MOD_ASSOC, k = pd$k_best))
        })

        gea_trait_colors  <- shiny::reactive(trait_color_map(gea_traits()))
        gea_method_shapes <- shiny::reactive(method_shape_map(gea_methods()))

        gea_combo_counts <- shiny::reactive({
            snps <- gea_method_sigsnps()
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
        gwas_methods <- shiny::reactive(find_assoc_methods(project_data()$name, MOD_PHENO))
        gwas_traits  <- shiny::reactive(find_assoc_traits(project_data()$name,  MOD_PHENO))

        gwas_method_sigsnps <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("all_sigsnps_", pd$name, "_", MOD_PHENO),
                        function() load_all_method_sigsnps(pd$name, MOD_PHENO, k = pd$k_best))
        })

        gwas_trait_colors  <- shiny::reactive({
            # GWAS uses the same color map; traits may overlap with GEA
            trait_color_map(gwas_traits())
        })
        gwas_method_shapes <- shiny::reactive(method_shape_map(gwas_methods()))

        gwas_combo_counts <- shiny::reactive({
            snps <- gwas_method_sigsnps()
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

        # ── Initialize per-side filter params from config/region_params ─────────
        # Each init_param: observe project_data, read saved param, fall back to config default
        .init_filter_param <- function(input_id, rp_key, config_default_fn) {
            shiny::observe({
                pd  <- project_data()
                rp  <- read_region_params(pd$name)
                val <- get_global_param(rp, module, rp_key)
                if (!is.null(val)) {
                    shiny::updateNumericInput(session, input_id, value = as.integer(val))
                } else {
                    shiny::updateNumericInput(session, input_id, value = config_default_fn(pd$config))
                }
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

        .init_filter_param("gea_region_distance", "gea_region_distance", function(cfg)
            config_get(cfg, "association", "region_distance", default = 1000000L))
        .init_filter_param("gea_combine_gap", "gea_combine_gap", function(cfg)
            config_get(cfg, "association", "combine_gap", default = 100000L))
        .init_filter_param("gwas_region_distance", "gwas_region_distance", function(cfg)
            config_get(cfg, "phenotype_association", "region_distance",
                       default = config_get(cfg, "association", "region_distance", default = 1000000L)))
        .init_filter_param("gwas_combine_gap", "gwas_combine_gap", function(cfg)
            config_get(cfg, "phenotype_association", "combine_gap",
                       default = config_get(cfg, "association", "combine_gap", default = 100000L)))

        gea_region_distance <- shiny::reactive({
            v <- input$gea_region_distance
            if (is.null(v) || is.na(v) || v < 1000L) 1000000L else as.integer(v)
        })
        gea_combine_gap <- shiny::reactive({
            v <- input$gea_combine_gap
            if (is.null(v) || is.na(v) || v < 0L) 100000L else as.integer(v)
        })
        gwas_region_distance <- shiny::reactive({
            v <- input$gwas_region_distance
            if (is.null(v) || is.na(v) || v < 1000L) 1000000L else as.integer(v)
        })
        gwas_combine_gap <- shiny::reactive({
            v <- input$gwas_combine_gap
            if (is.null(v) || is.na(v) || v < 0L) 100000L else as.integer(v)
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

        # ── GEA filter bar UI ──────────────────────────────────────────────────
        output$gea_filter_bar <- shiny::renderUI({
            build_filter_bar_ui(
                ns                     = ns,
                traits                 = gea_traits(),
                methods                = gea_methods(),
                trait_colors           = gea_trait_colors(),
                combo_counts           = gea_combo_counts(),
                default_strategy_value = "All",
                region_distance_value  = gea_region_distance(),
                combine_gap_value      = gea_combine_gap(),
                input_prefix           = "gea_"
            )
        })

        # ── GWAS filter bar UI ─────────────────────────────────────────────────
        output$gwas_filter_bar <- shiny::renderUI({
            build_filter_bar_ui(
                ns                     = ns,
                traits                 = gwas_traits(),
                methods                = gwas_methods(),
                trait_colors           = gwas_trait_colors(),
                combo_counts           = gwas_combo_counts(),
                default_strategy_value = "All",
                region_distance_value  = gwas_region_distance(),
                combine_gap_value      = gwas_combine_gap(),
                input_prefix           = "gwas_"
            )
        })

        # ── Interactive sig SNPs (independent per side) ────────────────────────
        gea_interactive_sigsnps <- shiny::reactive({
            compute_interactive_sigsnps(
                all_method_sigsnps = gea_method_sigsnps(),
                tm_selection_json  = input$gea_tm_selection,
                combo_counts       = gea_combo_counts(),
                known_traits       = gea_traits(),
                strategy           = input$gea_combine_strategy %||% "All",
                gap                = gea_combine_gap(),
                project_name       = project_data()$name,
                module             = MOD_ASSOC
            )
        })

        gwas_interactive_sigsnps <- shiny::reactive({
            compute_interactive_sigsnps(
                all_method_sigsnps = gwas_method_sigsnps(),
                tm_selection_json  = input$gwas_tm_selection,
                combo_counts       = gwas_combo_counts(),
                known_traits       = gwas_traits(),
                strategy           = input$gwas_combine_strategy %||% "All",
                gap                = gwas_combine_gap(),
                project_name       = project_data()$name,
                module             = MOD_PHENO
            )
        })

        # ── Independent regions (one set per source) ───────────────────────────
        gea_regions <- shiny::reactive({
            snps <- gea_interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(.empty_regions())
            compute_all_regions(snps, gea_region_distance())
        })

        gwas_regions <- shiny::reactive({
            snps <- gwas_interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(.empty_regions())
            compute_all_regions(snps, gwas_region_distance())
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
            override_regions    = overlap_regions
        )

        # ── Miami plot ─────────────────────────────────────────────────────────
        # Combine trait/method maps for both sources (unified legend)
        all_traits  <- shiny::reactive(union(gea_traits(),  gwas_traits()))
        all_methods <- shiny::reactive(union(gea_methods(), gwas_methods()))
        all_trait_colors  <- shiny::reactive(trait_color_map(all_traits()))
        all_method_shapes <- shiny::reactive(method_shape_map(all_methods()))

        miami_click <- mod_manhattan_overlay_server("miami",
            project_data         = project_data,
            module               = module,
            is_miami             = TRUE,
            combined             = TRUE,
            title_label          = shiny::reactive("Miami Plot (GEA \u2191 | GWAS \u2193)"),
            regions              = gea_regions,
            current_region_id    = explorer$selected_region_id,
            show_regions_control = FALSE,
            sig_snps_override    = unified_sigsnps,
            trait_colors         = all_trait_colors,
            method_shapes        = all_method_shapes,
            gwas_regions         = gwas_regions,
            overlap_regions      = overlap_regions
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

        # ── Pairwise Trait Overlap (unchanged, pipeline-computed) ──────────────
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
