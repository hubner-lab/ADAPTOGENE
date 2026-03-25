#' Overlapping Regions tab UI
#'
#' Miami Manhattan overlay, overlap pairs table, and region detail panel.
#'
#' @param id module namespace id
#' @noRd
mod_overlapping_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Miami Manhattan (always visible)
        mod_manhattan_overlay_ui(ns("miami"), height = "420px"),

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
        shiny::uiOutput(ns("region_detail_ui")),

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

        overlap_data <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("overlap_summary_", pd$name), function() {
                load_overlap_summary(pd$name)
            })
        })

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

        # ── Miami Manhattan ────────────────────────────────────────────────────
        miami_click <- mod_manhattan_overlay_server("miami",
            project_data         = project_data,
            module               = module,
            is_miami             = TRUE,
            combined             = TRUE,
            title_label          = shiny::reactive("Miami Plot (GEA \u2191 | GWAS \u2193)"),
            regions              = regions_data,
            show_regions_control = TRUE
        )

        shiny::observeEvent(miami_click(), {
            rid <- miami_click()
            if (!is.null(rid) && nzchar(rid)) {
                selected_region(rid)
                shiny::updateSelectInput(session, "region_id", selected = rid)
            }
        }, ignoreNULL = TRUE)

        # ── Region dropdown ────────────────────────────────────────────────────
        output$region_selector <- shiny::renderUI({
            rg <- regions_data()
            if (nrow(rg) == 0) {
                return(htmltools::p("No overlapping regions found.", class = "text-muted small"))
            }
            genes   <- genes_data()
            labels  <- build_region_labels(rg, genes)
            current <- selected_region()
            shiny::selectInput(ns("region_id"), "Select region",
                choices  = c(setNames("", ""), labels),
                selected = current %||% ""
            )
        })

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

        # ── Region detail panel ────────────────────────────────────────────────
        output$region_detail_ui <- shiny::renderUI({
            rid <- selected_region()
            if (is.null(rid)) return(NULL)
            mod_region_detail_ui(ns("region_detail"))
        })

        region_trait <- shiny::reactive({
            rid <- selected_region()
            rg  <- regions_data()
            if (is.null(rid) || nrow(rg) == 0) return(NULL)
            row <- rg[rg$region_id == rid, ]
            if (nrow(row) == 0) return(NULL)
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
