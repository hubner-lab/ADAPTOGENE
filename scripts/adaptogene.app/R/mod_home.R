#' Home tab UI
#'
#' Pipeline overview dashboard: value boxes, module status, summary table.
#'
#' @param id module namespace id
#' @noRd
mod_home_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Value box row
        bslib::layout_column_wrap(
            width = "200px",
            fill  = FALSE,
            bslib::value_box(
                title = "Samples",
                value = shiny::textOutput(ns("n_samples")),
                theme = "primary",
                showcase = bsicons::bs_icon("people")
            ),
            bslib::value_box(
                title = "SNPs (filtered)",
                value = shiny::textOutput(ns("n_snps")),
                theme = "info",
                showcase = bsicons::bs_icon("database")
            ),
            bslib::value_box(
                title = "K best",
                value = shiny::textOutput(ns("k_best")),
                theme = "success",
                showcase = bsicons::bs_icon("diagram-3")
            ),
            bslib::value_box(
                title = "Assoc. methods",
                value = shiny::textOutput(ns("n_methods")),
                theme = "warning",
                showcase = bsicons::bs_icon("bar-chart")
            ),
            bslib::value_box(
                title = "Regions",
                value = shiny::textOutput(ns("n_regions")),
                theme = "danger",
                showcase = bsicons::bs_icon("geo-alt")
            )
        ),

        shiny::br(),

        bslib::layout_column_wrap(
            width = 1 / 2,

            # Module status card
            bslib::card(
                bslib::card_header(
                    bsicons::bs_icon("check-circle"),
                    " Pipeline Status"
                ),
                bslib::card_body(
                    shiny::uiOutput(ns("module_status"))
                )
            ),

            # Summary table card
            bslib::card(
                full_screen = TRUE,
                bslib::card_header(
                    class = "d-flex justify-content-between align-items-center",
                    htmltools::span(
                        bsicons::bs_icon("table"),
                        " Pipeline Summary"
                    ),
                    shiny::downloadButton(ns("dl_summary"), "CSV",
                                          class = "btn-sm btn-outline-secondary")
                ),
                bslib::card_body(
                    DT::DTOutput(ns("summary_table"))
                )
            )
        )
    )
}

#' Home tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_home_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {

        summary_data <- shiny::reactive({
            pd <- project_data()
            load_pipeline_summary(pd$name)
        })

        # ── Value boxes ────────────────────────────────────────────────────────
        summary_val <- function(step, metric, default = "—") {
            shiny::reactive({
                s <- summary_data()
                if (nrow(s) == 0) return(default)
                row <- s[s$step == step & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            })
        }

        output$n_samples <- shiny::renderText(summary_val("processing", "n_samples")())
        output$n_snps    <- shiny::renderText(summary_val("processing", "n_snps_filtered")())
        output$k_best    <- shiny::renderText({
            pd <- project_data()
            if (!is.na(pd$k_best)) as.character(pd$k_best) else "—"
        })
        output$n_methods <- shiny::renderText({
            pd <- project_data()
            methods <- find_assoc_methods(pd$name)
            if (length(methods) > 0) as.character(length(methods)) else "—"
        })
        output$n_regions <- shiny::renderText(
            summary_val("association", "n_regions_combined")()
        )

        # ── Module status ──────────────────────────────────────────────────────
        output$module_status <- shiny::renderUI({
            pd  <- project_data()
            st  <- check_module_status(pd$name)
            labels <- c(
                processing  = "Processing",
                structure   = "Structure",
                structure_k = "Structure K",
                association = "Association",
                phenotype   = "Phenotype Association",
                overlapping = "Overlapping Regions",
                maladaptation = "Maladaptation",
                haplotype   = "Haplotype Analysis"
            )
            rows <- lapply(names(labels), function(nm) {
                mode_status_row(labels[[nm]], isTRUE(st[[nm]]))
            })
            htmltools::tagList(rows)
        })

        # ── Summary table ──────────────────────────────────────────────────────
        output$summary_table <- DT::renderDataTable({
            s <- summary_data()
            safe_datatable(
                as.data.frame(s),
                options = list(
                    dom       = "frtip",
                    scrollX   = TRUE,
                    pageLength = 25
                ),
                rownames = FALSE
            )
        })

        output$dl_summary <- shiny::downloadHandler(
            filename = function() {
                paste0("pipeline_summary_", project_data()$name, ".csv")
            },
            content = function(file) {
                s <- summary_data()
                utils::write.csv(s, file, row.names = FALSE)
            }
        )
    })
}
