#' Processing QC tab UI
#'
#' Displays all VCF QC diagnostic plots and tables from mode=processing.
#'
#' @param id module namespace id
#' @noRd
mod_processing_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(

        # ── Summary cards ──────────────────────────────────────────────────────
        bslib::layout_column_wrap(
            width = 1 / 4,
            fill  = FALSE,
            # Samples: raw → filtered
            bslib::value_box(
                title    = "Samples",
                value    = shiny::uiOutput(ns("samples_display")),
                theme    = "primary",
                showcase = bsicons::bs_icon("people-fill")
            ),
            # SNPs: raw → filtered
            bslib::value_box(
                title    = "SNPs",
                value    = shiny::uiOutput(ns("snps_display")),
                theme    = "info",
                showcase = bsicons::bs_icon("database-fill")
            ),
            # SNPs LD-pruned
            bslib::value_box(
                title    = "SNPs (LD-pruned)",
                value    = shiny::textOutput(ns("snps_ld")),
                theme    = "success",
                showcase = bsicons::bs_icon("scissors")
            ),
            # Ti/Tv Ratio
            bslib::value_box(
                title    = "Ti/Tv Ratio",
                value    = shiny::textOutput(ns("titv")),
                theme    = "warning",
                showcase = bsicons::bs_icon("arrow-left-right")
            )
        ),

        shiny::br(),

        # ── Row 2: Filtering overview ──────────────────────────────────────────
        bslib::layout_columns(
            col_widths = c(8, 4),
            mod_image_card_ui(ns("attrition")),
            bslib::card(
                full_screen = TRUE,
                bslib::card_header(
                    bsicons::bs_icon("table"),
                    " Filtering Summary"
                ),
                bslib::card_body(
                    DT::DTOutput(ns("filtering_table"))
                )
            )
        ),

        shiny::br(),

        # ── Row 3: Sample QC ───────────────────────────────────────────────────
        bslib::layout_column_wrap(
            width = 1 / 2,
            mod_image_card_ui(ns("sample_miss")),
            mod_image_card_ui(ns("het_miss"))
        ),

        shiny::br(),

        # ── Row 4: SNP QC ──────────────────────────────────────────────────────
        bslib::layout_column_wrap(
            width = 1 / 2,
            mod_image_card_ui(ns("snp_miss")),
            mod_image_card_ui(ns("maf"))
        ),

        shiny::br(),

        # ── Row 5: SNP density (full width) ───────────────────────────────────
        mod_image_card_ui(ns("snp_density")),

        shiny::br(),

        # ── Row 6: Depth (conditional) ─────────────────────────────────────────
        shiny::uiOutput(ns("depth_section")),

        # ── Row 7: Sample heterozygosity table ────────────────────────────────
        bslib::card(
            full_screen = TRUE,
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                htmltools::span(
                    bsicons::bs_icon("table"),
                    " Sample Heterozygosity"
                ),
                shiny::downloadButton(ns("dl_het"), "CSV",
                                      class = "btn-sm btn-outline-secondary")
            ),
            bslib::card_body(
                DT::DTOutput(ns("het_table"))
            )
        )
    )
}

#' Processing QC tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_processing_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        summary_data <- shiny::reactive({
            pd <- project_data()
            load_pipeline_summary(pd$name)
        })

        summary_val <- function(metric_key, default = "—") {
            shiny::reactive({
                s <- summary_data()
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "processing" & s$metric == metric_key, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            })
        }

        # ── Value boxes ────────────────────────────────────────────────────────
        arrow_display <- function(raw_id, filtered_id) {
            shiny::renderUI({
                raw      <- summary_val(raw_id)()
                filtered <- summary_val(filtered_id)()
                htmltools::div(
                    style = "display:flex; align-items:baseline; gap:0.4rem;",
                    htmltools::span(style = "opacity:0.7;", raw),
                    htmltools::span(style = "opacity:0.6; font-size:0.9em;", "\u2192"),
                    htmltools::strong(filtered)
                )
            })
        }
        output$samples_display <- arrow_display("samples_total", "samples_after_filtering")
        output$snps_display    <- arrow_display("snps_raw", "snps_after_filtering")
        output$snps_ld          <- shiny::renderText(summary_val("snps_after_ld_pruning")())
        output$titv             <- shiny::renderText(summary_val("titv_ratio")())

        # ── Image cards ────────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()

            mod_image_card_server("attrition",
                path    = shiny::reactive(qc_plot_path(pd$name, "filtering_attrition.png")),
                title   = shiny::reactive("Filtering Attrition"),
                dl_name = shiny::reactive(paste0(pd$name, "_filtering_attrition")),
                suggestion = shiny::reactive("Run mode=processing to generate QC plots")
            )
            mod_image_card_server("sample_miss",
                path    = shiny::reactive(qc_plot_path(pd$name, "sample_missingness_distribution.png")),
                title   = shiny::reactive("Sample Missingness"),
                dl_name = shiny::reactive(paste0(pd$name, "_sample_missingness"))
            )
            mod_image_card_server("het_miss",
                path    = shiny::reactive(qc_plot_path(pd$name, "het_vs_missingness.png")),
                title   = shiny::reactive("Heterozygosity vs Missingness"),
                dl_name = shiny::reactive(paste0(pd$name, "_het_vs_missingness"))
            )
            mod_image_card_server("snp_miss",
                path    = shiny::reactive(qc_plot_path(pd$name, "snp_missingness_distribution.png")),
                title   = shiny::reactive("SNP Missingness"),
                dl_name = shiny::reactive(paste0(pd$name, "_snp_missingness"))
            )
            mod_image_card_server("maf",
                path    = shiny::reactive(qc_plot_path(pd$name, "maf_distribution.png")),
                title   = shiny::reactive("MAF Distribution"),
                dl_name = shiny::reactive(paste0(pd$name, "_maf_distribution"))
            )
            mod_image_card_server("snp_density",
                path    = shiny::reactive(qc_plot_path(pd$name, "snp_density_by_chr.png")),
                title   = shiny::reactive("SNP Density by Chromosome"),
                dl_name = shiny::reactive(paste0(pd$name, "_snp_density"))
            )
        })

        # ── Filtering summary table ────────────────────────────────────────────
        output$filtering_table <- DT::renderDataTable({
            pd <- project_data()
            df <- as.data.frame(load_filtering_summary(pd$name))
            safe_datatable(
                df,
                options  = list(dom = "t", scrollX = TRUE),
                rownames = FALSE
            )
        })

        # ── Depth section (conditional) ────────────────────────────────────────
        depth_available <- shiny::reactive({
            s <- summary_data()
            if (nrow(s) == 0) return(FALSE)
            row <- s[s$step == "processing" & s$metric == "depth_qc", ]
            nrow(row) > 0 && row$value[1] == "available"
        })

        output$depth_section <- shiny::renderUI({
            if (!depth_available()) return(NULL)
            htmltools::tagList(
                mod_image_card_ui(ns("depth")),
                shiny::br()
            )
        })

        shiny::observe({
            shiny::req(depth_available())
            pd <- project_data()
            mod_image_card_server("depth",
                path    = shiny::reactive(qc_plot_path(pd$name, "depth_distribution.png")),
                title   = shiny::reactive("Depth Distribution"),
                dl_name = shiny::reactive(paste0(pd$name, "_depth_distribution"))
            )
        })

        # ── Sample heterozygosity table ────────────────────────────────────────
        het_data <- shiny::reactive({
            pd <- project_data()
            as.data.frame(load_sample_heterozygosity(pd$name))
        })

        output$het_table <- DT::renderDataTable({
            safe_datatable(
                het_data(),
                options  = list(dom = "frtip", scrollX = TRUE, pageLength = 25),
                rownames = FALSE
            )
        })

        output$dl_het <- shiny::downloadHandler(
            filename = function() paste0("sample_heterozygosity_", project_data()$name, ".csv"),
            content  = function(file) utils::write.csv(het_data(), file, row.names = FALSE)
        )
    })
}
