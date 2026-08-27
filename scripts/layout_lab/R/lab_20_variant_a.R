# lab_20_variant_a.R — Variant A: COMMAND CENTER
#
# At-a-glance operational health. Compact KPI strip with traffic-light
# amplification, ONE dominant hero figure, everything else demoted to grouped
# small multiples. Depth + het table live in a COLLAPSED accordion, which the
# Step 0 smoke test proved means the PNG is never transmitted until opened.
#
# Reuses mod_processing_server() unchanged: every id below is byte-identical.

#' Marks the one dominant figure so SCSS can enlarge it past the 45vh default.
#' Variant-private (A and C want it; B treats all figures equally).
lab_hero <- function(...) htmltools::div(class = "lab-hero", ...)

mod_processing_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "processing",

        lab_kpi_row(ns),

        htmltools::div(
            class = "lab-overview",
            bslib::layout_columns(
                col_widths = c(5, 7),
                # layout_columns defaults to fill=TRUE, which would stretch the
                # columns and fight the fraction-based heights below.
                fill = FALSE, fillable = FALSE,
                gap = "0.75rem",

                # LEFT — hero figure over the short filtering table
                htmltools::div(
                    class = "lab-hero-col",
                    lab_hero(mod_image_card_ui(ns("attrition"))),
                    lab_filtering_summary_card(ns)
                ),

                # RIGHT — grouped small multiples
                htmltools::div(
                    class = "lab-multiples-col",

                    lab_section_header("Sample QC", icon = "people-fill"),
                    lab_thumb_grid(
                        cols = 4,
                        lab_thumb(mod_image_card_ui(ns("sample_miss"))),
                        lab_thumb(mod_image_card_ui(ns("het_miss"))),
                        lab_thumb(mod_image_card_ui(ns("relatedness"))),
                        lab_thumb(mod_image_card_ui(ns("relatedness_mds")))
                    ),

                    lab_section_header("SNP QC", icon = "database-fill"),
                    # Same 4-column rhythm as Sample QC so every tile is the same
                    # width and the two groups read as one grid. A 2-col grid made
                    # these tiles ~520px wide while the height cap held the image
                    # to ~280px, stranding it in empty space.
                    lab_thumb_grid(
                        cols = 4,
                        lab_thumb(mod_image_card_ui(ns("snp_miss"))),
                        lab_thumb(mod_image_card_ui(ns("maf"))),
                        # a chromosome strip chart is unreadable in a square tile
                        lab_thumb(mod_image_card_ui(ns("snp_density")), span = 2)
                    )
                )
            )
        ),

        # ── Below the fold: both collapsed -> outputs suspended ───────────────
        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE, multiple = TRUE,
            bslib::accordion_panel(
                value = "depth",
                title = htmltools::tagList(
                    bsicons::bs_icon("layers"), " Depth Distribution ",
                    htmltools::span(class = "text-muted small",
                                    "(context — often empty for GBS)")
                ),
                shiny::uiOutput(ns("depth_section"))
            ),
            bslib::accordion_panel(
                value = "het",
                title = htmltools::tagList(bsicons::bs_icon("table"),
                                           " Sample Heterozygosity"),
                htmltools::div(
                    class = "d-flex justify-content-end mb-2",
                    shiny::downloadButton(ns("dl_het"), "CSV",
                                          class = "btn-sm btn-outline-secondary")
                ),
                DT::DTOutput(ns("het_table"))
            )
        )
    )
}
