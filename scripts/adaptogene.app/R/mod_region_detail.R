#' Region detail panel UI
#'
#' Shows enrichment plots, genes table, GO table, regionplot, and haplotype viz
#' for a selected region. Reused across Association, Phenotype, and Overlapping tabs.
#'
#' @param id module namespace id
#' @noRd
mod_region_detail_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        shiny::uiOutput(ns("region_info")),
        bslib::accordion(
            id       = ns("detail_sections"),
            open     = c("genes"),
            multiple = TRUE,

            bslib::accordion_panel(
                "Genes in Region",
                value = "genes",
                icon  = bsicons::bs_icon("list-ul"),
                DT::DTOutput(ns("genes_table"))
            ),

            bslib::accordion_panel(
                "GO Enrichment Plots",
                value = "enrichment_plots",
                icon  = bsicons::bs_icon("bar-chart"),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("dotplot")),
                    mod_image_card_ui(ns("emapplot"))
                )
            ),

            bslib::accordion_panel(
                "GO Enrichment Table",
                value = "go_table",
                icon  = bsicons::bs_icon("table"),
                DT::DTOutput(ns("enrichment_table"))
            ),

            bslib::accordion_panel(
                "Regional Manhattan (topr)",
                value = "regionplot",
                icon  = bsicons::bs_icon("graph-up"),
                mod_image_card_ui(ns("regionplot"))
            ),

        )
    )
}

#' Region detail panel server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param region_id reactive character: currently selected region_id (NULL = hidden)
#' @param module character: MOD_ASSOC, MOD_PHENO, or MOD_OVERLAP
#' @param trait reactive character: currently selected trait (for enrichment paths)
#' @param genes_data reactive data.table from load_genes()
#' @noRd
mod_region_detail_server <- function(id, project_data, region_id,
                                      module     = MOD_ASSOC,
                                      trait      = shiny::reactive(NULL),
                                      genes_data = shiny::reactive(data.table::data.table())) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Region info bar
        output$region_info <- shiny::renderUI({
            rid <- region_id()
            shiny::req(rid)
            genes <- genes_data()
            rg <- if (nrow(genes) > 0 && "region_id" %in% names(genes))
                genes[genes$region_id == rid, ] else data.table::data.table()
            n_exon <- if (nrow(rg) > 0) sum(rg$exon_snp_count,     na.rm = TRUE) else 0L
            n_prom <- if (nrow(rg) > 0) sum(rg$promoter_snp_count, na.rm = TRUE) else 0L
            region_info_bar(rid, n_exon = n_exon, n_promoter = n_prom)
        })

        # ── Enrichment images ──────────────────────────────────────────────────
        # Use the first available trait if trait() is NULL
        effective_trait <- shiny::reactive({
            t <- trait()
            if (!is.null(t) && nzchar(t)) return(t)
            # Try to find from genes table
            rid <- region_id()
            genes <- genes_data()
            if (!is.null(rid) && nrow(genes) > 0 && "trait" %in% names(genes)) {
                tr <- unique(genes[genes$region_id == rid, "trait"])
                if (length(tr) > 0) return(tr[[1]])
            }
            NULL
        })

        shiny::observe({
            rid <- region_id()
            shiny::req(rid)
            tr <- effective_trait()
            pd <- project_data()

            dotplot_p  <- if (!is.null(tr))
                enrichment_plot_path(pd$name, module, tr, rid, "dotplot")  else NULL
            emapplot_p <- if (!is.null(tr))
                enrichment_plot_path(pd$name, module, tr, rid, "emapplot") else NULL

            mod_image_card_server("dotplot",
                path        = shiny::reactive(dotplot_p),
                title       = shiny::reactive("GO Dotplot"),
                dl_name     = shiny::reactive(paste0("dotplot_", rid)),
                placeholder = shiny::reactive("No enriched GO terms found for this region")
            )
            mod_image_card_server("emapplot",
                path        = shiny::reactive(emapplot_p),
                title       = shiny::reactive("GO Emapplot"),
                dl_name     = shiny::reactive(paste0("emapplot_", rid)),
                placeholder = shiny::reactive("Emapplot requires \u22652 enriched terms (dotplot needs 1+)")
            )
        })

        # ── Genes table ────────────────────────────────────────────────────────
        output$genes_table <- DT::renderDataTable({
            rid   <- region_id()
            shiny::req(rid)
            genes <- genes_data()
            rg <- if (nrow(genes) > 0 && "region_id" %in% names(genes))
                genes[genes$region_id == rid, ] else data.table::data.table()
            safe_datatable(
                as.data.frame(rg),
                extensions = "Buttons",
                options    = list(
                    dom       = "Bfrtip",
                    buttons   = list("csv"),
                    scrollX   = TRUE,
                    pageLength = 10
                ),
                rownames = FALSE
            )
        })

        # ── GO enrichment table ────────────────────────────────────────────────
        output$enrichment_table <- DT::renderDataTable({
            rid <- region_id()
            shiny::req(rid)
            tr  <- effective_trait()
            pd  <- project_data()
            enrich <- if (!is.null(tr))
                load_enrichment_table(pd$name, module, tr, rid)
            else
                data.table::data.table()

            dt <- safe_datatable(
                as.data.frame(enrich),
                extensions = "Buttons",
                options    = list(
                    dom       = "Bfrtip",
                    buttons   = list("csv"),
                    scrollX   = TRUE,
                    pageLength = 10
                ),
                rownames = FALSE
            )
            dt_format_pvals(dt)
        })

        # ── Regional Manhattan (topr) ──────────────────────────────────────────
        shiny::observe({
            rid <- region_id()
            shiny::req(rid)
            tr  <- effective_trait()
            pd  <- project_data()
            rp  <- regionplot_path(pd$name, rid, tr)
            mod_image_card_server("regionplot",
                path     = shiny::reactive(rp),
                title    = shiny::reactive("Regional Manhattan"),
                dl_name  = shiny::reactive(paste0("regionplot_", rid))
            )
        })

    })
}
