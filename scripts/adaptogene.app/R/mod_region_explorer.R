#' Interactive Region Explorer UI
#'
#' Replaces the static region dropdown in the Association tab.
#' Regions are auto-computed from the current sig SNPs on every filter/strategy change.
#' Clicking a region row shows genes (on-the-fly), and on-demand enrichment / regionplot
#' via processx subprocesses.
#'
#' @param id module namespace id
#' @noRd
mod_region_explorer_ui <- function(id) {
    ns <- shiny::NS(id)

    htmltools::tagList(

        # ── Region list card ─────────────────────────────────────────────────
        bslib::card(
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                htmltools::span("Regions"),
                shiny::uiOutput(ns("region_count_badge"), inline = TRUE)
            ),
            bslib::card_body(
                class = "p-2",
                DT::DTOutput(ns("region_table"), width = "100%")
            )
        ),

        # ── Region detail (conditional on selection) ─────────────────────────
        shiny::uiOutput(ns("region_detail_ui"))
    )
}

#' Interactive Region Explorer Server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param module character: MOD_ASSOC or MOD_PHENO
#' @param interactive_sigsnps reactive data.table of current filtered sig SNPs
#' @return list with computed_regions (reactive) and selected_region_id (reactiveVal)
#' @noRd
mod_region_explorer_server <- function(id, project_data, module = MOD_ASSOC,
                                        interactive_sigsnps = shiny::reactive(NULL),
                                        region_distance     = shiny::reactive(2000000L)) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ── Computed regions ───────────────────────────────────────────────────
        computed_regions <- shiny::reactive({
            snps <- interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(.empty_regions())
            compute_all_regions(snps, region_distance())
        })

        # ── Selected region state ──────────────────────────────────────────────
        selected_region_id <- shiny::reactiveVal(NULL)

        shiny::observeEvent(computed_regions(), {
            rid  <- selected_region_id()
            regs <- computed_regions()
            if (!is.null(rid) && (nrow(regs) == 0 || !rid %in% regs$region_id)) {
                selected_region_id(NULL)
            }
        }, ignoreInit = TRUE)

        # ── GFF genes (cached per project) ─────────────────────────────────────
        gff_genes <- shiny::reactive({
            pd <- project_data()
            load_gff_genes(pd$name, pd$config)
        })

        # ── Region count badge ─────────────────────────────────────────────────
        output$region_count_badge <- shiny::renderUI({
            n   <- nrow(computed_regions())
            cls <- if (n == 0) "bg-secondary" else "bg-primary"
            htmltools::span(class = paste("badge rounded-pill", cls), n, "regions")
        })

        # ── Region list table ──────────────────────────────────────────────────
        output$region_table <- DT::renderDataTable({
            regs <- computed_regions()
            if (nrow(regs) == 0) {
                return(safe_datatable(
                    data.frame(Message = "No significant SNPs with current filters"),
                    options  = list(dom = "t"),
                    rownames = FALSE
                ))
            }

            rid <- selected_region_id()
            selected_row <- if (!is.null(rid)) which(regs$region_id == rid) else integer(0)

            disp <- data.frame(
                Region  = sapply(regs$region_id, format_region_id),
                Chr     = regs$chr,
                Length  = .format_bp(regs$length),
                SNPs    = regs$snp_count,
                Traits  = regs$traits,
                Methods = regs$methods,
                MinP    = signif(regs$min_pvalue, 3),
                stringsAsFactors = FALSE
            )

            safe_datatable(
                disp,
                selection  = list(mode = "single", selected = selected_row),
                options    = list(
                    dom        = "ftp",
                    pageLength = 10,
                    scrollX    = TRUE,
                    order      = list(list(3L, "desc"))
                ),
                rownames = FALSE
            )
        })

        shiny::observeEvent(input$region_table_rows_selected, {
            row  <- input$region_table_rows_selected
            regs <- computed_regions()
            if (length(row) == 1 && nrow(regs) >= row) {
                selected_region_id(regs$region_id[row])
            } else {
                selected_region_id(NULL)
            }
        }, ignoreNULL = FALSE)

        # ── Genes for selected region ──────────────────────────────────────────
        region_genes <- shiny::reactive({
            rid <- selected_region_id()
            shiny::req(rid)
            regs <- computed_regions()
            row  <- regs[region_id == rid]
            if (nrow(row) == 0) return(data.table::data.table())
            gff  <- gff_genes()
            if (is.null(gff) || nrow(gff) == 0) return(data.table::data.table())
            pd   <- project_data()
            plen <- config_get(pd$config, "association", "promoter_length", default = 10000L)
            find_genes_in_region(gff, row, as.integer(plen))
        })

        # ── SNPs for selected region ───────────────────────────────────────────
        region_snps <- shiny::reactive({
            rid  <- selected_region_id()
            shiny::req(rid)
            regs <- computed_regions()
            row  <- regs[region_id == rid]
            if (nrow(row) == 0) return(data.table::data.table())
            snps <- interactive_sigsnps()
            if (is.null(snps) || nrow(snps) == 0) return(data.table::data.table())
            snps[as.character(chr) == as.character(row$chr[1]) &
                 as.integer(pos) >= row$start[1] &
                 as.integer(pos) <= row$end[1]]
        })

        # ── Enrichment state ───────────────────────────────────────────────────
        enrichment_cache <- shiny::reactiveValues()
        regionplot_cache <- shiny::reactiveValues()

        # Handles for background processes: keyed by cache key
        # Each entry: list(process, log_file, type, ...) or NULL when done
        subprocess_handles <- shiny::reactiveValues()

        # ── Derived keys (region + selected trait) ────────────────────────────
        current_region_row <- shiny::reactive({
            rid  <- selected_region_id()
            shiny::req(rid)
            regs <- computed_regions()
            row  <- regs[region_id == rid]
            if (nrow(row) == 0) return(NULL)
            row
        })

        current_trait <- shiny::reactive({
            tr <- input$detail_trait
            if (is.null(tr) || tr == "") return(NULL)
            tr
        })

        # ── Update trait selector when region changes (smart default: lowest p) ─
        shiny::observe({
            row <- current_region_row()
            if (is.null(row)) return()
            traits <- trimws(strsplit(as.character(row$traits[1]), ",")[[1]])
            traits <- traits[nzchar(traits)]
            if (length(traits) == 0) return()

            # Pick default: trait with lowest min p-value among sig SNPs in region
            snps <- region_snps()
            default_trait <- if (!is.null(snps) && nrow(snps) > 0 && "trait" %in% names(snps)) {
                per_trait <- snps[, .(min_p = min(pvalue, na.rm = TRUE)), by = "trait"]
                per_trait <- per_trait[trait %in% traits]
                if (nrow(per_trait) > 0) per_trait[which.min(min_p), trait] else traits[1]
            } else {
                traits[1]
            }

            shiny::updateSelectInput(session, "detail_trait",
                choices  = traits,
                selected = default_trait)
        })

        enrich_key <- shiny::reactive({
            rid <- selected_region_id()
            tr  <- current_trait()
            if (is.null(rid) || is.null(tr)) return(NULL)
            paste0(rid, "___", tr)
        })

        rplot_key <- shiny::reactive({
            rid <- selected_region_id()
            if (is.null(rid)) return(NULL)
            paste0(rid, "___rplot")
        })

        # ── Helper: is a subprocess currently running for this key? ────────────
        is_running <- function(key) {
            if (is.null(key)) return(FALSE)
            h <- subprocess_handles[[key]]
            !is.null(h) && h$process$is_alive()
        }

        # ── Polling observer: checks background processes every second ──────────
        shiny::observe({
            handles_list <- shiny::reactiveValuesToList(subprocess_handles)
            active <- Filter(Negate(is.null), handles_list)
            if (length(active) == 0) return()

            shiny::invalidateLater(1000)

            for (key in names(active)) {
                h <- active[[key]]
                if (h$process$is_alive()) next  # still running

                status <- tryCatch(h$process$get_exit_status(), error = function(e) 1L)

                if (h$type == "enrichment") {
                    result <- tryCatch(
                        .collect_enrichment_result(h),
                        error = function(e) list(error = conditionMessage(e))
                    )
                    if (!is.null(result$error) && status != 0L) {
                        log_lines <- tryCatch(readLines(h$log_file, warn = FALSE), error = function(e) character())
                        result$log_tail <- paste(utils::tail(log_lines, 10), collapse = "\n")
                    }
                    enrichment_cache[[key]] <- result

                } else if (h$type == "regionplot") {
                    result <- tryCatch(
                        .collect_regionplot_result(h),
                        error = function(e) list(error = conditionMessage(e))
                    )
                    if (!is.null(result$error) && status != 0L) {
                        log_lines <- tryCatch(readLines(h$log_file, warn = FALSE), error = function(e) character())
                        result$log_tail <- paste(utils::tail(log_lines, 10), collapse = "\n")
                    }
                    regionplot_cache[[key]] <- result
                }

                subprocess_handles[[key]] <- NULL
            }
        })

        # ── Detail panel (conditional) ─────────────────────────────────────────
        output$region_detail_ui <- shiny::renderUI({
            rid <- selected_region_id()
            if (is.null(rid)) return(NULL)
            htmltools::tagList(
                shiny::uiOutput(ns("detail_info_bar")),
                htmltools::div(
                    class = "d-flex align-items-center gap-2 px-1 py-2",
                    htmltools::strong("Trait:", class = "text-nowrap small"),
                    shiny::selectInput(ns("detail_trait"), label = NULL,
                        choices = NULL, width = "220px")
                ),
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
                        "SNPs in Region",
                        value = "snps",
                        icon  = bsicons::bs_icon("dot"),
                        DT::DTOutput(ns("snps_table"))
                    ),

                    bslib::accordion_panel(
                        "GO Enrichment",
                        value = "enrichment",
                        icon  = bsicons::bs_icon("bar-chart"),
                        shiny::uiOutput(ns("enrichment_ui"))
                    ),

                    bslib::accordion_panel(
                        "Regional Manhattan",
                        value = "regionplot",
                        icon  = bsicons::bs_icon("graph-up"),
                        shiny::uiOutput(ns("regionplot_ui"))
                    )
                )
            )
        })

        # ── Detail info bar ────────────────────────────────────────────────────
        output$detail_info_bar <- shiny::renderUI({
            row <- current_region_row()
            if (is.null(row)) return(NULL)
            region_info_bar(row$region_id[1], n_snps = row$snp_count[1])
        })

        # ── Genes table ────────────────────────────────────────────────────────
        output$genes_table <- DT::renderDataTable({
            shiny::req(selected_region_id())
            genes <- region_genes()
            if (nrow(genes) == 0) {
                return(safe_datatable(
                    data.frame(Message = "No genes found (check GFF config)"),
                    options = list(dom = "t"), rownames = FALSE))
            }
            exclude_cols <- c("attributes", "source", "feature", "score", "strand", "phase")
            show_cols <- setdiff(names(genes), exclude_cols)
            safe_datatable(
                as.data.frame(genes[, show_cols, with = FALSE]),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                  scrollX = TRUE, pageLength = 10),
                rownames = FALSE
            )
        })

        # ── SNPs table ─────────────────────────────────────────────────────────
        output$snps_table <- DT::renderDataTable({
            shiny::req(selected_region_id())
            snps <- region_snps()
            if (nrow(snps) == 0) {
                return(safe_datatable(data.frame(Message = "No SNPs in region"),
                    options = list(dom = "t"), rownames = FALSE))
            }
            keep <- intersect(c("SNPID", "chr", "pos", "pvalue", "method", "trait"), names(snps))
            safe_datatable(
                as.data.frame(snps[, keep, with = FALSE]),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                  scrollX = TRUE, pageLength = 15),
                rownames = FALSE
            )
        })

        # ── Enrichment UI ──────────────────────────────────────────────────────
        output$enrichment_ui <- shiny::renderUI({
            rid <- selected_region_id()
            if (is.null(rid)) return(NULL)
            key    <- enrich_key()
            cached <- if (!is.null(key)) enrichment_cache[[key]] else NULL
            running <- is_running(key)

            if (running) {
                return(htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2 py-1",
                    bsicons::bs_icon("arrow-repeat",
                        class = "text-warning spin-icon flex-shrink-0", size = "0.85em"),
                    htmltools::span(class = "small", "Running GO enrichment...")
                ))
            }

            if (!is.null(cached)) {
                if (!is.null(cached$error)) {
                    return(htmltools::div(
                        class = "text-muted small py-2",
                        bsicons::bs_icon("exclamation-triangle"), " ", cached$error,
                        htmltools::br(),
                        shiny::actionButton(ns("run_enrichment"),
                            "Re-run Enrichment",
                            class = "btn-sm btn-outline-secondary mt-2",
                            icon  = shiny::icon("rotate-right"))
                    ))
                }
                return(htmltools::tagList(
                    bslib::layout_column_wrap(
                        width = 1 / 2,
                        if (!is.null(cached$dotplot_path) && file_ok(cached$dotplot_path))
                            mod_image_card_ui(ns("enrich_dotplot"))
                        else
                            bslib::card(bslib::card_body(
                                plot_placeholder("Dotplot not available"))),
                        if (!is.null(cached$emapplot_path) && file_ok(cached$emapplot_path))
                            mod_image_card_ui(ns("enrich_emapplot"))
                        else
                            bslib::card(bslib::card_body(
                                plot_placeholder("Emapplot not available (\u22652 GO terms needed)")))
                    ),
                    DT::DTOutput(ns("enrich_table"))
                ))
            }

            shiny::actionButton(ns("run_enrichment"),
                "Run Enrichment",
                class = "btn-sm btn-outline-primary",
                icon  = shiny::icon("play"))
        })

        # Image card servers (re-initialized when new enrichment result cached)
        shiny::observe({
            key    <- enrich_key()
            cached <- if (!is.null(key)) enrichment_cache[[key]] else NULL
            if (is.null(cached) || !is.null(cached$error)) return()
            mod_image_card_server("enrich_dotplot",
                path    = shiny::reactive(cached$dotplot_path),
                title   = shiny::reactive("GO Dotplot"),
                dl_name = shiny::reactive("enrichment_dotplot")
            )
            mod_image_card_server("enrich_emapplot",
                path    = shiny::reactive(cached$emapplot_path),
                title   = shiny::reactive("GO Emapplot"),
                dl_name = shiny::reactive("enrichment_emapplot")
            )
        })

        output$enrich_table <- DT::renderDataTable({
            key    <- enrich_key()
            cached <- if (!is.null(key)) enrichment_cache[[key]] else NULL
            if (is.null(cached) || !is.null(cached$error) || is.null(cached$table))
                return(safe_datatable(data.frame(), options = list(dom = "t"), rownames = FALSE))
            dt_format_pvals(safe_datatable(
                as.data.frame(cached$table),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                  scrollX = TRUE, pageLength = 10),
                rownames = FALSE
            ))
        })

        shiny::observeEvent(input$run_enrichment, {
            key <- enrich_key()
            if (is.null(key) || is_running(key)) return()
            row   <- current_region_row()
            if (is.null(row)) return()
            tr    <- current_trait()
            if (is.null(tr)) return()

            handle <- launch_enrichment_subprocess(
                region_genes(), row$region_id[1], tr, project_data()
            )
            if (!is.null(handle$error)) {
                enrichment_cache[[key]] <- handle   # store error immediately
            } else {
                subprocess_handles[[key]] <- handle # start polling
            }
        })

        # ── Regionplot UI ──────────────────────────────────────────────────────
        output$regionplot_ui <- shiny::renderUI({
            rid    <- selected_region_id()
            if (is.null(rid)) return(NULL)
            key    <- rplot_key()
            cached <- if (!is.null(key)) regionplot_cache[[key]] else NULL
            running <- is_running(key)

            if (running) {
                return(htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2 py-1",
                    bsicons::bs_icon("arrow-repeat",
                        class = "text-warning spin-icon flex-shrink-0", size = "0.85em"),
                    htmltools::span(class = "small", "Generating regional Manhattan...")
                ))
            }

            if (!is.null(cached)) {
                if (!is.null(cached$error)) {
                    return(htmltools::div(
                        class = "text-muted small py-2",
                        bsicons::bs_icon("exclamation-triangle"), " ", cached$error,
                        htmltools::br(),
                        shiny::actionButton(ns("run_regionplot"),
                            "Retry Region Plot",
                            class = "btn-sm btn-outline-secondary mt-2",
                            icon  = shiny::icon("rotate-right"))
                    ))
                }
                path <- cached$path
                if (!is.null(path) && file_ok(path)) {
                    return(mod_image_card_ui(ns("regionplot_img")))
                }
            }

            shiny::actionButton(ns("run_regionplot"),
                "Generate Region Plot",
                class = "btn-sm btn-outline-primary",
                icon  = shiny::icon("chart-line"))
        })

        shiny::observe({
            key    <- rplot_key()
            cached <- if (!is.null(key)) regionplot_cache[[key]] else NULL
            if (is.null(cached) || !is.null(cached$error)) return()
            path <- cached$path
            if (is.null(path) || !file_ok(path)) return()
            mod_image_card_server("regionplot_img",
                path    = shiny::reactive(path),
                title   = shiny::reactive("Regional Manhattan"),
                dl_name = shiny::reactive("regionplot")
            )
        })

        shiny::observeEvent(input$run_regionplot, {
            key <- rplot_key()
            if (is.null(key) || is_running(key)) return()
            row  <- current_region_row()
            if (is.null(row)) return()
            snps <- region_snps()
            if (is.null(snps) || nrow(snps) == 0) return()

            handle <- launch_regionplot_subprocess(row, snps, project_data(), module)
            if (!is.null(handle$error)) {
                regionplot_cache[[key]] <- handle   # store error immediately
            } else {
                subprocess_handles[[key]] <- handle # start polling
            }
        })

        # ── Return value for parent module ─────────────────────────────────────
        list(
            computed_regions   = computed_regions,
            selected_region_id = selected_region_id
        )
    })
}

# ── Helper ────────────────────────────────────────────────────────────────────

#' Format a bp distance as Mb/kb/bp string
#' @noRd
.format_bp <- function(bp) {
    bp <- as.numeric(bp)
    ifelse(bp >= 1e6, paste0(round(bp / 1e6, 1), " Mb"),
    ifelse(bp >= 1e3, paste0(round(bp / 1e3,  1), " kb"),
                      paste0(bp, " bp")))
}
