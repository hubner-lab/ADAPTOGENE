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
        hap_scan_cache   <- shiny::reactiveValues()  # keyed by region_id___tag
        hap_viz_cache    <- shiny::reactiveValues()  # keyed by region_id___tag

        # Handles for background processes: keyed by cache key
        # Each entry: list(process, log_file, type, ...) or NULL when done
        subprocess_handles <- shiny::reactiveValues()

        # ── Haplotype tag for this module ──────────────────────────────────────
        matching_hap_tag <- shiny::reactive({
            pd <- project_data()
            find_matching_hap_tag(pd$name, module)
        })

        # expected_hap_tag: returns existing matching tag when found, otherwise
        # constructs one from config default so hap_key is stable before first scan.
        expected_hap_tag <- shiny::reactive({
            tag <- matching_hap_tag()
            if (!is.null(tag)) return(tag)
            src <- .module_to_hap_source(module)
            if (is.null(src)) return(NULL)
            pd   <- project_data()
            meta <- config_get(pd$config, "haplotype", "scan", "metadata_type", default = "site")
            paste0(meta, "_", src)
        })

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

        # combo_hash: deterministic 8-char hash of sorted trait/method pairs in region.
        # Different filter state → different hash → different cache key and disk path.
        combo_hash <- shiny::reactive({
            snps <- region_snps()
            if (is.null(snps) || nrow(snps) == 0) return(NULL)
            arg <- .build_trait_method_arg(snps)
            if (is.null(arg)) return(NULL)
            substr(digest::digest(arg, algo = "md5"), 1, 8)
        })

        rplot_key <- shiny::reactive({
            rid  <- selected_region_id()
            hash <- combo_hash()
            if (is.null(rid) || is.null(hash)) return(NULL)
            paste0(rid, "___", hash, "___rplot")
        })

        hap_key <- shiny::reactive({
            rid <- selected_region_id()
            tag <- expected_hap_tag()
            if (is.null(rid) || is.null(tag)) return(NULL)
            paste0(rid, "___", tag)
        })

        # ── Disk-cache hydration: load persisted regionplot on region/filter change ─
        shiny::observe({
            key <- rplot_key()
            if (is.null(key) || !is.null(regionplot_cache[[key]]) || is_running(key)) return()
            rid  <- selected_region_id()
            hash <- combo_hash()
            pd   <- project_data()
            if (is.null(rid) || is.null(hash) || is.null(pd)) return()
            region_safe <- gsub("[:\\-]", "_", rid)
            disk_dir  <- ondemand_regionplot_dir(pd$name, module, hash)
            disk_path <- file.path(disk_dir,
                                   paste0("regionplot_custom_", region_safe, ".png"))
            if (file_ok(disk_path)) {
                regionplot_cache[[key]] <- list(path = disk_path)
            }
        })

        # ── Disk-cache hydration: load persisted enrichment on trait/region change ─
        shiny::observe({
            key <- enrich_key()
            if (is.null(key) || !is.null(enrichment_cache[[key]]) || is_running(key)) return()
            rid <- selected_region_id()
            tr  <- current_trait()
            pd  <- project_data()
            if (is.null(rid) || is.null(tr) || is.null(pd)) return()
            tsv <- enrichment_table_path(pd$name, module, tr, rid)
            if (!file_ok(tsv)) return()
            tbl <- tryCatch(
                data.table::fread(tsv, sep = "\t", header = TRUE),
                error = function(e) NULL
            )
            if (is.null(tbl) || nrow(tbl) == 0) return()
            dotplot_p  <- enrichment_plot_path(pd$name, module, tr, rid, "dotplot")
            emapplot_p <- enrichment_plot_path(pd$name, module, tr, rid, "emapplot")
            enrichment_cache[[key]] <- list(
                table         = tbl,
                dotplot_path  = if (file_ok(dotplot_p))  dotplot_p  else NULL,
                emapplot_path = if (file_ok(emapplot_p)) emapplot_p else NULL
            )
        })

        # ── Disk-cache hydration: load pre-computed haplotype scan results ────
        shiny::observe({
            key <- hap_key()
            if (is.null(key) || !is.null(hap_scan_cache[[key]]) || is_running(paste0(key, "__scan"))) return()
            rid <- selected_region_id()
            tag <- expected_hap_tag()
            pd  <- project_data()
            if (is.null(rid) || is.null(tag) || is.null(pd)) return()
            if (check_hap_scan_results(pd$name, tag, rid)) {
                mg_p  <- hap_clustree_path(pd$name, tag, rid, "MG")
                hap_p <- hap_clustree_path(pd$name, tag, rid, "hap")
                hap_scan_cache[[key]] <- list(
                    mg_path  = if (file_ok(mg_p))  mg_p  else NULL,
                    hap_path = if (file_ok(hap_p)) hap_p else NULL
                )
            }
        })

        # ── Disk-cache hydration: load pre-computed haplotype viz results ────
        shiny::observe({
            key <- hap_key()
            if (is.null(key) || !is.null(hap_viz_cache[[key]]) || is_running(paste0(key, "__viz"))) return()
            rid <- selected_region_id()
            tag <- expected_hap_tag()
            pd  <- project_data()
            if (is.null(rid) || is.null(tag) || is.null(pd)) return()
            traits <- check_hap_viz_results(pd$name, tag, rid)
            if (length(traits) > 0) {
                hap_viz_cache[[key]] <- list(traits = traits)
            }
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

                } else if (h$type == "hap_scan") {
                    result <- tryCatch(
                        .collect_hap_scan_result(h),
                        error = function(e) list(error = conditionMessage(e))
                    )
                    if (!is.null(result$error) && status != 0L) {
                        log_lines <- tryCatch(readLines(h$log_file, warn = FALSE), error = function(e) character())
                        result$log_tail <- paste(utils::tail(log_lines, 10), collapse = "\n")
                    }
                    # key for hap_scan is stored with __scan suffix
                    base_key <- sub("__scan$", "", key)
                    hap_scan_cache[[base_key]] <- result

                } else if (h$type == "hap_viz") {
                    result <- tryCatch(
                        .collect_hap_viz_result(h),
                        error = function(e) list(error = conditionMessage(e))
                    )
                    if (!is.null(result$error) && status != 0L) {
                        log_lines <- tryCatch(readLines(h$log_file, warn = FALSE), error = function(e) character())
                        result$log_tail <- paste(utils::tail(log_lines, 10), collapse = "\n")
                    }
                    base_key <- sub("__viz$", "", key)
                    hap_viz_cache[[base_key]] <- result
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
                    open     = c("genes", "enrichment", "regionplot"),
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
                    ),

                    bslib::accordion_panel(
                        "Haplotype Analysis",
                        value = "haplotype",
                        icon  = bsicons::bs_icon("diagram-3"),
                        shiny::uiOutput(ns("haplotype_ui"))
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
                region_genes(), row$region_id[1], tr, project_data(), module
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

        # ── Haplotype UI ───────────────────────────────────────────────────────
        output$haplotype_ui <- shiny::renderUI({
            rid <- selected_region_id()
            if (is.null(rid)) return(NULL)
            key      <- hap_key()
            tag      <- expected_hap_tag()
            scan_res <- if (!is.null(key)) hap_scan_cache[[key]] else NULL
            viz_res  <- if (!is.null(key)) hap_viz_cache[[key]]  else NULL
            scan_run <- is_running(paste0(key %||% "", "__scan"))
            viz_run  <- is_running(paste0(key %||% "", "__viz"))

            # Config defaults for parameter inputs
            pd <- project_data()
            cfg_eps   <- config_get(pd$config, "haplotype", "epsilon_selected", default = NULL)
            cfg_range <- paste(config_get(pd$config, "haplotype", "scan",
                                          "epsilon_range", default = "0.3,0.5,0.7,0.9"),
                               collapse = ",")
            if (is.list(cfg_range)) cfg_range <- paste(unlist(cfg_range), collapse = ",")
            cfg_mgmin  <- config_get(pd$config, "haplotype", "scan", "min_group_size",  default = 50L)
            cfg_minhap <- config_get(pd$config, "haplotype", "scan", "min_haplotype_size", default = 15L)
            cfg_minsnp <- config_get(pd$config, "haplotype", "scan", "min_snps", default = 3L)
            cfg_meta   <- config_get(pd$config, "haplotype", "scan", "metadata_type", default = "site")

            # ── Sub-box 1: Scan ──────────────────────────────────────────────
            # Shared input form — used in both initial and error states
            scan_inputs <- htmltools::tagList(
                bslib::layout_column_wrap(
                    width = 1/2,
                    shiny::textInput(ns("hap_epsilon_range"), "Epsilon values",
                        value = cfg_range, placeholder = "e.g. 0.3,0.5,0.7,0.9"),
                    shiny::numericInput(ns("hap_mgmin"), "Min group size (MGmin)",
                        value = as.integer(cfg_mgmin), min = 1L, step = 1L)
                ),
                bslib::layout_column_wrap(
                    width = 1/2,
                    shiny::numericInput(ns("hap_minhap"), "Min haplotype size",
                        value = as.integer(cfg_minhap), min = 1L, step = 1L),
                    shiny::numericInput(ns("hap_min_snps"), "Min SNPs in region",
                        value = as.integer(cfg_minsnp), min = 1L, step = 1L)
                ),
                shiny::selectInput(ns("hap_meta_type"), "Metadata grouping",
                    choices  = c("site", paste0("cluster_K", pd$k_best)),
                    selected = cfg_meta, width = "220px"),
                htmltools::tags$details(
                    htmltools::tags$summary(
                        class = "text-muted small mb-1",
                        bsicons::bs_icon("info-circle"), " How to choose epsilon"
                    ),
                    htmltools::p(class = "text-muted small ms-3",
                        "Run the scan to explore haplotype clustering across epsilon values. ",
                        "Review the clustree plot \u2014 look for stable clusters where lines ",
                        "don\u2019t cross. Fewer crossings = more stable partitioning. ",
                        "Then enter your chosen epsilon in the Visualization box below."
                    )
                )
            )

            scan_content <- if (scan_run) {
                htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2 py-1",
                    bsicons::bs_icon("arrow-repeat",
                        class = "text-warning spin-icon flex-shrink-0", size = "0.85em"),
                    htmltools::span(class = "small", "Running haplotype scan...")
                )
            } else if (!is.null(scan_res) && is.null(scan_res$error)) {
                htmltools::tagList(
                    if (!is.null(scan_res$traits) && length(scan_res$traits) > 0) {
                        shiny::selectInput(ns("hap_scan_trait"), "Trait (clustree)",
                            choices = scan_res$traits, width = "220px")
                    },
                    shiny::uiOutput(ns("hap_clustree_ui"))
                )
            } else if (!is.null(scan_res) && !is.null(scan_res$error)) {
                htmltools::tagList(
                    htmltools::div(
                        class = "alert alert-warning small mb-3",
                        bsicons::bs_icon("exclamation-triangle"), " ", scan_res$error
                    ),
                    scan_inputs,
                    shiny::actionButton(ns("run_hap_scan"), "Retry Scan",
                        class = "btn-sm btn-outline-secondary",
                        icon  = shiny::icon("rotate-right"))
                )
            } else {
                htmltools::tagList(
                    scan_inputs,
                    shiny::actionButton(ns("run_hap_scan"), "Run Scan",
                        class = "btn-sm btn-outline-primary",
                        icon  = shiny::icon("magnifying-glass"))
                )
            }

            # ── Sub-box 2: Visualization ────────────────────────────────────
            viz_content <- if (viz_run) {
                htmltools::div(
                    class = "pipeline-rule-item d-flex align-items-center gap-2 py-1",
                    bsicons::bs_icon("arrow-repeat",
                        class = "text-warning spin-icon flex-shrink-0", size = "0.85em"),
                    htmltools::span(class = "small", "Running haplotype visualization...")
                )
            } else if (!is.null(viz_res) && is.null(viz_res$error)) {
                shiny::uiOutput(ns("hap_viz_plots_ui"))
            } else if (!is.null(viz_res) && !is.null(viz_res$error)) {
                htmltools::div(
                    class = "text-muted small py-2",
                    bsicons::bs_icon("exclamation-triangle"), " ", viz_res$error,
                    htmltools::br(),
                    shiny::actionButton(ns("run_hap_viz"), "Retry Visualization",
                        class = "btn-sm btn-outline-secondary mt-2",
                        icon  = shiny::icon("rotate-right"))
                )
            } else {
                htmltools::tagList(
                    shiny::numericInput(ns("hap_epsilon_selected"), "Epsilon (from clustree)",
                        value = if (!is.null(cfg_eps) && !is.na(as.numeric(cfg_eps)))
                                    as.numeric(cfg_eps) else NA_real_,
                        min = 0.01, max = 1, step = 0.1, width = "160px"),
                    htmltools::p(class = "text-muted small",
                        bsicons::bs_icon("arrow-up"), " Pick a stable epsilon from the clustree above."),
                    shiny::actionButton(ns("run_hap_viz"), "Run Visualization",
                        class = "btn-sm btn-outline-primary",
                        icon  = shiny::icon("chart-bar"))
                )
            }

            bslib::layout_column_wrap(
                width = 1/2,
                bslib::card(
                    bslib::card_header(
                        bsicons::bs_icon("search"), " Step 1: Scan"
                    ),
                    bslib::card_body(scan_content)
                ),
                bslib::card(
                    bslib::card_header(
                        bsicons::bs_icon("palette"), " Step 2: Visualization"
                    ),
                    bslib::card_body(viz_content)
                )
            )
        })

        # Clustree images UI container (rendered inside scan success state)
        output$hap_clustree_ui <- shiny::renderUI({
            bslib::layout_column_wrap(
                width = 1/2,
                mod_image_card_ui(ns("hap_clustree_mg")),
                mod_image_card_ui(ns("hap_clustree_hap"))
            )
        })

        # Image servers for clustree plots (re-initialized when scan cache or trait updates)
        shiny::observe({
            key      <- hap_key()
            scan_res <- if (!is.null(key)) hap_scan_cache[[key]] else NULL
            if (is.null(scan_res) || !is.null(scan_res$error)) return()

            # Pick per-trait or default paths
            trait <- input$hap_scan_trait
            use_trait <- !is.null(trait) && length(scan_res$traits) > 0 &&
                         trait %in% scan_res$traits
            if (use_trait) {
                rid   <- scan_res$region_id
                pdir  <- scan_res$plots_dir
                mg_p  <- file.path(pdir, paste0("Region_", rid, "_clustree_MG_",  trait, ".png"))
                hap_p <- file.path(pdir, paste0("Region_", rid, "_clustree_hap_", trait, ".png"))
                lbl   <- paste0(" (", trait, ")")
            } else {
                mg_p  <- scan_res$mg_path
                hap_p <- scan_res$hap_path
                lbl   <- ""
            }

            if (!is.null(mg_p) && file_ok(mg_p))
                mod_image_card_server("hap_clustree_mg",
                    path    = shiny::reactive(mg_p),
                    title   = shiny::reactive(paste0("Clustree \u2014 Marker Groups", lbl)),
                    dl_name = shiny::reactive(paste0("clustree_MG", if (use_trait) paste0("_", trait) else "")))
            if (!is.null(hap_p) && file_ok(hap_p))
                mod_image_card_server("hap_clustree_hap",
                    path    = shiny::reactive(hap_p),
                    title   = shiny::reactive(paste0("Clustree \u2014 Haplotypes", lbl)),
                    dl_name = shiny::reactive(paste0("clustree_hap", if (use_trait) paste0("_", trait) else "")))
        })

        # Viz plots: rendered after viz completes (trait selector + 3 image cards)
        output$hap_viz_plots_ui <- shiny::renderUI({
            key     <- hap_key()
            viz_res <- if (!is.null(key)) hap_viz_cache[[key]] else NULL
            if (is.null(viz_res) || !is.null(viz_res$error)) return(NULL)
            traits  <- viz_res$traits
            if (length(traits) == 0) return(plot_placeholder("No traits found in visualization output"))

            htmltools::tagList(
                shiny::selectInput(ns("hap_viz_trait"), "Trait",
                    choices = traits, width = "220px"),
                mod_image_card_ui(ns("hap_crosshap_viz")),
                bslib::layout_column_wrap(
                    width = 1/2,
                    mod_image_card_ui(ns("hap_boxplot")),
                    mod_image_card_ui(ns("hap_piemap"))
                )
            )
        })

        shiny::observe({
            key     <- hap_key()
            viz_res <- if (!is.null(key)) hap_viz_cache[[key]] else NULL
            if (is.null(viz_res) || !is.null(viz_res$error)) return()
            rid   <- selected_region_id()
            tag   <- expected_hap_tag()
            pd    <- project_data()
            trait <- input$hap_viz_trait %||% viz_res$traits[1]
            if (is.null(trait) || is.null(rid) || is.null(tag) || is.null(pd)) return()

            viz_p <- crosshap_viz_path(pd$name, tag, rid, trait)
            box_p <- hap_boxplot_path(pd$name, tag, rid, trait)
            pie_p <- hap_piemap_path(pd$name, tag, rid, trait)

            if (file_ok(viz_p))
                mod_image_card_server("hap_crosshap_viz",
                    path    = shiny::reactive(viz_p),
                    title   = shiny::reactive("Haplotype Visualization"),
                    dl_name = shiny::reactive(paste0("crosshap_", rid, "_", trait)))
            if (file_ok(box_p))
                mod_image_card_server("hap_boxplot",
                    path    = shiny::reactive(box_p),
                    title   = shiny::reactive("Haplotype Phenotype Boxplot"),
                    dl_name = shiny::reactive(paste0("hap_boxplot_", rid, "_", trait)))
            if (file_ok(pie_p))
                mod_image_card_server("hap_piemap",
                    path    = shiny::reactive(pie_p),
                    title   = shiny::reactive("Haplotype Piemap"),
                    dl_name = shiny::reactive(paste0("hap_piemap_", rid, "_", trait)))
        })

        # ── Run Scan button ────────────────────────────────────────────────────
        shiny::observeEvent(input$run_hap_scan, {
            row <- current_region_row()
            if (is.null(row)) return()
            src <- .module_to_hap_source(module)
            if (is.null(src)) return()
            meta <- input$hap_meta_type %||% "site"
            tag  <- paste0(meta, "_", src)
            # Use a key consistent with expected_hap_tag (which also constructs from config)
            key  <- hap_key()
            if (is.null(key) || is_running(paste0(key, "__scan"))) return()

            params <- list(
                epsilon_range = input$hap_epsilon_range %||% "0.3,0.5,0.7,0.9",
                mgmin         = input$hap_mgmin  %||% 50L,
                minhap        = input$hap_minhap %||% 15L,
                min_snps      = input$hap_min_snps %||% 3L,
                meta_type     = meta
            )
            handle <- launch_hap_scan_subprocess(row, project_data(), tag, params)
            if (!is.null(handle$error)) {
                hap_scan_cache[[key]] <- handle
            } else {
                subprocess_handles[[paste0(key, "__scan")]] <- handle
            }
        })

        # ── Run Viz button ─────────────────────────────────────────────────────
        shiny::observeEvent(input$run_hap_viz, {
            row <- current_region_row()
            if (is.null(row)) return()
            src <- .module_to_hap_source(module)
            if (is.null(src)) return()
            meta <- input$hap_meta_type %||% "site"
            tag  <- paste0(meta, "_", src)
            key  <- hap_key()
            if (is.null(key) || is_running(paste0(key, "__viz"))) return()

            params <- list(
                epsilon_selected = input$hap_epsilon_selected,
                meta_type        = meta
            )
            handle <- launch_hap_viz_subprocess(row, project_data(), tag, params)
            if (!is.null(handle$error)) {
                hap_viz_cache[[key]] <- handle
            } else {
                subprocess_handles[[paste0(key, "__viz")]] <- handle
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
