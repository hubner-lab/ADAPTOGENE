#' Processing QC tab UI
#'
#' Displays all VCF QC diagnostic plots and tables from mode=processing.
#'
#' @param id module namespace id
#' @noRd
mod_processing_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(

        # ── Summary cards (dynamic themes based on data quality) ──────────────
        shiny::uiOutput(ns("summary_boxes")),

        # ── LEA conversion alert (shown only if SNPs were dropped) ─────────────
        shiny::uiOutput(ns("lea_losses_alert")),

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
        # Relatedness caveat/counts live in the hover note badge next to the plot
        # title (mod_image_card's `note` slot) rather than a static caption below —
        # keeps the tab compact and matches the "minimal on-screen text" rule.
        bslib::layout_column_wrap(
            width = 1 / 3,
            mod_image_card_ui(ns("sample_miss")),
            mod_image_card_ui(ns("het_miss")),
            mod_image_card_ui(ns("relatedness"))
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

        # ── Value boxes (dynamic themes + conditional het/depth boxes) ─────────
        output$summary_boxes <- shiny::renderUI({
            pd  <- project_data()
            # Convert to data.frame to avoid data.table column-name shadowing:
            # inside data.table's [, bare names like `metric` resolve to the column
            # named `metric` in the table, not the function parameter.
            s <- as.data.frame(summary_data())

            sv <- function(metric, default = NA) {
                row <- s[s$step == "processing" & s$metric == metric, ]
                if (nrow(row) == 0) default else suppressWarnings(as.numeric(row$value[1]))
            }
            sv_str <- function(metric, default = "\u2014") {
                row <- s[s$step == "processing" & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            }

            # Arrow display helper
            arrow_val <- function(raw, filtered) {
                htmltools::div(
                    style = "display:flex; align-items:baseline; gap:0.4rem;",
                    htmltools::span(style = "opacity:0.7;", raw),
                    htmltools::span(style = "opacity:0.6; font-size:0.9em;", "\u2192"),
                    htmltools::strong(filtered)
                )
            }

            # Samples — warning if >10% removed
            n_samp_raw <- sv("samples_total")
            n_samp_flt <- sv("samples_after_filtering")
            samp_theme <- if (!is.na(n_samp_raw) && !is.na(n_samp_flt) && n_samp_raw > 0) {
                if ((n_samp_raw - n_samp_flt) / n_samp_raw > 0.1) "warning" else "success"
            } else "primary"
            samp_box <- bslib::value_box(
                title    = "Samples",
                value    = arrow_val(sv_str("samples_total"), sv_str("samples_after_filtering")),
                theme    = samp_theme,
                showcase = bsicons::bs_icon("people-fill")
            )

            # SNPs — warning if >50% removed
            n_snp_raw <- sv("snps_raw")
            n_snp_flt <- sv("snps_after_filtering")
            snp_theme <- if (!is.na(n_snp_raw) && !is.na(n_snp_flt) && n_snp_raw > 0) {
                if ((n_snp_raw - n_snp_flt) / n_snp_raw > 0.5) "warning" else "success"
            } else "info"
            snp_box <- bslib::value_box(
                title    = "SNPs",
                value    = arrow_val(sv_str("snps_raw"), sv_str("snps_after_filtering")),
                theme    = snp_theme,
                showcase = bsicons::bs_icon("database-fill")
            )

            # SNPs LD-pruned — always info
            ld_box <- bslib::value_box(
                title    = "SNPs (LD-pruned)",
                value    = sv_str("snps_after_ld_pruning"),
                theme    = "info",
                showcase = bsicons::bs_icon("scissors")
            )

            # Ti/Tv — color-coded by quality
            titv_val <- sv("titv_ratio")
            titv_theme <- if (is.na(titv_val)) "secondary"
                          else if (titv_val >= 2.0) "success"
                          else if (titv_val >= 1.5) "warning"
                          else "danger"
            titv_box <- bslib::value_box(
                title    = "Ti/Tv Ratio",
                value    = sv_str("titv_ratio"),
                theme    = titv_theme,
                showcase = bsicons::bs_icon("arrow-left-right")
            )

            # Heterozygosity outliers — conditional on config
            het_sd <- config_get(pd$config, "Filter", "het_outlier_sd", default = NULL)
            het_box <- if (!is.null(het_sd) && !is.na(het_sd)) {
                n_het <- sv("samples_het_outliers_removed", default = 0)
                het_theme <- if (is.na(n_het) || n_het == 0) "success" else "warning"
                n_het_str <- if (is.na(n_het)) "\u2014" else as.character(as.integer(n_het))
                bslib::value_box(
                    title    = "Het Outliers",
                    value    = paste0(n_het_str, " removed"),
                    theme    = het_theme,
                    showcase = bsicons::bs_icon("activity")
                )
            } else NULL

            # Depth filter — conditional on config
            min_dp <- config_get(pd$config, "Filter", "min_depth", default = NULL)
            max_dp <- config_get(pd$config, "Filter", "max_depth", default = NULL)
            depth_box <- if (!is.null(min_dp) || !is.null(max_dp)) {
                dp_parts <- c(
                    if (!is.null(min_dp)) paste0("min=", min_dp),
                    if (!is.null(max_dp)) paste0("max=", max_dp)
                )
                bslib::value_box(
                    title    = "Depth Filter",
                    value    = paste(dp_parts, collapse = ", "),
                    theme    = "info",
                    showcase = bsicons::bs_icon("layers")
                )
            } else NULL

            # Relatedness (IBD) removal — only shown when actually removing samples.
            # In 'keep' mode (default) nothing is removed; the preview count lives in
            # the hover note next to the relatedness plot instead of this value box.
            pihat_thresh  <- config_get(pd$config, "Filter", "relatedness", default = NULL)
            rel_action    <- tolower(config_get(pd$config, "Filter", "relatedness_action", default = "keep"))
            rel_box <- if (!is.null(pihat_thresh) && !is.na(pihat_thresh) && rel_action == "remove") {
                n_rel <- sv("samples_removed_relatedness", default = 0)
                rel_theme <- if (is.na(n_rel) || n_rel == 0) "success" else "warning"
                n_rel_str <- if (is.na(n_rel)) "—" else as.character(as.integer(n_rel))
                bslib::value_box(
                    title    = "Related Samples",
                    value    = paste0(n_rel_str, " removed"),
                    theme    = rel_theme,
                    showcase = bsicons::bs_icon("people")
                )
            } else NULL

            extra_boxes <- Filter(Negate(is.null), list(het_box, depth_box, rel_box))
            n_extra <- length(extra_boxes)
            width <- if (n_extra == 0) 1/4 else if (n_extra == 1) 1/5 else 1/6

            bslib::layout_column_wrap(
                width = width,
                fill  = FALSE,
                samp_box, snp_box, ld_box, titv_box,
                !!!extra_boxes
            )
        })

        # ── LEA conversion losses alert ────────────────────────────────────────
        # Converted from a full-width alert banner to a compact hover badge (same
        # wording, now in the tooltip) — tidy FYI status, not a page-wide banner.
        output$lea_losses_alert <- shiny::renderUI({
            pd <- project_data()
            p  <- removed_snps_path(pd$name)
            if (!nzchar(p) || !file.exists(p)) return(NULL)
            lines <- tryCatch(readLines(p, warn = FALSE), error = function(e) character(0))
            n <- length(lines)
            if (n == 0) return(NULL)
            htmltools::div(
                class = "d-flex justify-content-end mb-2",
                filter_note(
                    paste0(n, " removed"),
                    htmltools::p(
                        htmltools::tags$strong(n, if (n == 1) "SNP" else "SNPs",
                                           "removed during LEA format conversion"),
                        " (multi-allelic or invariant after VCF\u2192LFMM conversion). ",
                        "These are excluded from PCA, sNMF, and association analyses."
                    ),
                    class = "bg-secondary"
                )
            )
        })

        # Keep legacy outputs in case anything still references them (safe no-ops)
        output$samples_display <- shiny::renderUI(NULL)
        output$snps_display    <- shiny::renderUI(NULL)
        output$snps_ld         <- shiny::renderText("\u2014")
        output$titv            <- shiny::renderText("\u2014")

        # ── Image cards ────────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()

            mod_image_card_server("attrition",
                path    = shiny::reactive(qc_plot_path(pd$name, "filtering_attrition.png")),
                title   = shiny::reactive("Filtering Attrition"),
                dl_name = shiny::reactive(paste0(pd$name, "_filtering_attrition")),
                suggestion = shiny::reactive("Run mode=processing to generate QC plots"),
                note    = shiny::reactive(filter_note(
                    "stages",
                    "SNP counts surviving each filter stage: MAF → missingness → LD pruning."
                ))
            )
            mod_image_card_server("sample_miss",
                path    = shiny::reactive(qc_plot_path(pd$name, "sample_missingness_distribution.png")),
                title   = shiny::reactive("Sample Missingness"),
                dl_name = shiny::reactive(paste0(pd$name, "_sample_missingness")),
                note    = shiny::reactive({
                    t <- config_get(pd$config, "Filter", "sample_miss", default = NULL)
                    if (is.null(t)) return(NULL)
                    filter_note(
                        paste0("miss > ", t),
                        paste0("Samples missing more than ", t, " of genotypes are dropped.")
                    )
                })
            )
            mod_image_card_server("het_miss",
                path    = shiny::reactive(qc_plot_path(pd$name, "het_vs_missingness.png")),
                title   = shiny::reactive("Heterozygosity vs Missingness"),
                dl_name = shiny::reactive(paste0(pd$name, "_het_vs_missingness")),
                note    = shiny::reactive({
                    sd_thresh <- config_get(pd$config, "Filter", "het_outlier_sd", default = NULL)
                    miss_thresh <- config_get(pd$config, "Filter", "sample_miss", default = NULL)
                    if (is.null(sd_thresh) && is.null(miss_thresh)) return(NULL)
                    txt <- if (!is.null(sd_thresh)) {
                        paste0("Samples beyond ±", sd_thresh,
                              " SD heterozygosity or above the missingness cutoff are flagged as outliers.")
                    } else {
                        "Samples above the missingness cutoff are flagged as outliers."
                    }
                    filter_note(
                        if (!is.null(sd_thresh)) paste0("±", sd_thresh, " SD") else "outliers",
                        txt
                    )
                })
            )
            mod_image_card_server("relatedness",
                path    = shiny::reactive(qc_plot_path(pd$name, "relatedness_distribution.png")),
                title   = shiny::reactive("Relatedness (IBS allele-sharing)"),
                dl_name = shiny::reactive(paste0(pd$name, "_relatedness_distribution")),
                note    = shiny::reactive({
                    thresh <- config_get(pd$config, "Filter", "relatedness", default = NULL)
                    if (is.null(thresh) || is.na(thresh)) return(NULL)
                    action <- config_get(pd$config, "Filter", "relatedness_action", default = "keep")
                    pairs   <- load_relatedness_pairs(pd$name)
                    removed <- load_relatedness_removed(pd$name)
                    n_pairs_above <- if (nrow(pairs) == 0) NA_integer_
                                      else sum(pairs$IBS > thresh, na.rm = TRUE)
                    n_would_remove <- if (nrow(removed) == 0) NA_integer_ else nrow(removed)
                    relatedness_note(thresh, action, n_pairs_above, n_would_remove)
                })
            )
            mod_image_card_server("snp_miss",
                path    = shiny::reactive(qc_plot_path(pd$name, "snp_missingness_distribution.png")),
                title   = shiny::reactive("SNP Missingness"),
                dl_name = shiny::reactive(paste0(pd$name, "_snp_missingness")),
                note    = shiny::reactive({
                    t <- config_get(pd$config, "Filter", "snp_miss", default = NULL)
                    if (is.null(t)) return(NULL)
                    filter_note(
                        paste0("miss > ", t),
                        paste0("SNPs missing in more than ", t, " of samples are dropped.")
                    )
                })
            )
            mod_image_card_server("maf",
                path    = shiny::reactive(qc_plot_path(pd$name, "maf_distribution.png")),
                title   = shiny::reactive("MAF Distribution"),
                dl_name = shiny::reactive(paste0(pd$name, "_maf_distribution")),
                note    = shiny::reactive({
                    t <- config_get(pd$config, "Filter", "maf", default = NULL)
                    if (is.null(t)) return(NULL)
                    filter_note(
                        paste0("MAF < ", t),
                        paste0("Variants with minor-allele frequency below ", t,
                              " are removed before analysis.")
                    )
                })
            )
            mod_image_card_server("snp_density",
                path    = shiny::reactive(qc_plot_path(pd$name, "snp_density_by_chr.png")),
                title   = shiny::reactive("SNP Density by Chromosome"),
                dl_name = shiny::reactive(paste0(pd$name, "_snp_density")),
                note    = shiny::reactive(help_note("snp_density"))
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

        # ── Depth section (always shown — placeholder explains unavailable states) ──
        # Three states written by write_summary.R: not_provided / declared_empty / available.
        depth_status <- shiny::reactive({
            s <- summary_data()
            if (nrow(s) == 0) return("not_provided")
            row <- s[s$step == "processing" & s$metric == "depth_qc", ]
            if (nrow(row) == 0) return("not_provided")
            row$value[1]
        })

        output$depth_section <- shiny::renderUI({
            htmltools::tagList(
                mod_image_card_ui(ns("depth")),
                shiny::br()
            )
        })

        shiny::observe({
            pd <- project_data()
            mod_image_card_server("depth",
                path    = shiny::reactive(qc_plot_path(pd$name, "depth_distribution.png")),
                title   = shiny::reactive("Depth Distribution"),
                dl_name = shiny::reactive(paste0(pd$name, "_depth_distribution")),
                placeholder = shiny::reactive(switch(depth_status(),
                    not_provided   = "VCF has no FORMAT/DP field",
                    declared_empty = "FORMAT/DP is declared but genotype DP values are empty",
                    "Plot not available"
                )),
                suggestion = shiny::reactive(switch(depth_status(),
                    not_provided   = "Depth QC needs per-genotype DP in the VCF.",
                    declared_empty = "Typical for GBS — no usable depth values. Re-call with DP populated to enable depth QC.",
                    NULL
                )),
                note    = shiny::reactive({
                    min_d <- config_get(pd$config, "Filter", "min_depth", default = NULL)
                    max_d <- config_get(pd$config, "Filter", "max_depth", default = NULL)
                    if (is.null(min_d) && is.null(max_d)) return(NULL)
                    lo <- if (!is.null(min_d)) min_d else "–"
                    hi <- if (!is.null(max_d)) max_d else "–"
                    filter_note(
                        paste0(lo, "–", hi, "×"),
                        paste0("Sites outside the ", lo, "–", hi,
                              "× mean-depth window are removed.")
                    )
                })
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
