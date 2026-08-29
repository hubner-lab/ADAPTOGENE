# lab_50_prestructure_a.R — PreStructure in the approved Variant A pattern.
#
# Anchor (left, 5/12): PCA -- the "what does this dataset look like" figure, and a
# perfect 1:1 for a narrow column -- with the Q-matrix docked beneath it, the way
# Filtering Summary sat under the attrition chart in Processing.
# Right (7/12): ANCESTRY, and only ancestry -- PCA (Clustered) + the structure
# barplot, given the whole column. The K selector rides the section divider and
# the K-choice diagnostics sit in a collapsed panel underneath, so the two
# figures a reader actually comes here for are the two that get the space.
#
# DIFFERENT FROM PROCESSING: mod_prestructure_server() has no summary_boxes and no
# clusters table output -- they do not exist to reuse. So this variant ships a small
# ADDITIVE server (lab_prestructure_extras) that runs alongside the untouched one and
# only ADDS outputs. Nothing in scripts/adaptogene.app/ is modified.

mod_prestructure_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "prestructure",

        # PreStructure has no badge output, hence alert_id = NULL
        lab_kpi_row(ns, alert_id = NULL),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── LEFT: the anchor ──────────────────────────────────────────────
            htmltools::div(
                class = "lab-hero-col",
                lab_hero(mod_image_card_ui(ns("pca"))),
                bslib::card(
                    class = "lab-dock-card lab-dock-qmatrix",
                    bslib::card_header(
                        bsicons::bs_icon("table"),
                        " Ancestry coefficients (Q-matrix)"
                    ),
                    bslib::card_body(DT::DTOutput(ns("clusters_table")))
                )
            ),

            # ── RIGHT: plenty of plots ────────────────────────────────────────
            htmltools::div(
                class = "lab-multiples-col",

                # The K selector drives structure_bar / pca_structure / pop_diff,
                # all of which live in this column. PCA on the left is
                # K-independent, so the control belongs here, not page-wide --
                # and it rides the Ancestry divider rather than claiming a
                # `.control-bar` row, which cost 64px to hold a 128px select.
                lab_section_header("Ancestry", icon = "diagram-3",
                                   aside = shiny::uiOutput(ns("k_selector"))),
                # PCA (Clustered) then the barplot -- the same K-dependent story
                # told two ways -- STACKED, not side by side.
                #
                # They cannot share a row at any useful size. pca_structure is
                # 1:1 and structure_bar is 4:1, so filling one row of height H
                # needs H + 4H of width: in this 904px column that pins H at
                # ~180px. Side by side at 1+3 columns they measured 202px and
                # 166px inside 466px tiles -- both width-bound, each stranded in
                # ~270px of vertical whitespace. Stacked, each tile is sized to
                # its own aspect and the figures roughly double.
                htmltools::div(
                    class = "lab-prestr-ancestry",
                    htmltools::div(class = "lab-prestr-pca",
                                   mod_image_card_ui(ns("pca_structure"))),
                    htmltools::div(class = "lab-prestr-bar",
                                   mod_image_card_ui(ns("structure_bar")))
                ),

                # Diagnostics are read once, when choosing K -- not on every
                # visit. Collapsed, they cost one 35px header instead of a full
                # tile row, and that reclaimed height goes to the Ancestry pair
                # above, which is what this column is FOR. (A closed panel is
                # display:none, so the three PNGs are not transmitted until it
                # is opened -- see lab_depth_below_fold()'s note.)
                bslib::accordion(
                    class = "lab-below-fold lab-prestr-diag",
                    open = FALSE,
                    bslib::accordion_panel(
                        value = "diag",
                        title = htmltools::tagList(bsicons::bs_icon("graph-up"),
                                                   " Diagnostics (Tracy-Widom, cross-entropy, population differentiation)"),
                        # tracy_widom and cross_entropy are 1:1 -> one column each.
                        # pop_diff is 1.33 -> two columns, filling the row exactly.
                        lab_thumb_grid(
                            cols = 4,
                            lab_thumb(mod_image_card_ui(ns("tracy_widom"))),
                            lab_thumb(mod_image_card_ui(ns("cross_entropy"))),
                            lab_thumb(mod_image_card_ui(ns("pop_diff")), span = 2)
                        )
                    )
                )
            )
        )
    )
}

#' Additive server: supplies ONLY the two outputs the real server does not have.
#'
#' Runs as a second moduleServer on the SAME id, so it shares the namespace and can
#' read `input$k` from the real server's own K selector.
lab_prestructure_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {

        selected_k <- shiny::reactive({
            k <- input$k
            if (!is.null(k)) return(as.integer(k))
            pd <- project_data()
            ks <- find_k_values(pd$name)
            if (length(ks) == 0) return(NULL)
            if (!is.na(pd$k_best) && pd$k_best %in% ks) pd$k_best else ks[1]
        })

        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))

            sv <- function(step, metric, default = "—") {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == step & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            }

            # NOT find_k_range(): it parses whichever cross_entropy_K*.png sorts
            # first, and a stale K2-6 file sitting beside the current K2-7.0 one
            # makes it report the WRONG range (logged in
            # docs/pipeline_improvement_requests.md). The K* directories on disk
            # are the K values actually computed, so derive from those and fall
            # back to the summary table, which is authoritative.
            ks <- find_k_values(pd$name)
            k_range_txt <- if (length(ks) > 0) {
                paste0(min(ks), " – ", max(ks))
            } else sv("prestructure", "K_range")

            k_best_txt <- if (!is.null(pd$k_best) && !is.na(pd$k_best))
                as.character(pd$k_best) else "—"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("K best",  k_best_txt,  "primary", "diagram-3"),
                lab_value_box("K tested", k_range_txt, "info",   "list-ol"),
                lab_value_box("Samples",
                              sv("processing", "samples_after_filtering"),
                              "success", "people-fill"),
                lab_value_box("SNPs (LD-pruned)",
                              sv("processing", "snps_after_ld_pruning"),
                              "info", "scissors")
            )
        })

        output$clusters_table <- DT::renderDataTable({
            pd <- project_data()
            k  <- selected_k()
            shiny::req(k)
            df <- as.data.frame(load_clusters(pd$name, k))
            shiny::validate(shiny::need(nrow(df) > 0, "No Q-matrix for this K"))
            # C1..CK are ancestry proportions -- 3 dp is plenty and keeps the
            # table narrow enough for a 5/12 column.
            num <- vapply(df, is.numeric, logical(1))
            df[num] <- lapply(df[num], function(x) round(x, 3))
            safe_datatable(
                df,
                options  = list(dom = "tp", scrollX = TRUE, pageLength = 6),
                rownames = FALSE
            )
        })
    })
}
