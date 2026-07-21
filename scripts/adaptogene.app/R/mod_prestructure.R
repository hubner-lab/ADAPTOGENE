#' PreStructure tab UI
#'
#' PCA, Tracy-Widom, cross-entropy, and K-dependent plots
#' (structure barplot, PCA-structure, pop differentiation).
#'
#' @param id module namespace id
#' @noRd
mod_prestructure_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # K-independent plots: PCA + Tracy-Widom + Cross-Entropy in one row
        bslib::layout_column_wrap(
            width = 1 / 3,
            mod_image_card_ui(ns("pca")),
            mod_image_card_ui(ns("tracy_widom")),
            mod_image_card_ui(ns("cross_entropy"))
        ),

        # K selector inline above K-dependent plots
        htmltools::div(
            class = "control-bar",
            shiny::uiOutput(ns("k_selector"))
        ),

        # K-dependent plots: Structure + PCA-structure + Pop Diff in one row
        bslib::layout_column_wrap(
            width = 1 / 3,
            mod_image_card_ui(ns("structure_bar")),
            mod_image_card_ui(ns("pca_structure")),
            mod_image_card_ui(ns("pop_diff"))
        )
    )
}

#' PreStructure tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_prestructure_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        k_values <- shiny::reactive({
            find_k_values(project_data()$name)
        })

        output$k_selector <- shiny::renderUI({
            ks <- k_values()
            if (length(ks) == 0) return(shiny::p("No K values found.", class = "text-muted small"))
            pd <- project_data()
            selected <- if (!is.na(pd$k_best) && pd$k_best %in% ks) pd$k_best else ks[1]
            shiny::selectInput(ns("k"), "K value", choices = ks, selected = selected)
        })

        selected_k <- shiny::reactive({
            k <- input$k
            if (is.null(k)) {
                ks <- k_values()
                if (length(ks) == 0) return(NULL)
                pd <- project_data()
                if (!is.na(pd$k_best) && pd$k_best %in% ks) pd$k_best else ks[1]
            } else {
                as.integer(k)
            }
        })

        # ── K-independent images ───────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            kr <- find_k_range(pd$name)
            k_best <- pd$k_best

            n_samples <- {
                dt <- if (!is.null(k_best) && !is.na(k_best)) load_clusters(pd$name, k_best)
                      else data.table::data.table()
                if (nrow(dt) == 0) NULL else nrow(dt)
            }

            mod_image_card_server("pca",
                path    = shiny::reactive(pca_path(pd$name)),
                title   = shiny::reactive("PCA"),
                dl_name = shiny::reactive("pca"),
                note    = shiny::reactive(help_note("pca",
                    results = if (!is.null(n_samples)) paste0("N = ", n_samples, " samples.") else NULL
                ))
            )
            mod_image_card_server("tracy_widom",
                path    = shiny::reactive(tracy_widom_path(pd$name)),
                title   = shiny::reactive("Tracy-Widom Test"),
                dl_name = shiny::reactive("tracy_widom"),
                note    = shiny::reactive(help_note("tracy_widom"))
            )
            mod_image_card_server("cross_entropy",
                path    = shiny::reactive(kr$path),
                title   = shiny::reactive("Cross-Entropy"),
                dl_name = shiny::reactive("cross_entropy"),
                note    = shiny::reactive({
                    if (is.null(k_best) || is.na(k_best)) return(NULL)
                    range_txt <- if (!is.na(kr$k_start) && !is.na(kr$k_end))
                        paste0(" Tested K ", kr$k_start, "–", kr$k_end, ".") else ""
                    # config_k_best() is NA when sNMF.k_best is unset — resolve_k_best()
                    # then falls back to the highest K found on disk, which is NOT the
                    # lowest-cross-entropy K. Word the note accordingly so it never
                    # overclaims a "best" that wasn't actually computed from this plot.
                    results_txt <- if (!is.na(config_k_best(pd$config))) {
                        paste0("Selected K = ", k_best, " (sNMF.k_best).", range_txt)
                    } else {
                        paste0("sNMF.k_best is not set — showing the highest K found on disk (",
                              k_best, "), not the lowest cross-entropy value.", range_txt,
                              " Read the curve above and set sNMF.k_best to the K where it flattens.")
                    }
                    help_note("cross_entropy", results = results_txt, label = paste0("K=", k_best))
                })
            )
        })

        # ── K-dependent images ─────────────────────────────────────────────────
        shiny::observe({
            k  <- shiny::req(selected_k())
            pd <- project_data()

            pop_summary <- cluster_pop_summary(load_clusters(pd$name, k), k)
            pop_results <- if (!is.null(pop_summary)) {
                rng <- if (pop_summary$min_n == pop_summary$max_n) as.character(pop_summary$min_n)
                       else paste0(pop_summary$min_n, "–", pop_summary$max_n)
                paste0(pop_summary$n_pops, " populations at K=", k, " — ",
                      rng, " samples/pop (N=", pop_summary$n_samples, ").")
            } else NULL

            mod_image_card_server("structure_bar",
                path    = shiny::reactive(structure_k_path(pd$name, k)),
                title   = shiny::reactive("Structure"),
                dl_name = shiny::reactive(paste0("structure_K", k)),
                note    = shiny::reactive(help_note("structure_bar", results = pop_results,
                                                    label = paste0("K=", k)))
            )
            mod_image_card_server("pca_structure",
                path    = shiny::reactive(pca_structure_path(pd$name, k)),
                title   = shiny::reactive("PCA (Clustered)"),
                dl_name = shiny::reactive(paste0("pca_structure_K", k)),
                note    = shiny::reactive(help_note("pca_structure", results = pop_results,
                                                    label = paste0("K=", k)))
            )
            mod_image_card_server("pop_diff",
                path    = shiny::reactive(pop_diff_path(pd$name, k)),
                title   = shiny::reactive("Pop Differentiation"),
                dl_name = shiny::reactive(paste0("pop_diff_K", k)),
                note    = shiny::reactive(help_note("pop_diff", results = pop_results,
                                                    label = paste0("K=", k)))
            )
        })
    })
}
