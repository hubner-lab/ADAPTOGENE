# lab_85_gea_a.R — GEA in the Variant A pattern, adapted for an interactive module.
#
# The pattern bends here, deliberately:
#   * The combined Manhattan is a genome-wide x-axis. It gets the FULL width as
#     the hero rather than a 5/12 anchor column -- compressing it by 40% would
#     cost real resolution on the one plot that needs it.
#   * There is NO KPI strip. Every number a GEA value box could show
#     (significant SNPs, regions, genes) is read from pipeline_summary.tsv,
#     which records what the PIPELINE computed at ITS configured threshold.
#     The tab's whole point is that the user re-derives those numbers live from
#     the threshold bar and the trait x method matrix, so the boxes contradicted
#     the screen the moment anything was touched (measured 2026-08-27:
#     `gea/selected_snps_total = 8` while the on-screen matrix summed to more)
#     and never updated. Static numbers next to live controls are worse than no
#     numbers. Live counts are an app-level change, not a layout one -- see the
#     note on lab_gea_extras() below.
#   * The threshold bar, WZA note, trait x method matrix and Save button live in
#     a COLLAPSED accordion above the plot. They are set once and then read
#     rarely, but in the old layout they occupied ~230px of permanent vertical
#     budget between the card header and the plot. Collapsing them (plus
#     dropping the KPI strip) is what buys the combined Manhattan AND the
#     per-method detail row on one screen.
#   * Region Explorer is left exactly as-is, full width at the bottom. Its server
#     is 1342 lines; there is nothing to gain from re-housing it.

mod_gea_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "gea",

        shiny::uiOutput(ns("config_badges")),

        # ── Controls: collapsed by default, above the plot they drive ─────────
        # Safe to collapse: mod_gea_server() has exactly one req() in 729 lines
        # (on project_data(), inside the save observer), and every control input
        # it reads is guarded by `%||%` fallback -- input$threshold_type %||%
        # "bonf", input$combine_strategy %||% default_strategy(), and so on. A
        # suspended control bar therefore yields pipeline defaults rather than a
        # blocked plot. lab_gea_extras() additionally un-suspends the outputs so
        # the values saved in region_params.json still reach them.
        bslib::accordion(
            id    = ns("params_accordion"),
            class = "lab-gea-params",
            open  = FALSE,
            bslib::accordion_panel(
                value = "params",
                title = htmltools::tagList(
                    bsicons::bs_icon("sliders"),
                    htmltools::span(" Significance, methods & strategy")
                ),
                htmltools::div(
                    class = "lab-gea-controls",
                    shiny::uiOutput(ns("threshold_bar")),
                    shiny::uiOutput(ns("wza_collapse_note")),
                    shiny::uiOutput(ns("filter_bar")),
                    htmltools::div(
                        class = "d-flex justify-content-end mt-1",
                        shiny::actionButton(
                            ns("save_snp_set"),
                            label = htmltools::tagList(
                                bsicons::bs_icon("save"),
                                " Save SNP set for maladaptation"
                            ),
                            class = "btn-sm btn-success"
                        )
                    )
                )
            )
        ),

        # ── HERO: combined Manhattan, full width, no injected filter_ui ───────
        htmltools::div(
            class = "lab-gea-hero",
            mod_manhattan_overlay_ui(ns("combined_manhattan"), height = "100%")
        ),

        # ── Per-method detail ─────────────────────────────────────────────────
        lab_section_header("Per-method detail", icon = "layers"),
        htmltools::div(
            class = "lab-gea-methodbar d-flex align-items-center gap-3 flex-wrap",
            shiny::uiOutput(ns("method_tabs_ui")),
            shiny::uiOutput(ns("per_method_trait_ui"))
        ),
        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",
            # QQ is square-ish -> the narrow column; the per-method Manhattan is
            # wide -> the broad one. Inverted vs the display modules on purpose.
            htmltools::div(class = "lab-gea-qq", mod_image_card_ui(ns("qq_plot"))),
            htmltools::div(class = "lab-gea-methodman",
                           mod_manhattan_overlay_ui(ns("method_manhattan"),
                                                    height = "100%"))
        ),

        # ── RDA diagnostics: collapsed (suspended until opened) ───────────────
        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE,
            bslib::accordion_panel(
                value = "rda",
                title = htmltools::tagList(bsicons::bs_icon("compass"),
                                           " RDA diagnostics & candidates"),
                mod_rda_details_ui(ns("rda_details"))
            )
        ),

        # ── Region Explorer: unchanged, full width, bottom ────────────────────
        mod_region_explorer_ui(ns("region_explorer"))
    )
}

#' Additive server for the GEA lab layout.
#'
#' There is no KPI strip to feed any more (see the header note). What is left is
#' one correctness fix that the collapsed control accordion makes necessary.
#'
#' Shiny suspends an output whose DOM node is hidden, and a closed accordion
#' panel is `display: none`. `threshold_bar` and `filter_bar` are `renderUI`
#' outputs, so collapsing them by default would mean they never render on load
#' -- and the `observe()` in mod_gea_server() that replays the user's saved
#' threshold from region_params.json via updateSelectInput/updateNumericInput
#' would push into inputs that do not exist yet. The messages are dropped, and
#' when the panel is finally opened the bar rebuilds from
#' `input$threshold_type %||% "bonf"` -- silently discarding the saved value.
#' That is exactly the "never silently discard a user's chosen parameter value"
#' failure the app's UI rules call out.
#'
#' Un-suspending the two outputs keeps them rendering while hidden, so the
#' inputs exist from load and the replay lands. This runs AFTER
#' .lab_orig_mod_gea_server() (see lab_90_dispatch.R), so the outputs are
#' already registered on the session by the time we set the option.
lab_gea_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        for (out_id in c("threshold_bar", "filter_bar", "wza_collapse_note")) {
            ok <- tryCatch({
                shiny::outputOptions(output, out_id, suspendWhenHidden = FALSE)
                TRUE
            }, error = function(e) {
                message("[lab] could not un-suspend '", out_id, "': ", conditionMessage(e))
                FALSE
            })
            if (!ok) next
        }
    })
}
