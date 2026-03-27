#' Build the application UI
#' @noRd
app_ui <- function(request) {
    pipeline_path <- get_pipeline_path()
    projects <- find_projects(pipeline_path)

    # Serve pipeline output files as static resources
    shiny::addResourcePath("pipeline", pipeline_path)

    bslib::page_navbar(
        title = htmltools::strong("ADAPTOGENE", style = "letter-spacing: 0.05em;"),
        theme    = app_theme(),
        id       = "main_tabs",
        fillable = c("Home"),

        # Enable shinyjs and load config-dirty.js
        header = htmltools::tagList(
            shinyjs::useShinyjs(),
            htmltools::tags$script(src = "www/config-dirty.js")
        ),

        # ── Navbar right: project selector + dark mode toggle ──────────────────
        nav_item(
            htmltools::div(
                class = "d-flex align-items-center gap-3",
                shiny::selectInput(
                    "project_selector",
                    label    = NULL,
                    choices  = if (length(projects) > 0) projects else "No projects found",
                    selected = if (length(projects) > 0) projects[1] else NULL,
                    width    = "220px"
                ),
                bslib::input_dark_mode(id = "dark_mode")
            )
        ),

        # ── Tabs ───────────────────────────────────────────────────────────────

        bslib::nav_panel("Home",
            class = "bslib-page-dashboard",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_home", "Processing Config",
                    runner_ui = mod_pipeline_runner_ui("runner_home",
                                                       btn_label = "Run Processing")
                ),
                fill = TRUE, fillable = TRUE,
                mod_home_ui("home")
            )
        ),

        bslib::nav_panel("Structure",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_structure", "Structure Config",
                    runner_ui = mod_pipeline_runner_ui("runner_structure",
                                                       btn_label = "Run Structure")
                ),
                fill = FALSE, fillable = FALSE,
                mod_structure_ui("structure")
            )
        ),

        bslib::nav_panel("Structure K",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_structure_k", "Structure K Config",
                    runner_ui = mod_pipeline_runner_ui("runner_structure_k",
                                                       btn_label = "Run Structure K")
                ),
                fill = FALSE, fillable = FALSE,
                mod_structure_k_ui("structure_k")
            )
        ),

        bslib::nav_panel("Association",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_association", "GEA Config",
                    runner_ui = mod_pipeline_runner_ui("runner_association",
                                                       btn_label = "Run GEA")
                ),
                fill = FALSE, fillable = FALSE,
                mod_association_ui("association")
            )
        ),

        bslib::nav_panel("Phenotype Association",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_phenotype", "GWAS Config",
                    runner_ui = mod_pipeline_runner_ui("runner_phenotype",
                                                       btn_label = "Run GWAS")
                ),
                fill = FALSE, fillable = FALSE,
                mod_phenotype_ui("phenotype")
            )
        ),

        bslib::nav_panel("Overlapping Regions",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_overlapping", "Overlap Config",
                    runner_ui = mod_pipeline_runner_ui("runner_overlapping",
                                                       btn_label = "Run Overlapping")
                ),
                fill = FALSE, fillable = FALSE,
                mod_overlapping_ui("overlapping")
            )
        ),

        bslib::nav_panel("Haplotype Analysis",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_haplotype", "Haplotype Config",
                    runner_ui = htmltools::div(
                        class = "d-flex gap-2 flex-wrap",
                        mod_pipeline_runner_ui("runner_hap_scan", btn_label = "Run Scan"),
                        mod_pipeline_runner_ui("runner_hap_viz",  btn_label = "Run Haplotype")
                    )
                ),
                fill = FALSE, fillable = FALSE,
                mod_haplotype_ui("haplotype")
            )
        ),

        bslib::nav_panel("Maladaptation",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_maladaptation", "Maladaptation Config",
                    runner_ui = mod_pipeline_runner_ui("runner_maladaptation",
                                                       btn_label = "Run Maladaptation")
                ),
                fill = FALSE, fillable = FALSE,
                mod_maladaptation_ui("maladaptation")
            )
        )
    )
}

#' Wrapper for nav_item (bslib helper)
#' @noRd
nav_item <- function(...) {
    bslib::nav_item(...)
}
