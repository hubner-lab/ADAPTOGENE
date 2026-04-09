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
                    choices  = c(
                        "+ New Project" = "__new__",
                        if (length(projects) > 0) setNames(projects, projects) else character(0)
                    ),
                    selected = if (length(projects) > 0) projects[1] else "__new__",
                    width    = "220px"
                ),
                bslib::input_dark_mode(id = "dark_mode", mode = "dark")
            )
        ),

        # ── Tabs ───────────────────────────────────────────────────────────────

        bslib::nav_panel("Home",
            class = "bslib-page-dashboard",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_home", "Project Files",
                    runner_ui = shiny::actionButton(
                        "save_project_files",
                        label = htmltools::tagList(
                            bsicons::bs_icon("floppy", size = "0.9em"),
                            " Save Project Files"
                        ),
                        class = "btn btn-primary btn-sm w-100 mt-1"
                    )
                ),
                fill = TRUE, fillable = TRUE,
                mod_home_ui("home")
            )
        ),

        bslib::nav_panel("Processing",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_processing", "Processing Config",
                    runner_ui = mod_pipeline_runner_ui("runner_processing",
                                                       btn_label = "Run Processing")
                ),
                fill = FALSE, fillable = FALSE,
                mod_processing_ui("processing")
            )
        ),

        bslib::nav_panel("PreStructure",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_prestructure", "PreStructure Config",
                    runner_ui = mod_pipeline_runner_ui("runner_prestructure",
                                                       btn_label = "Run PreStructure")
                ),
                fill = FALSE, fillable = FALSE,
                mod_prestructure_ui("prestructure")
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

        bslib::nav_panel("GEA",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_gea", "GEA Config",
                    runner_ui = mod_pipeline_runner_ui("runner_gea",
                                                       btn_label = "Run GEA")
                ),
                fill = FALSE, fillable = FALSE,
                mod_gea_ui("gea")
            )
        ),

        bslib::nav_panel("GWAS",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_gwas", "GWAS Config",
                    runner_ui = mod_pipeline_runner_ui("runner_gwas",
                                                       btn_label = "Run GWAS")
                ),
                fill = FALSE, fillable = FALSE,
                mod_gwas_ui("gwas")
            )
        ),

        bslib::nav_panel("GEAxGWAS",
            bslib::layout_sidebar(
                sidebar = mod_config_sidebar_ui(
                    "config_gea_x_gwas", "GEAxGWAS Config",
                    runner_ui = mod_pipeline_runner_ui("runner_gea_x_gwas",
                                                       btn_label = "Run GEAxGWAS")
                ),
                fill = FALSE, fillable = FALSE,
                mod_gea_x_gwas_ui("gea_x_gwas")
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
