#' Build the application UI
#' @noRd
app_ui <- function(request) {
    pipeline_path <- get_pipeline_path()
    projects <- find_projects(pipeline_path)

    # Serve pipeline output files as static resources
    shiny::addResourcePath("pipeline", pipeline_path)

    bslib::page_navbar(
        title = htmltools::span(
            class = "d-flex align-items-center gap-2",
            htmltools::strong("ADAPTOGENE", style = "letter-spacing: 0.05em;"),
            htmltools::span("Results Viewer",
                            style = "font-weight: 400; opacity: 0.75; font-size: 0.9em;")
        ),
        theme    = app_theme(),
        id       = "main_tabs",
        fillable = c("Home"),   # only Home tab fills viewport; others scroll

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
            mod_home_ui("home")
        ),
        bslib::nav_panel("Structure",
            mod_structure_ui("structure")
        ),
        bslib::nav_panel("Structure K",
            mod_structure_k_ui("structure_k")
        ),
        bslib::nav_panel("Association",
            mod_association_ui("association")
        ),
        bslib::nav_panel("Phenotype Association",
            mod_phenotype_ui("phenotype")
        ),
        bslib::nav_panel("Overlapping Regions",
            mod_overlapping_ui("overlapping")
        ),
        bslib::nav_panel("Maladaptation",
            mod_maladaptation_ui("maladaptation")
        ),
        bslib::nav_panel("Haplotype Analysis",
            mod_haplotype_ui("haplotype")
        )
    )
}

#' Wrapper for nav_item (bslib helper)
#' @noRd
nav_item <- function(...) {
    bslib::nav_item(...)
}
