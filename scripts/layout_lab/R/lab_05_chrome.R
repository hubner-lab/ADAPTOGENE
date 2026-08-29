# lab_05_chrome.R — single bottom bar instead of two top panels.
#
# The app stacks nav.navbar (brand + benchmark + project + dark mode) above
# uiOutput("module_bar"), costing ~104px of every module's vertical budget before
# any content renders. This folds the upper panel DOWN INTO the lower one, so
# there is ONE bar, still at the top:
#
#     CLINE-GO | (project selector) | module tabs | ... | (dark toggle)
#
# It is pinned: the tab content scrolls beneath it, so navigation never scrolls
# away, and the merge still buys back ~56px on every module.
#
# The Benchmark switch is dropped entirely. `is_simulated` is declared per project
# by `Simulation.enabled` in that project's own config (fct_discovery.R:7), so it
# is a property of the project, not a session mode. Removing the input is safe:
# app_server reads it as `isTRUE(input$benchmark_mode)` -> FALSE when absent, and
# app_ui's own find_projects() call already defaults to non-simulated.
#
# Structural note: the real nav.navbar is still rendered but VISUALLY hidden with
# the clip technique, never display:none -- module-bar.js drives tab switching by
# firing synthetic clicks on those anchors, and they must stay in the layout tree.

.lab_orig_app_ui <- app_ui

# value, nav label, config/runner label, module UI function name.
# The config label is NOT always the nav label: "Environmental" is configured as
# "Climate", "Phenotypic" as "Traits".
.LAB_PANELS <- list(
    list("processing",    "Processing",    "Processing",    "mod_processing_ui"),
    list("prestructure",  "PreStructure",  "PreStructure",  "mod_prestructure_ui"),
    list("structure",     "Structure",     "Structure",     "mod_structure_ui"),
    list("climate",       "Environmental", "Climate",       "mod_climate_ui"),
    list("traits",        "Phenotypic",    "Traits",        "mod_traits_ui"),
    list("pregea",        "PreGEA",        "PreGEA",        "mod_pregea_ui"),
    list("gea",           "GEA",           "GEA",           "mod_gea_ui"),
    list("gwas",          "GWAS",          "GWAS",          "mod_gwas_ui"),
    list("gea_x_gwas",    "GEAxGWAS",      "GEAxGWAS",      "mod_gea_x_gwas_ui"),
    list("maladaptation", "Maladaptation", "Maladaptation", "mod_maladaptation_ui")
)

.lab_nav_panel <- function(spec) {
    value <- spec[[1]]; nav_label <- spec[[2]]; cfg <- spec[[3]]; ui_fn <- spec[[4]]
    bslib::nav_panel(
        nav_label, value = value,
        bslib::layout_sidebar(
            sidebar = mod_config_sidebar_ui(
                paste0("config_", value), paste0(cfg, " Config"),
                runner_ui = mod_pipeline_runner_ui(paste0("runner_", value),
                                                   btn_label = paste("Run", cfg))
            ),
            fill = FALSE, fillable = FALSE,
            # base::get, NOT get: dev_layouts.R does library(config), whose
            # config::get() masks base::get() and tries to read a config.yml.
            # Resolved at call time so the LAB's masked mod_*_ui is used.
            base::get(ui_fn)(value)
        )
    )
}

app_ui <- function(request) {
    if (identical(lab_chrome(), "top")) return(.lab_orig_app_ui(request))

    pipeline_path <- get_pipeline_path()
    projects <- find_projects(pipeline_path)
    shiny::addResourcePath("pipeline", pipeline_path)

    bslib::page_navbar(
        title    = NULL,
        theme    = app_theme(),
        id       = "main_tabs",
        fillable = c("home"),

        header = htmltools::tagList(
            shinyjs::useShinyjs(),
            shiny::useBusyIndicators(spinners = TRUE),
            htmltools::tags$script(src = "www/config-dirty.js"),
            htmltools::tags$script(src = "www/module-bar.js"),

            # ONE bar, first inside container-fluid, pinned at the top.
            htmltools::div(
                class = "lab-chrome-bar",
                htmltools::div(
                    class = "lab-bar-left",
                    brand_logo(),
                    shiny::selectInput(
                        "project_selector", label = NULL,
                        choices = c(
                            "+ New Project" = "__new__",
                            if (length(projects) > 0) setNames(projects, projects)
                            else character(0)
                        ),
                        selected = if (length(projects) > 0) projects[1] else "__new__",
                        width = "200px"
                    )
                ),
                htmltools::div(class = "lab-bar-modules",
                               shiny::uiOutput("module_bar")),
                # Far right of the bar. The toggle is a session-wide display
                # preference, not part of the identity block (brand + project)
                # nor of the navigation, so it sits apart from both rather than
                # wedged between the wordmark and the project selector.
                htmltools::div(class = "lab-bar-right",
                               bslib::input_dark_mode(id = "dark_mode", mode = "dark"))
            )
        ),

        bslib::nav_panel("Home", value = "home",
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
        !!!lapply(.LAB_PANELS, .lab_nav_panel)
    )
}
