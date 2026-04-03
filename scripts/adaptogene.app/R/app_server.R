#' Main application server
#' @noRd
app_server <- function(input, output, session) {

    # ── Config state: working copy + saved (YAML) ──────────────────────────────
    # config_state$saved   = last-written YAML config (matches pipeline results)
    # config_state$working = user's edits (lost on reload unless Run is clicked)
    # config_state$project = current project name
    config_state <- shiny::reactiveValues(
        saved   = list(),
        working = list(),
        project = NULL
    )

    # Populate config_state on project switch
    shiny::observeEvent(input$project_selector, {
        shiny::req(input$project_selector)
        proj   <- input$project_selector
        cfg    <- read_project_config(proj)
        config_state$project <- proj
        config_state$saved   <- cfg
        config_state$working <- cfg
    }, ignoreInit = FALSE)

    # ── Global run lock (prevents concurrent pipeline runs) ───────────────────
    pipeline_running <- shiny::reactiveVal(FALSE)

    # ── Trigger to invalidate project_data after a successful run ─────────────
    project_data_trigger <- shiny::reactiveVal(0L)

    # ── Reactive project data bundle ───────────────────────────────────────────
    # Depends on project selector AND project_data_trigger (incremented post-run)
    project_data <- shiny::reactive({
        shiny::req(input$project_selector)
        project_data_trigger()  # declare dependency
        proj <- input$project_selector
        make_project_data(proj)
    }) |> shiny::bindCache(input$project_selector, project_data_trigger())

    # ── Config sidebar servers (one per tab) ───────────────────────────────────
    mod_config_sidebar_server("config_home",          config_state, "home")
    mod_config_sidebar_server("config_processing",    config_state, "processing")
    mod_config_sidebar_server("config_structure",     config_state, "structure")
    mod_config_sidebar_server("config_structure_k",   config_state, "structure_k")
    mod_config_sidebar_server("config_association",   config_state, "association")
    mod_config_sidebar_server("config_phenotype",     config_state, "phenotype")
    mod_config_sidebar_server("config_overlapping",   config_state, "overlapping")
    mod_config_sidebar_server("config_maladaptation", config_state, "maladaptation")

    # ── Pipeline runner servers (one per tab / mode) ───────────────────────────
    mod_pipeline_runner_server("runner_home",         config_state, pipeline_running,
                                project_data_trigger, tab_name = "home")
    mod_pipeline_runner_server("runner_processing",   config_state, pipeline_running,
                                project_data_trigger, tab_name = "processing")
    mod_pipeline_runner_server("runner_structure",    config_state, pipeline_running,
                                project_data_trigger, tab_name = "structure")
    mod_pipeline_runner_server("runner_structure_k",  config_state, pipeline_running,
                                project_data_trigger, tab_name = "structure_k")
    mod_pipeline_runner_server("runner_association",  config_state, pipeline_running,
                                project_data_trigger, tab_name = "association")
    mod_pipeline_runner_server("runner_phenotype",    config_state, pipeline_running,
                                project_data_trigger, tab_name = "phenotype")
    mod_pipeline_runner_server("runner_overlapping",  config_state, pipeline_running,
                                project_data_trigger, tab_name = "overlapping")
    mod_pipeline_runner_server("runner_maladaptation",config_state, pipeline_running,
                                project_data_trigger, tab_name = "maladaptation")

    # ── Module servers ─────────────────────────────────────────────────────────
    mod_home_server("home",               project_data = project_data)
    mod_processing_server("processing",   project_data = project_data)
    mod_structure_server("structure",     project_data = project_data)
    mod_structure_k_server("structure_k", project_data = project_data)
    mod_association_server("association", project_data = project_data)
    mod_phenotype_server("phenotype",     project_data = project_data)
    mod_overlapping_server("overlapping", project_data = project_data)
    mod_maladaptation_server("maladaptation", project_data = project_data)
}
