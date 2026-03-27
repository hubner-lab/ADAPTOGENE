#' Main application server
#' @noRd
app_server <- function(input, output, session) {

    # ── Reactive project data bundle (lazy, cached per project) ───────────────
    project_data <- shiny::reactive({
        shiny::req(input$project_selector)
        proj <- input$project_selector
        load_cached(paste0("project_data_", proj), function() {
            make_project_data(proj)
        })
    }) |> shiny::bindCache(input$project_selector)

    # ── Config sidebar servers (one per tab) ───────────────────────────────────
    mod_config_sidebar_server("config_home",          project_data, "home")
    mod_config_sidebar_server("config_structure",     project_data, "structure")
    mod_config_sidebar_server("config_structure_k",   project_data, "structure_k")
    mod_config_sidebar_server("config_association",   project_data, "association")
    mod_config_sidebar_server("config_phenotype",     project_data, "phenotype")
    mod_config_sidebar_server("config_overlapping",   project_data, "overlapping")
    mod_config_sidebar_server("config_haplotype",     project_data, "haplotype")
    mod_config_sidebar_server("config_maladaptation", project_data, "maladaptation")

    # ── Module servers ─────────────────────────────────────────────────────────
    mod_home_server("home",               project_data = project_data)
    mod_structure_server("structure",     project_data = project_data)
    mod_structure_k_server("structure_k", project_data = project_data)
    mod_association_server("association", project_data = project_data)
    mod_phenotype_server("phenotype",     project_data = project_data)
    mod_overlapping_server("overlapping", project_data = project_data)
    mod_maladaptation_server("maladaptation", project_data = project_data)
    mod_haplotype_server("haplotype",     project_data = project_data)
}
