# dev_layouts.R — Layout Lab entry point.
#
# Boots the REAL ADAPTOGENE app, then masks only `mod_processing_ui()` with a
# variant selected by the ?v= query string. Nothing under scripts/adaptogene.app/
# is edited; masking works because dev.R-style sourcing puts every definition in
# the global env, so a later source() shadows an earlier one.
#
#   ?v=a  Command Center   ?v=b  Analyst Workbench   ?v=c  Narrative Flow
#   (no ?v=, or ?v=orig)   the untouched production layout, as the control
#
# Run:
#   docker run --user $(id -u):$(id -g) --rm -e USER=pipeline -p 3839:3838 \
#     -v $PWD:/pipeline adaptogene:latest \
#     Rscript /pipeline/scripts/layout_lab/dev_layouts.R

options(shiny.autoreload = TRUE)

.dev_user_lib <- "/pipeline/.R_libs_dev"
dir.create(.dev_user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(.dev_user_lib, .libPaths()))

for (.p in c("shiny", "bslib", "bsicons", "golem", "htmltools", "data.table",
             "plotly", "DT", "yaml", "jsonlite", "base64enc", "cachem", "config",
             "processx", "shinyjs", "shinyFiles")) {
    library(.p, character.only = TRUE)
}

APP_DIR <- "/pipeline/scripts/adaptogene.app"
LAB_DIR <- "/pipeline/scripts/layout_lab"

# 1. the real app, verbatim
for (f in sort(list.files(file.path(APP_DIR, "R"), pattern = "\\.R$", full.names = TRUE))) {
    source(f, local = FALSE)
}

# 2. lab overrides — masks app_theme() and mod_processing_ui()
for (f in sort(list.files(file.path(LAB_DIR, "R"), pattern = "\\.R$", full.names = TRUE))) {
    source(f, local = FALSE)
}

shiny::addResourcePath("www", file.path(APP_DIR, "inst/app/www"))
shiny::addResourcePath("lab", file.path(LAB_DIR, "www"))

message("=== ADAPTOGENE Layout Lab ===")
message("  ?v=a  Command Center | ?v=b  Analyst Workbench | ?v=c  Narrative Flow")
message("  no ?v= -> production layout (control)")

shiny::shinyApp(
    ui      = lab_app_ui,     # sets the variant from ?v=, then calls the REAL app_ui
    server  = app_server,     # untouched
    options = list(host = "0.0.0.0", port = 3838)
)
