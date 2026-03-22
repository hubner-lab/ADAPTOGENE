# dev.R — Development runner for ADAPTOGENE Shiny app
# Sources R files from mounted volume instead of installed package.
# Avoids Docker rebuilds: edit files and shiny.autoreload restarts the app.
#
# Usage:
#   docker run --user $(id -u):$(id -g) --rm -p 3838:3838 -v $PWD:/pipeline adaptogene:latest \
#     Rscript /pipeline/scripts/adaptogene.app/dev.R

options(shiny.autoreload = TRUE)

# Load all dependencies (from DESCRIPTION Imports, minus base R packages)
library(shiny)
library(bslib)
library(bsicons)
library(golem)
library(htmltools)
library(data.table)
library(plotly)
library(DT)
library(yaml)
library(jsonlite)
library(base64enc)
library(cachem)
library(config)

# Source all app R files from mounted volume
app_r_dir <- "/pipeline/scripts/adaptogene.app/R"
for (f in sort(list.files(app_r_dir, pattern = "\\.R$", full.names = TRUE))) {
    source(f, local = FALSE)
}

# Serve static assets (CSS/JS/SCSS from inst/app/www/)
shiny::addResourcePath("www", "/pipeline/scripts/adaptogene.app/inst/app/www")

# Run app
shiny::shinyApp(
    ui     = app_ui,
    server = app_server,
    options = list(host = "0.0.0.0", port = 3838)
)
