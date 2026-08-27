# smoke_suspend.R — Step 0 of the layout lab.
#
# Question: does a nested renderUI -> imageOutput chain resume when its
# container is revealed?  Variants A and C both bet their viewport budget on
# "yes"; suspendWhenHidden cannot be turned off from outside the module.
#
# Case 1 (2 levels): image card in an INACTIVE navset panel      -> Variant C
# Case 2 (3 levels): image card in a CLOSED accordion via uiOutput -> Variant A + depth_section
#
# Run:  docker run ... Rscript /pipeline/scripts/layout_lab/smoke_suspend.R

options(shiny.autoreload = FALSE)
.dev_user_lib <- "/pipeline/.R_libs_dev"
if (dir.exists(.dev_user_lib)) .libPaths(c(.dev_user_lib, .libPaths()))

for (p in c("shiny","bslib","bsicons","htmltools","data.table","DT","yaml",
            "jsonlite","cachem","config","shinyjs")) library(p, character.only = TRUE)

app_r_dir <- "/pipeline/scripts/adaptogene.app/R"
for (f in sort(list.files(app_r_dir, pattern = "\\.R$", full.names = TRUE))) source(f, local = FALSE)

shiny::addResourcePath("pipeline", "/pipeline")
shiny::addResourcePath("www", "/pipeline/scripts/adaptogene.app/inst/app/www")

PNG1 <- "/pipeline/SIMDATA_results/Processing/plots/maf_distribution.png"
PNG2 <- "/pipeline/SIMDATA_results/Processing/plots/snp_missingness_distribution.png"
cat("PNG1 exists:", file.exists(PNG1), " PNG2 exists:", file.exists(PNG2), "\n")

ui <- bslib::page_fluid(
  theme = app_theme(),
  shinyjs::useShinyjs(),
  htmltools::h4("Case 1 — image card inside an INACTIVE navset panel (2 levels)"),
  bslib::navset_card_underline(
    id = "nav",
    bslib::nav_panel("Stop 1 (active)", value = "s1", htmltools::p("nothing here")),
    bslib::nav_panel("Stop 2 (starts hidden)", value = "s2",
                     htmltools::div(id = "case1", mod_image_card_ui("t1")))
  ),
  htmltools::hr(),
  htmltools::h4("Case 2 — image card inside uiOutput inside a CLOSED accordion (3 levels)"),
  bslib::accordion(
    id = "acc", open = FALSE,
    bslib::accordion_panel("Depth-like section", value = "p1",
                           htmltools::div(id = "case2", shiny::uiOutput("section")))
  )
)

server <- function(input, output, session) {
  mod_image_card_server("t1",
    path  = shiny::reactive(PNG1),
    title = shiny::reactive("Case 1 — MAF"))

  output$section <- shiny::renderUI({
    message(">>> case2 renderUI EXECUTED at ", format(Sys.time(), "%H:%M:%S"))
    mod_image_card_ui("t2")
  })
  mod_image_card_server("t2",
    path  = shiny::reactive(PNG2),
    title = shiny::reactive("Case 2 — SNP missingness"))
}

shiny::shinyApp(ui, server, options = list(host = "0.0.0.0", port = 3838))
