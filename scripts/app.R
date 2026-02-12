#TODO change pipeline_path to /pipeline

library(shiny)
library(shinydashboard)
library(DT)
#library(tidyverse)
library(data.table)
library(stringr)
library(qs)  # for reading qs files
library(ggplot2)  # for rendering ggplot objects
library(ggplotify) # for as.ggplot object load
library(gridExtra) # grid.arrange
library(plotly)
library(htmlwidgets)
library(scattermore)
# Function to scan directory and find projects

# Function to find available K values from filenames
find_k_values <- function(project, pipeline_path = '/pipeline/') {
  # Get all files in the intermediate directory
  files <- list.files(paste0(pipeline_path, project, '_results/intermediate'))

  # Extract K values from different plot types
  k_values <- unique(as.numeric(stringr::str_extract(
    grep("^(pca_structure_K|pop_diff_K|structure_K)[0-9]+\\.qs$",
         files, value = TRUE),
    "(?<=K)[0-9]+"
  )))

  # Sort K values
  k_values <- sort(k_values[!is.na(k_values)])

  return(k_values)
}

find_bio_values <- function(project, pipeline_path = '/pipeline/') {
  # Get all files in the intermediate directory
  files <- list.files(paste0(pipeline_path, project, '_results/intermediate'))

  # Extract K values from different plot types
  bio_values <- unique(as.numeric(stringr::str_extract(
    grep("PieMap_bio_[0-9]+_TajimaD.qs",
         files, value = TRUE),
    "PieMap_bio_([0-9]+)_TajimaD.qs", group = 1
  )))

  # Sort K values
  bio_values <- sort(bio_values)

  return(bio_values)
}

find_metric_values <- function(project, pipeline_path = '/pipeline/') {
  # Get all files in the intermediate directory
  files <- list.files(paste0(pipeline_path, project, '_results/intermediate'))

  # Extract K values from different plot types
  metric_values <- unique(stringr::str_extract(
    grep("^PieMap_bio_\\d+_.*\\.qs$",
         files, value = TRUE),
    "^PieMap_bio_\\d+_(.*)\\.qs$", group = 1
  ))
							message(metric_values)
  # Sort K values
  metric_values <- metric_values

  return(metric_values)
}

find_projects <- function(pipeline_path = '/pipeline') {
  # List all directories in the pipeline path
  all_dirs <- list.dirs(pipeline_path, 
			full.names = FALSE, recursive = FALSE) %>%
		grep('_results$', value = T, .)
  
  # Find unique project names by removing _results and _logs suffixes
  projects <- all_dirs %>% gsub('_results', '', .)
  
  # Validate each project has both results and logs directories
  # valid_projects <- projects[sapply(projects, function(proj) {
    # all(file.exists(file.path(pipeline_path, paste0(proj, c("_results", "_logs"))))
  # })]
  
  return(projects)
}

###############################################################################
# MALADAPTATION HELPER FUNCTIONS
###############################################################################

# Discover available Gradient Forest suffixes from intermediate .qs files
find_gf_suffixes <- function(project, pipeline_path = '/pipeline/') {
  inter <- paste0(pipeline_path, project, '_results/intermediate')
  if (!dir.exists(inter)) return(character(0))
  files <- list.files(inter, pattern = '^GeneticOffsetPieMap_.*\\.qs$')
  # Strip prefix and .qs, remove variant suffixes (PiDiversity, TajimaD, zoom_*)
  base <- sub('\\.qs$', '', sub('^GeneticOffsetPieMap_', '', files))
  base <- base[!grepl('_(PiDiversity|TajimaD|zoom_)', base)]
  sort(unique(base))
}

# Discover available zoom region tags for a GF suffix
find_gf_zooms <- function(project, suffix, pipeline_path = '/pipeline/') {
  inter <- paste0(pipeline_path, project, '_results/intermediate')
  files <- list.files(inter,
    pattern = paste0('^GeneticOffsetPieMap_', suffix, '_zoom_.*\\.qs$'))
  sub('\\.qs$', '', sub(paste0('^GeneticOffsetPieMap_', suffix, '_zoom_'), '', files))
}

###############################################################################
# ASSOCIATION ANALYSIS HELPER FUNCTIONS
###############################################################################

# Find K value(s) from association output files
find_assoc_k <- function(project, pipeline_path = '/pipeline/') {
  base <- paste0(pipeline_path, project, '_results/tables/association')
  if (!dir.exists(base)) return(NULL)

  method_dirs <- list.dirs(base, recursive = FALSE, full.names = FALSE)
  method_dirs <- method_dirs[!method_dirs %in% c('enrichment')]
  if (length(method_dirs) == 0) return(NULL)

  files <- list.files(file.path(base, method_dirs[1]),
                      pattern = '_pvalues_K\\d+\\.tsv$')
  k <- unique(as.integer(stringr::str_extract(files, '(?<=_K)\\d+')))
  sort(k[!is.na(k)])
}

# Load all association data for a project
load_assoc_data <- function(project, k, pipeline_path = '/pipeline/') {
  base <- paste0(pipeline_path, project, '_results/tables/association')

  safe_fread <- function(path) {
    if (file.exists(path)) data.table::fread(path) else NULL
  }

  method_dirs <- list.dirs(base, recursive = FALSE, full.names = FALSE)
  method_dirs <- method_dirs[!method_dirs %in% c('enrichment')]

  pvals <- list()
  for (m in method_dirs) {
    f <- file.path(base, m, paste0(m, '_pvalues_K', k, '.tsv'))
    pvals[[m]] <- safe_fread(f)
  }

  list(
    pvals = pvals,
    methods = method_dirs,
    selected = safe_fread(file.path(base, 'Selected_SNPs.tsv')),
    regions_combined = safe_fread(file.path(base, 'Regions_climate_combined.tsv')),
    regions_trait = safe_fread(file.path(base, 'Regions_per_trait.tsv')),
    genes = safe_fread(file.path(base, 'Genes_per_region.tsv')),
    genes_collapsed = safe_fread(file.path(base, 'Genes_per_region_collapsed.tsv'))
  )
}

# Prepare Manhattan plot data: reshape to long format, compute cumulative positions
prep_manhattan <- function(pvals_list, selected_snps = NULL) {
  all_data <- data.table::rbindlist(lapply(names(pvals_list), function(method) {
    dt <- pvals_list[[method]]
    if (is.null(dt) || nrow(dt) == 0) return(NULL)

    trait_cols <- setdiff(names(dt), c('SNPID', 'chr', 'pos'))
    if (length(trait_cols) == 0) return(NULL)

    melted <- data.table::melt(dt, id.vars = c('SNPID', 'chr', 'pos'),
                               measure.vars = trait_cols,
                               variable.name = 'trait',
                               value.name = 'pvalue')
    melted[, method := method]
    melted
  }))

  if (is.null(all_data) || nrow(all_data) == 0) return(NULL)

  all_data <- all_data[!is.na(pvalue) & pvalue > 0]
  all_data[, neglog10p := -log10(pvalue)]
  all_data[, series := paste0(trait, ' (', method, ')')]

  # Sort chromosomes numerically
  all_data[, chr_num := suppressWarnings(as.integer(gsub('[^0-9]', '', as.character(chr))))]
  chr_order <- unique(all_data[order(chr_num, chr), chr])
  all_data[, chr := factor(chr, levels = chr_order)]

  # Cumulative positions with gaps between chromosomes
  chr_info <- all_data[, .(max_pos = max(pos)), by = chr]
  chr_info <- chr_info[match(chr_order, chr)]
  chr_gap <- max(chr_info$max_pos) * 0.03
  chr_info[, cum_offset := cumsum(c(0, max_pos[-.N] + chr_gap))]
  chr_info[, midpoint := cum_offset + max_pos / 2]

  all_data <- merge(all_data, chr_info[, .(chr, cum_offset)], by = 'chr')
  all_data[, pos_cum := pos + cum_offset]

  # Mark significant SNPs per trait+method (not just by SNPID)
  if (!is.null(selected_snps) && nrow(selected_snps) > 0) {
    # Build sig lookup: (SNPID, trait, method) combinations
    method_cols <- intersect(names(selected_snps),
                             setdiff(names(selected_snps),
                                     c('SNPID', 'chr', 'pos', 'min_pvalue')))
    sig_keys <- character(0)
    for (m in method_cols) {
      ss <- selected_snps[, .(SNPID, traits = get(m))]
      ss <- ss[!is.na(traits) & traits != "" & traits != '""']
      if (nrow(ss) > 0) {
        expanded <- ss[, .(trait = unlist(strsplit(as.character(traits), ","))),
                       by = SNPID]
        sig_keys <- c(sig_keys,
                      paste0(expanded$SNPID, '_', expanded$trait, ' (', m, ')'))
      }
    }
    all_data[, is_sig := paste0(SNPID, '_', series) %in% sig_keys]
  } else {
    all_data[, is_sig := FALSE]
  }

  list(data = all_data, chr_info = chr_info)
}

# Generate distinct colors per trait
generate_trait_colors <- function(trait_names) {
  n <- length(trait_names)
  if (n <= 8) {
    pal <- RColorBrewer::brewer.pal(max(3, n), 'Set1')[1:n]
  } else if (n <= 12) {
    pal <- RColorBrewer::brewer.pal(n, 'Set3')
  } else {
    pal <- viridis::turbo(n)
  }
  setNames(pal, trait_names)
}

# Map methods to plotly marker symbols
generate_method_symbols <- function(method_names) {
  symbols <- c('circle', 'triangle-up', 'square', 'diamond',
               'cross', 'x', 'star', 'hexagon')
  setNames(symbols[seq_along(method_names)], method_names)
}

# Read a config value for a project (searches config*.yaml files)
read_project_config <- function(project, key, pipeline_path = '/pipeline/') {
  configs <- list.files(pipeline_path, pattern = '^config.*\\.ya?ml$', full.names = TRUE)
  for (cf in configs) {
    lines <- readLines(cf, warn = FALSE)
    proj_match <- grep(paste0('PROJECT_NAME:\\s*["\']?', project, '["\']?'), lines)
    if (length(proj_match) > 0) {
      val_line <- grep(paste0('^', key, ':'), lines, value = TRUE)
      if (length(val_line) > 0) {
        val <- trimws(sub(paste0('^', key, ':\\s*'), '', val_line[1]))
        val <- trimws(sub('#.*$', '', val))  # Strip YAML comments
        return(gsub('["\']', '', val))
      }
    }
  }
  return(NULL)
}

# Load enrichment data for a region (handles both combined and per-trait region IDs)
load_region_enrichment <- function(project, rid, pipeline_path = '/pipeline/') {
  base <- paste0(pipeline_path, project, '_results/tables/association/enrichment')
  if (!dir.exists(base)) return(NULL)

  trait_dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)

  results <- data.table::rbindlist(lapply(trait_dirs, function(tdir) {
    trait_name <- basename(tdir)

    # Exact match (per-trait region_id already includes trait suffix)
    f <- file.path(tdir, paste0('Region_', rid, '_enrichment.tsv'))
    if (file.exists(f)) return(data.table::fread(f))

    # Combined region_id: add trait suffix
    f2 <- file.path(tdir, paste0('Region_', rid, '_', trait_name, '_enrichment.tsv'))
    if (file.exists(f2)) return(data.table::fread(f2))

    NULL
  }), fill = TRUE)

  if (nrow(results) > 0) return(results)
  NULL
}

###############################################################################
# PHENOTYPE ASSOCIATION HELPER FUNCTIONS
###############################################################################

# Find K value(s) from phenotype association output files
find_pheno_assoc_k <- function(project, pipeline_path = '/pipeline/') {
  base <- paste0(pipeline_path, project, '_results/tables/association_phenotypes')
  if (!dir.exists(base)) return(NULL)

  # Combined pvalues files are at parent level: {METHOD}_phenotypes_pvalues_K{K}.tsv
  files <- list.files(base, pattern = '_phenotypes_pvalues_K\\d+\\.tsv$')
  if (length(files) == 0) return(NULL)

  k <- unique(as.integer(stringr::str_extract(files, '(?<=_K)\\d+')))
  sort(k[!is.na(k)])
}

# Load all phenotype association data for a project
load_pheno_assoc_data <- function(project, k, pipeline_path = '/pipeline/') {
  base <- paste0(pipeline_path, project, '_results/tables/association_phenotypes')

  safe_fread <- function(path) {
    if (file.exists(path)) data.table::fread(path) else NULL
  }

  # Combined pvalues are at parent level: {METHOD}_phenotypes_pvalues_K{K}.tsv
  pval_files <- list.files(base, pattern = paste0('_phenotypes_pvalues_K', k, '\\.tsv$'))
  methods <- stringr::str_extract(pval_files, '^[A-Z]+')

  pvals <- list()
  for (i in seq_along(methods)) {
    pvals[[methods[i]]] <- safe_fread(file.path(base, pval_files[i]))
  }

  list(
    pvals = pvals,
    methods = methods,
    selected = safe_fread(file.path(base, 'Selected_SNPs.tsv')),
    regions_combined = safe_fread(file.path(base, 'Regions_phenotype_combined.tsv')),
    regions_trait = safe_fread(file.path(base, 'Regions_per_trait.tsv')),
    genes = safe_fread(file.path(base, 'Genes_per_region.tsv')),
    genes_collapsed = safe_fread(file.path(base, 'Genes_per_region_collapsed.tsv')),
    genes_combined = safe_fread(file.path(base, 'Genes_per_region_combined.tsv'))
  )
}

# Load enrichment data for a phenotype association region
load_region_enrichment_pheno <- function(project, rid, pipeline_path = '/pipeline/') {
  base <- paste0(pipeline_path, project, '_results/tables/association_phenotypes/enrichment')
  if (!dir.exists(base)) return(NULL)

  trait_dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)

  results <- data.table::rbindlist(lapply(trait_dirs, function(tdir) {
    trait_name <- basename(tdir)

    # Exact match (per-trait region_id already includes trait suffix)
    f <- file.path(tdir, paste0('Region_', rid, '_enrichment.tsv'))
    if (file.exists(f)) return(data.table::fread(f))

    # Combined region_id: add trait suffix
    f2 <- file.path(tdir, paste0('Region_', rid, '_', trait_name, '_enrichment.tsv'))
    if (file.exists(f2)) return(data.table::fread(f2))

    NULL
  }), fill = TRUE)

  if (nrow(results) > 0) return(results)
  NULL
}

###############################################################################
# OVERLAPPING ANALYSIS HELPER FUNCTIONS
###############################################################################

# Load overlap analysis data
load_overlap_data <- function(project, pipeline_path = '/pipeline/') {
  base <- paste0(pipeline_path, project, '_results/tables/overlapping')

  safe_fread <- function(path) {
    if (file.exists(path)) data.table::fread(path) else NULL
  }

  list(
    selected = safe_fread(file.path(base, 'Selected_SNPs_all.tsv')),
    regions_combined = safe_fread(file.path(base, 'Regions_all_combined.tsv')),
    regions_trait = safe_fread(file.path(base, 'Regions_per_trait_all.tsv')),
    overlap_summary = safe_fread(file.path(base, 'Overlap_summary.tsv')),
    genes = safe_fread(file.path(base, 'Genes_per_region.tsv')),
    genes_collapsed = safe_fread(file.path(base, 'Genes_per_region_collapsed.tsv')),
    genes_combined = safe_fread(file.path(base, 'Genes_per_region_combined.tsv'))
  )
}

# Load enrichment data for an overlap region
load_region_enrichment_overlap <- function(project, rid, pipeline_path = '/pipeline/') {
  base <- paste0(pipeline_path, project, '_results/tables/overlapping/enrichment')
  if (!dir.exists(base)) return(NULL)

  trait_dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)

  results <- data.table::rbindlist(lapply(trait_dirs, function(tdir) {
    trait_name <- basename(tdir)

    f <- file.path(tdir, paste0('Region_', rid, '_enrichment.tsv'))
    if (file.exists(f)) return(data.table::fread(f))

    f2 <- file.path(tdir, paste0('Region_', rid, '_', trait_name, '_enrichment.tsv'))
    if (file.exists(f2)) return(data.table::fread(f2))

    NULL
  }), fill = TRUE)

  if (nrow(results) > 0) return(results)
  NULL
}

# Function to load plot from .qs file (plotname is a regex pattern)
load_plot <- function(project, plotname, pipeline_path = '/pipeline/') {
  tryCatch({
    plot_paths <- list.files(paste0(pipeline_path, project, '_results/intermediate'),
                             pattern = plotname, full.names = TRUE)
    if (length(plot_paths) == 0) {
      stop('NO FILE matching: ', plotname)
    }
    plot_path <- plot_paths[1]  # Use first match
    message("Loading: ", plot_path)
    qread(plot_path)
  }, error = function(e) {
    warning(paste("Error loading", plotname, "for project", project, ":", e$message))
    return(NULL)
  })
}

render_plot <- function(local_proj, plotname){

	renderPlot({
		   gPlot <- load_plot(local_proj, plotname)
		   if (!is.null(gPlot)) {
			gPlot
		   } else {
		   # Return an empty plot if plot couldn't be loaded
		    ggplot() +
			theme_void() +
			theme(plot.background = element_rect(fill = "white"))
			   }
			          }, res = 96)
} # end render_plot

render_noplot_message <- function(local_proj, plotname){ #TODO further add the name of the plot(extract from local_proj)
	renderUI({
          gPlot <- load_plot(local_proj, plotname)

          if (is.null(gPlot)) {
            div(
              style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);",
              h4("No plot available for this project", style = "color: #666;")
            )
          }
        })
}
# UI
ui <- dashboardPage(
  dashboardHeader(title = "Pipeline Monitor"),
  
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      menuItem("Home", tabName = "home", icon = icon("home")),
      tags$div(id = "project_menu")  # Placeholder for dynamic menu items
    )
  ),
  
  dashboardBody(
    tabItems(
      # Home tab
      tabItem(
        tabName = "home",
        h2("Welcome to Pipeline Monitor"),
        p("Select a project from the sidebar to begin.")
      ),
      
      # Dynamic project content will be rendered here
      uiOutput("dynamic_content")
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive expression to get projects
  projects <- reactive({
    find_projects()
  })
  
  # Update sidebar menu with projects
  observe({
    proj_list <- projects()
    
    # Remove existing project menu items
    removeUI(selector = "#project_menu > *")
    
    # Add new project menu items
    for (proj in proj_list) {
      insertUI(
        selector = "#project_menu",
        ui = menuItem(proj, tabName = paste0("project_", proj), icon = icon("folder"))
      )
    }
  })
  
  # Create dynamic content based on selected tab
  output$dynamic_content <- renderUI({
    # Get current selected tab
    current_tab <- input$sidebar
    
    # If it's a project tab
    if (grepl("^project_", current_tab)) {
      proj <- sub("^project_", "", current_tab)
      k_values <- find_k_values(proj)
      bio_values <- find_bio_values(proj)
      metric_values <- find_metric_values(proj)
      assoc_k <- find_assoc_k(proj)
      has_assoc <- !is.null(assoc_k) && length(assoc_k) > 0
      gf_suffixes <- find_gf_suffixes(proj)
      has_malad <- length(gf_suffixes) > 0

      # Read top N regions default from project config
      assoc_top_n_default <- 10L
      if (has_assoc) {
        cfg_val <- read_project_config(proj, 'ASSOC_TOP_REGIONS')
        if (!is.null(cfg_val)) {
          assoc_top_n_default <- suppressWarnings(as.integer(cfg_val))
          if (is.na(assoc_top_n_default)) assoc_top_n_default <- 10L
        }
      }

      # Phenotype association detection
      pheno_assoc_k <- find_pheno_assoc_k(proj)
      has_pheno_assoc <- !is.null(pheno_assoc_k) && length(pheno_assoc_k) > 0

      # Read PHENO config defaults
      pheno_top_n_default <- 10L
      if (has_pheno_assoc) {
        cfg_val <- read_project_config(proj, 'PHENO_TOP_REGIONS')
        if (!is.null(cfg_val)) {
          pheno_top_n_default <- suppressWarnings(as.integer(cfg_val))
          if (is.na(pheno_top_n_default)) pheno_top_n_default <- 10L
        }
      }

      # Overlapping analysis detection
      has_overlap <- has_assoc && has_pheno_assoc &&
        file.exists(paste0('/pipeline/', proj, '_results/tables/overlapping/Overlap_summary.tsv'))

      # Read OVERLAP config defaults
      ov_top_n_default <- 10L
      if (has_overlap) {
        cfg_val <- read_project_config(proj, 'OVERLAP_TOP_REGIONS')
        if (!is.null(cfg_val)) {
          ov_top_n_default <- suppressWarnings(as.integer(cfg_val))
          if (is.na(ov_top_n_default)) ov_top_n_default <- 10L
        }
      }

      tabItem(
        tabName = paste0("project_", proj),
        fluidRow(
          box(
            width = 12,
            title = paste("Project:", proj),
            status = "primary",
            solidHeader = TRUE
          )
        ),
        fluidRow(
          box(
            width = 12,
            tabBox(
              id = paste0("tabbox_", proj),
              width = NULL,
              tabPanel(
                "Structure",
                h3("Structure Analysis"),
			fluidRow(
				 column(
					width = 4,
			                div(
		                        	style = "position: relative; min-height: 200px;",
	         	     			plotOutput(paste0("pcaTW_plot_", proj)),
			                	uiOutput(paste0("pcaTW_message_", proj))
			               	   )
					),
				column(
				       width = 4,
					div(
                                                style = "position: relative; min-height: 200px;",
                                                plotOutput(paste0("CrossEntropy_plot_", proj)),
                                                uiOutput(paste0("CrossEntropy_message_", proj))
                                           )

					),
				 column(
					width = 4,
					div(
			                  	style = "position: relative; min-height: 200px;",
                 			 	plotOutput(paste0("PCA_plot_", proj)),
                  				uiOutput(paste0("PCA_message_", proj))
                			   )
					)
              			),
		        # New section for K-dependent plots
                fluidRow(
                  box(
                    width = 12,
                    title = "K-dependent Analysis",
                    sliderInput(
	                      paste0("k_slider_", proj),
	                      "Select K value:",
	                      min = min(k_values),
	                      max = max(k_values),
	                      value = min(k_values),
	                      step = 1
        	               ),
                    fluidRow(
                      column(
                        width = 6,
                        div(
                          style = "position: relative; min-height: 200px;",
                          plotOutput(paste0("structurePie_plot_", proj)),
                          uiOutput(paste0("structurePie_message_", proj))
                        )
                      ),
                      column(
                        width = 6,
                        div(
                          style = "position: relative; min-height: 200px;",
                          plotOutput(paste0("popDiffTest_plot_", proj)),
                          uiOutput(paste0("popDiffTest_message_", proj))
                        )
                      )
                    ),
		      fluidRow(
                        width = 12,
                        div(
                          style = "position: relative; min-height: 300px;",
                          plotOutput(paste0("structure_plot_", proj)),
                          uiOutput(paste0("structure_message_", proj))
                        )
                      ),
                  )
                )
		      ),
###############################################################################
              tabPanel(
                "Structure K",
                h3("Structure K Analysis"),
                fluidRow(
                         column(
                                width = 12,
                                div(
                                    style = "position: relative; min-height: 300px;",
                                    plotOutput(paste0("CorHM_plot_", proj)),
                                    uiOutput(paste0("CorHM_message_", proj))
                                   )
                               ) # add new column here
			), # add new rows here
		fluidRow(
			 box(
			     width = 12,
			     title = 'Bio dependent analysis',
			     selectInput(
       				 inputId = paste0("bio_choose_", proj),
			         label = "Choose a bio:",
  			         choices = bio_values,
			         selected = bio_values[1]
      					),
			     selectInput(
                                 inputId = paste0("metric_choose_", proj),
                                 label = "Choose a metric used for control pie chart size:",
                                 choices = metric_values,
                                 selected = metric_values[1]
                                        ),
				     fluidRow(
					column(
                                                width = 4,
                                                div(
                                                   style = "position: relative; min-height: 200px;",
                                                   plotOutput(paste0("BioDensity_plot_", proj)),
                                                   uiOutput(paste0("BioDensity_message_", proj))
                                                    )
                                              ),
                      			column(
                        			width = 4,
                        			div(
			                           style = "position: relative; min-height: 200px;",
			                           plotOutput(paste0("PieMap_plot_", proj)),
                        			   uiOutput(paste0("PieMap_message_", proj))
                           			    )
                      	    		      )
                    		     	     )	
			    )
			)
	              ),
              tabPanel(
                "Association",
                h3("Association Analysis"),
                if (!has_assoc) {
                  div(style = "text-align: center; padding: 60px; color: #888;",
                    icon("bar-chart"), h4("Association analysis not yet run."),
                    p("Run the pipeline with mode=association to generate results."))
                } else { tagList(
                  # Summary value boxes
                  fluidRow(
                    valueBoxOutput(paste0("assoc_vb_snps_unique_", proj), width = 2),
                    valueBoxOutput(paste0("assoc_vb_snps_hits_", proj), width = 2),
                    valueBoxOutput(paste0("assoc_vb_regions_", proj), width = 2),
                    valueBoxOutput(paste0("assoc_vb_genes_", proj), width = 2),
                    valueBoxOutput(paste0("assoc_vb_enrich_", proj), width = 2),
                    valueBoxOutput(paste0("assoc_vb_methods_", proj), width = 2)
                  ),
                  # Controls
                  fluidRow(
                    box(width = 12, title = "Controls", collapsible = TRUE,
                        status = "primary", solidHeader = TRUE,
                      fluidRow(
                        column(4,
                          h5("Trait x Method layers"),
                          uiOutput(paste0("assoc_layer_ui_", proj)),
                          fluidRow(
                            column(6, actionButton(paste0("assoc_all_", proj), "All",
                                                   class = "btn-xs btn-default", width = "100%")),
                            column(6, actionButton(paste0("assoc_none_", proj), "None",
                                                   class = "btn-xs btn-default", width = "100%"))
                          )
                        ),
                        column(4,
                          h5("Regions"),
                          checkboxInput(paste0("assoc_show_regions_", proj),
                                        "Show top regions", value = TRUE),
                          sliderInput(paste0("assoc_top_n_", proj),
                                      "Top N regions:", min = 1,
                                      max = assoc_top_n_default,
                                      value = assoc_top_n_default, step = 1),
                          radioButtons(paste0("assoc_region_type_", proj),
                                       "Region type:",
                                       choices = c("Combined" = "combined",
                                                   "Per-trait" = "per_trait"),
                                       selected = "combined", inline = TRUE)
                        ),
                        column(4,
                          h5("Region Selection"),
                          selectInput(paste0("assoc_region_select_", proj),
                                      "Jump to region:",
                                      choices = c("(click SNP or select)" = ""),
                                      width = "100%"),
                          actionButton(paste0("assoc_clear_region_", proj),
                                       "Clear selection", class = "btn-xs btn-default",
                                       icon = icon("times"))
                        )
                      )
                    )
                  ),
                  # Manhattan plot
                  fluidRow(
                    box(width = 12, title = "Interactive Manhattan Plot",
                        status = "success", solidHeader = TRUE,
                      plotlyOutput(paste0("assoc_manhattan_", proj), height = "550px")
                    )
                  ),
                  # Region info bar
                  uiOutput(paste0("assoc_region_info_", proj)),
                  # Gene and Enrichment tables
                  fluidRow(
                    box(width = 6, title = "Genes in Selected Region",
                        status = "info", solidHeader = TRUE,
                      div(style = "min-height: 280px;",
                        uiOutput(paste0("assoc_genes_prompt_", proj)),
                        DT::dataTableOutput(paste0("assoc_genes_dt_", proj))
                      )
                    ),
                    box(width = 6, title = "GO Enrichment",
                        status = "warning", solidHeader = TRUE,
                      div(style = "min-height: 280px;",
                        uiOutput(paste0("assoc_enrich_prompt_", proj)),
                        DT::dataTableOutput(paste0("assoc_enrich_dt_", proj))
                      )
                    )
                  )
                )} # end tagList
              ),
###############################################################################
# PHENOTYPE ASSOCIATION TAB
###############################################################################
              tabPanel(
                "Phenotype Association",
                h3("Phenotype Association Analysis (GWAS)"),
                if (!has_pheno_assoc) {
                  div(style = "text-align: center; padding: 60px; color: #888;",
                    icon("bar-chart"), h4("Phenotype association analysis not yet run."),
                    p("Run the pipeline with mode=association_phenotypes to generate results."))
                } else { tagList(
                  # Summary value boxes
                  fluidRow(
                    valueBoxOutput(paste0("pheno_vb_snps_unique_", proj), width = 2),
                    valueBoxOutput(paste0("pheno_vb_snps_hits_", proj), width = 2),
                    valueBoxOutput(paste0("pheno_vb_regions_", proj), width = 2),
                    valueBoxOutput(paste0("pheno_vb_genes_", proj), width = 2),
                    valueBoxOutput(paste0("pheno_vb_enrich_", proj), width = 2),
                    valueBoxOutput(paste0("pheno_vb_methods_", proj), width = 2)
                  ),
                  # Controls
                  fluidRow(
                    box(width = 12, title = "Controls", collapsible = TRUE,
                        status = "primary", solidHeader = TRUE,
                      fluidRow(
                        column(3,
                          h5("Trait Selector"),
                          uiOutput(paste0("pheno_trait_selector_ui_", proj))
                        ),
                        column(3,
                          h5("Trait x Method layers"),
                          uiOutput(paste0("pheno_layer_ui_", proj)),
                          fluidRow(
                            column(6, actionButton(paste0("pheno_all_", proj), "All",
                                                   class = "btn-xs btn-default", width = "100%")),
                            column(6, actionButton(paste0("pheno_none_", proj), "None",
                                                   class = "btn-xs btn-default", width = "100%"))
                          )
                        ),
                        column(3,
                          h5("Regions"),
                          checkboxInput(paste0("pheno_show_regions_", proj),
                                        "Show top regions", value = TRUE),
                          sliderInput(paste0("pheno_top_n_", proj),
                                      "Top N regions:", min = 1,
                                      max = pheno_top_n_default,
                                      value = pheno_top_n_default, step = 1),
                          radioButtons(paste0("pheno_region_type_", proj),
                                       "Region type:",
                                       choices = c("Combined" = "combined",
                                                   "Per-trait" = "per_trait"),
                                       selected = "combined", inline = TRUE)
                        ),
                        column(3,
                          h5("Region Selection"),
                          selectInput(paste0("pheno_region_select_", proj),
                                      "Jump to region:",
                                      choices = c("(click SNP or select)" = ""),
                                      width = "100%"),
                          actionButton(paste0("pheno_clear_region_", proj),
                                       "Clear selection", class = "btn-xs btn-default",
                                       icon = icon("times"))
                        )
                      )
                    )
                  ),
                  # Trait exploration row: distribution + site map
                  fluidRow(
                    box(width = 6, title = "Trait Distribution",
                        status = "info", solidHeader = TRUE,
                      plotlyOutput(paste0("pheno_trait_dist_", proj), height = "350px")
                    ),
                    box(width = 6, title = "Trait PieMap",
                        status = "info", solidHeader = TRUE,
                      div(style = "position: relative; min-height: 350px;",
                        plotOutput(paste0("pheno_site_map_", proj), height = "450px")
                      )
                    )
                  ),
                  # Manhattan plot
                  fluidRow(
                    box(width = 12, title = "Interactive Manhattan Plot (Phenotype GWAS)",
                        status = "success", solidHeader = TRUE,
                      plotlyOutput(paste0("pheno_manhattan_", proj), height = "550px")
                    )
                  ),
                  # Region info bar
                  uiOutput(paste0("pheno_region_info_", proj)),
                  # Gene and Enrichment tables
                  fluidRow(
                    box(width = 6, title = "Genes in Selected Region",
                        status = "info", solidHeader = TRUE,
                      div(style = "min-height: 280px;",
                        uiOutput(paste0("pheno_genes_prompt_", proj)),
                        DT::dataTableOutput(paste0("pheno_genes_dt_", proj))
                      )
                    ),
                    box(width = 6, title = "GO Enrichment",
                        status = "warning", solidHeader = TRUE,
                      div(style = "min-height: 280px;",
                        uiOutput(paste0("pheno_enrich_prompt_", proj)),
                        DT::dataTableOutput(paste0("pheno_enrich_dt_", proj))
                      )
                    )
                  )
                )} # end tagList
              ),
###############################################################################
# HAPLOTYPE ANALYSIS TAB (PLACEHOLDER)
###############################################################################
              tabPanel(
                "Haplotype Analysis",
                h3("Haplotype Analysis"),
                div(style = "text-align: center; padding: 60px; color: #999;",
                  icon("dna", "fa-3x"),
                  h4("Coming Soon"),
                  p("Haplotype analysis will be available in a future update.")
                )
              ),
###############################################################################
# OVERLAPPING REGIONS TAB
###############################################################################
              tabPanel(
                "Overlapping Regions",
                h3("Overlapping Regions (GEA vs GWAS)"),
                if (!has_overlap) {
                  div(style = "text-align: center; padding: 60px; color: #888;",
                    icon("exchange"), h4("Overlapping analysis not yet run."),
                    p("Run mode=association, mode=association_phenotypes, and mode=overlapping."))
                } else { tagList(
                  # Summary boxes
                  fluidRow(
                    valueBoxOutput(paste0("ov_box_gea_", proj), width = 2),
                    valueBoxOutput(paste0("ov_box_gwas_", proj), width = 2),
                    valueBoxOutput(paste0("ov_box_pairs_", proj), width = 2),
                    valueBoxOutput(paste0("ov_box_merged_", proj), width = 2),
                    valueBoxOutput(paste0("ov_box_genes_", proj), width = 2),
                    valueBoxOutput(paste0("ov_box_go_", proj), width = 2)
                  ),
                  # Controls (same pattern as Association tab)
                  fluidRow(
                    box(width = 12, title = "Controls", collapsible = TRUE,
                        status = "primary", solidHeader = TRUE,
                      fluidRow(
                        column(4,
                          h5("Trait x Method layers"),
                          uiOutput(paste0("ov_layer_ui_", proj)),
                          fluidRow(
                            column(6, actionButton(paste0("ov_all_", proj), "All",
                                                   class = "btn-xs btn-default", width = "100%")),
                            column(6, actionButton(paste0("ov_none_", proj), "None",
                                                   class = "btn-xs btn-default", width = "100%"))
                          )
                        ),
                        column(4,
                          h5("Regions"),
                          checkboxInput(paste0("ov_show_regions_", proj),
                                        "Show top regions", value = TRUE),
                          sliderInput(paste0("ov_top_n_", proj),
                                      "Top N regions:", min = 1,
                                      max = ov_top_n_default,
                                      value = ov_top_n_default, step = 1),
                          radioButtons(paste0("ov_region_type_", proj),
                                       "Region type:",
                                       choices = c("Combined" = "combined",
                                                   "Per-trait" = "per_trait"),
                                       selected = "combined", inline = TRUE)
                        ),
                        column(4,
                          h5("Region Selection"),
                          selectInput(paste0("ov_region_select_", proj),
                                      "Jump to region:",
                                      choices = c("(click SNP or select)" = ""),
                                      width = "100%"),
                          actionButton(paste0("ov_clear_region_", proj),
                                       "Clear selection", class = "btn-xs btn-default",
                                       icon = icon("times"))
                        )
                      )
                    )
                  ),
                  # Miami Manhattan plot
                  fluidRow(
                    box(width = 12, title = "Miami Manhattan Plot (GEA top / GWAS bottom)",
                        status = "success", solidHeader = TRUE,
                      plotlyOutput(paste0("ov_miami_", proj), height = "550px")
                    )
                  ),
                  # Region info bar
                  uiOutput(paste0("ov_region_info_", proj)),
                  # Genes + Enrichment tables
                  fluidRow(
                    box(width = 6, title = "Genes in Selected Region",
                        status = "info", solidHeader = TRUE,
                      div(style = "min-height: 280px;",
                        uiOutput(paste0("ov_genes_prompt_", proj)),
                        DT::dataTableOutput(paste0("ov_genes_dt_", proj))
                      )
                    ),
                    box(width = 6, title = "GO Enrichment",
                        status = "warning", solidHeader = TRUE,
                      div(style = "min-height: 280px;",
                        uiOutput(paste0("ov_enrich_prompt_", proj)),
                        DT::dataTableOutput(paste0("ov_enrich_dt_", proj))
                      )
                    )
                  ),
                  # Overlap pairs table
                  fluidRow(
                    box(width = 12, title = "Overlap Pairs (GEA vs GWAS regions)",
                        status = "success", solidHeader = TRUE, collapsible = TRUE,
                      div(style = "min-height: 200px;",
                        DT::dataTableOutput(paste0("ov_overlap_dt_", proj))
                      )
                    )
                  )
                )} # end tagList
              ),
###############################################################################
# MALADAPTATION TAB
###############################################################################
              tabPanel(
                "Maladaptation",
                h3("Maladaptation Analysis"),
                if (!has_malad) {
                  div(style = "text-align: center; padding: 60px; color: #888;",
                    icon("leaf"), h4("Maladaptation analysis not yet run."),
                    p("Run the pipeline with mode=maladaptation to generate results."))
                } else { tagList(
                  # GF suffix selector (PCNM variant)
                  fluidRow(
                    box(width = 12, title = "Model Selection", collapsible = TRUE,
                        status = "primary", solidHeader = TRUE,
                      fluidRow(
                        column(4,
                          selectInput(paste0("gf_suffix_", proj),
                                      "Gradient Forest model:",
                                      choices = gf_suffixes,
                                      selected = gf_suffixes[1])
                        ),
                        column(8,
                          uiOutput(paste0("gf_site_info_", proj))
                        )
                      )
                    )
                  ),
                  # Row 1: Importance plots
                  fluidRow(
                    box(width = 6, title = "Overall Importance",
                        status = "info", solidHeader = TRUE,
                      div(style = "position: relative; min-height: 300px;",
                        plotOutput(paste0("gf_importance_", proj)),
                        uiOutput(paste0("gf_importance_msg_", proj))
                      )
                    ),
                    box(width = 6, title = "Cumulative Importance",
                        status = "info", solidHeader = TRUE,
                      div(style = "position: relative; min-height: 300px;",
                        plotOutput(paste0("gf_cumimp_", proj)),
                        uiOutput(paste0("gf_cumimp_msg_", proj))
                      )
                    )
                  ),
                  # Row 2: Genetic Offset PieMaps
                  fluidRow(
                    box(width = 4, title = "Genetic Offset Map",
                        status = "success", solidHeader = TRUE,
                      div(style = "position: relative; min-height: 350px;",
                        plotOutput(paste0("gf_piemap_", proj)),
                        uiOutput(paste0("gf_piemap_msg_", proj))
                      )
                    ),
                    box(width = 4, title = "Offset + Tajima's D",
                        status = "success", solidHeader = TRUE,
                      div(style = "position: relative; min-height: 350px;",
                        plotOutput(paste0("gf_piemap_tajima_", proj)),
                        uiOutput(paste0("gf_piemap_tajima_msg_", proj))
                      )
                    ),
                    box(width = 4, title = "Offset + Pi Diversity",
                        status = "success", solidHeader = TRUE,
                      div(style = "position: relative; min-height: 350px;",
                        plotOutput(paste0("gf_piemap_diversity_", proj)),
                        uiOutput(paste0("gf_piemap_diversity_msg_", proj))
                      )
                    )
                  ),
                  # Row 3: Zoom maps (if available)
                  uiOutput(paste0("gf_zoom_ui_", proj)),
                  # Row 4: Genetic Offset site values table
                  fluidRow(
                    box(width = 12, title = "Genetic Offset per Site",
                        status = "warning", solidHeader = TRUE, collapsible = TRUE,
                      DT::dataTableOutput(paste0("gf_site_dt_", proj))
                    )
                  )
                )} # end tagList
              )
            )
          )
        )
      )
    }
  })
  
  # Create plot outputs for the currently selected project only
  observe({
    current_tab <- input$sidebar
    
    if (grepl("^project_", current_tab)) {
      proj <- sub("^project_", "", current_tab)
      
      local({
        local_proj <- proj
        
	# PCA
	output[[paste0("PCA_plot_", local_proj)]] <- render_plot(local_proj, '^pca\\.qs$')
        output[[paste0("PCA_message_", local_proj)]] <- render_noplot_message(local_proj, '^pca\\.qs$')

        # pcaTW
        output[[paste0("pcaTW_plot_", local_proj)]] <- render_plot(local_proj, 'tracy_widom\\.qs$')
        output[[paste0("pcaTW_message_", local_proj)]] <- render_noplot_message(local_proj, 'tracy_widom\\.qs$')

	#   	Cross_entropy_K
	output[[paste0("CrossEntropy_plot_", local_proj)]] <- render_plot(local_proj, 'cross_entropy_K[0-9]*-[0-9]*\\.qs$')
        output[[paste0("CrossEntropy_message_", local_proj)]] <- render_noplot_message(local_proj, 'cross_entropy_K[0-9]*-[0-9]*\\.qs$')
	
	 # New reactive plot outputs for K-dependent plots
        observeEvent(input[[paste0("k_slider_", local_proj)]], {
          k_value <- input[[paste0("k_slider_", local_proj)]]

          # Structure Pie plot (PCA + sNMF ancestry pie)
          output[[paste0("structurePie_plot_", local_proj)]] <-
            render_plot(local_proj, paste0('pca_structure_K', k_value, '\\.qs$'))
          output[[paste0("structurePie_message_", local_proj)]] <-
            render_noplot_message(local_proj, paste0('pca_structure_K', k_value, '\\.qs$'))

          # Population Difference Test plot
          output[[paste0("popDiffTest_plot_", local_proj)]] <-
            render_plot(local_proj, paste0('pop_diff_K', k_value, '\\.qs$'))
          output[[paste0("popDiffTest_message_", local_proj)]] <-
            render_noplot_message(local_proj, paste0('pop_diff_K', k_value, '\\.qs$'))

          # Structure plot (SNMF barplot)
          output[[paste0("structure_plot_", local_proj)]] <-
            render_plot(local_proj, paste0('^structure_K', k_value, '\\.qs$'))
          output[[paste0("structure_message_", local_proj)]] <-
            render_noplot_message(local_proj, paste0('^structure_K', k_value, '\\.qs$'))
        })

################################################################################################
	# STRUCTURE K

        # CorHM
	output[[paste0("CorHM_plot_", local_proj)]] <- render_plot(local_proj, 'CorrelationHeatmap.*\\.qs$')
        output[[paste0("CorHM_message_", local_proj)]] <- render_noplot_message(local_proj, 'CorrelationHeatmap.*\\.qs$')

	# PieMaps + BioDensity
	observeEvent({
		     input[[paste0("bio_choose_", local_proj)]]
		     input[[paste0("metric_choose_", local_proj)]]
		     }, {
          bio_value <- input[[paste0("bio_choose_", local_proj)]]
	  metric_value <- input[[paste0("metric_choose_", local_proj)]]
	  # Density Plot
	  # DensityPlot_bio_5.qs
	  output[[paste0("BioDensity_plot_", local_proj)]] <-
            render_plot(local_proj, paste0('DensityPlot_bio_', bio_value, '.qs')) 
          output[[paste0("BioDensity_message_", local_proj)]] <-
            render_noplot_message(local_proj, paste0('DensityPlot_bio_', bio_value, '.qs')) 

          # Structure Pie plot Tajima
          output[[paste0("PieMap_plot_", local_proj)]] <-
            render_plot(local_proj, paste0('PieMap_bio_', bio_value, '_', metric_value, '.qs'))
          output[[paste0("PieMap_message_", local_proj)]] <-
            render_noplot_message(local_proj, paste0('PieMap_bio_', bio_value, '_', metric_value, '.qs'))

        })

        
        ###############################################################
        # ASSOCIATION TAB SERVER LOGIC
        ###############################################################
        assoc_k_val <- find_assoc_k(local_proj)
        has_assoc_data <- !is.null(assoc_k_val) && length(assoc_k_val) > 0

        # Initialize to NULL (needed by Overlapping Regions tab)
        ad <- NULL; mpd <- NULL

        if (has_assoc_data) {
          # Load data eagerly (pipeline outputs don't change during session)
          ad <- load_assoc_data(local_proj, assoc_k_val[1])
          mpd <- prep_manhattan(ad$pvals, ad$selected)
          series_names <- if (!is.null(mpd)) sort(unique(mpd$data$series)) else character(0)
          trait_names <- if (!is.null(mpd)) sort(unique(as.character(mpd$data$trait))) else character(0)
          method_names <- if (!is.null(mpd)) sort(unique(mpd$data$method)) else character(0)
          t_colors <- if (length(trait_names) > 0) generate_trait_colors(trait_names) else character(0)
          m_symbols <- if (length(method_names) > 0) generate_method_symbols(method_names) else character(0)

          # Reactive value for selected region
          sel_region <- reactiveVal(NULL)

          # Layer checkboxes
          output[[paste0("assoc_layer_ui_", local_proj)]] <- renderUI({
            if (length(series_names) == 0) return(p("No data"))
            checkboxGroupInput(paste0("assoc_layers_", local_proj),
                               label = NULL, choices = series_names,
                               selected = series_names)
          })

          # Select All / None
          observeEvent(input[[paste0("assoc_all_", local_proj)]], {
            updateCheckboxGroupInput(session, paste0("assoc_layers_", local_proj),
                                     selected = series_names)
          })
          observeEvent(input[[paste0("assoc_none_", local_proj)]], {
            updateCheckboxGroupInput(session, paste0("assoc_layers_", local_proj),
                                     selected = character(0))
          })

          # Populate region dropdown based on region type toggle
          observe({
            rt <- input[[paste0("assoc_region_type_", local_proj)]]
            if (is.null(rt)) return()
            regions <- if (rt == "combined") ad$regions_combined else ad$regions_trait
            if (is.null(regions) || nrow(regions) == 0) return()
            regions <- regions[order(min_pvalue)]
            # Region length (already extended by REGION_DISTANCE in pipeline)
            fmt_len <- ifelse(regions$length >= 1e6,
                              paste0(round(regions$length / 1e6, 1), ' Mb'),
                       ifelse(regions$length >= 1e3,
                              paste0(round(regions$length / 1e3, 1), ' kb'),
                              paste0(regions$length, ' bp')))
            ch <- c("(click SNP or select)" = "",
                    setNames(regions$region_id,
                             paste0(regions$region_id, ' (', regions$snp_count,
                                    ' SNPs, ', fmt_len, ')')))
            updateSelectInput(session, paste0("assoc_region_select_", local_proj), choices = ch)
          })

          # Region selection from dropdown
          observeEvent(input[[paste0("assoc_region_select_", local_proj)]], {
            val <- input[[paste0("assoc_region_select_", local_proj)]]
            if (!is.null(val) && val != "") sel_region(val)
          }, ignoreInit = TRUE)

          # Clear region
          observeEvent(input[[paste0("assoc_clear_region_", local_proj)]], {
            sel_region(NULL)
            updateSelectInput(session, paste0("assoc_region_select_", local_proj), selected = "")
          })

          # Click on Manhattan plot -> find region for clicked SNP
          observeEvent(event_data("plotly_click", source = paste0("manh_", local_proj)), {
            click <- event_data("plotly_click", source = paste0("manh_", local_proj))
            if (is.null(click) || is.null(click$customdata)) return()
            snpid <- click$customdata[1]
            if (is.null(mpd)) return()
            snp <- mpd$data[SNPID == snpid][1]
            if (nrow(snp) == 0) return()
            rt <- input[[paste0("assoc_region_type_", local_proj)]]
            regions <- if (rt == "combined") ad$regions_combined else ad$regions_trait
            if (is.null(regions)) return()
            # Match using region boundaries (already extended by REGION_DISTANCE in pipeline)
            matching <- regions[as.character(chr) == as.character(snp$chr) &
                                start <= snp$pos &
                                end >= snp$pos]
            if (nrow(matching) > 0) {
              best <- matching[which.min(min_pvalue), region_id]
              sel_region(best)
              updateSelectInput(session, paste0("assoc_region_select_", local_proj),
                                selected = best)
            }
          })

          # ---- MANHATTAN PLOT ----
          # Strategy: non-sig SNPs as plotly scatter (hoverinfo disabled),
          # sig SNPs as interactive plotly traces.
          # For large datasets (>50k non-sig), downsample non-sig for performance.
          output[[paste0("assoc_manhattan_", local_proj)]] <- renderPlotly({
            if (is.null(mpd)) {
              return(plot_ly() %>%
                add_annotations(text = "No association data", x = 0.5, y = 0.5,
                                xref = "paper", yref = "paper", showarrow = FALSE))
            }

            active <- input[[paste0("assoc_layers_", local_proj)]]
            show_reg <- input[[paste0("assoc_show_regions_", local_proj)]]
            top_n <- input[[paste0("assoc_top_n_", local_proj)]]
            if (is.null(active)) active <- series_names
            if (is.null(top_n)) top_n <- 10

            plot_dt <- mpd$data[series %in% active]
            chr_info <- mpd$chr_info

            if (nrow(plot_dt) == 0) {
              return(plot_ly() %>%
                add_annotations(text = "No series selected", x = 0.5, y = 0.5,
                                xref = "paper", yref = "paper", showarrow = FALSE))
            }

            nonsig <- plot_dt[is_sig == FALSE]
            sig <- plot_dt[is_sig == TRUE]
            y_max <- max(plot_dt$neglog10p) * 1.1

            # Downsample non-sig if very large (>50k points per series)
            max_nonsig <- 50000
            if (nrow(nonsig) > max_nonsig) {
              nonsig <- nonsig[sample(.N, max_nonsig)]
            }

            p <- plot_ly(source = paste0("manh_", local_proj))

            # Non-significant points: color=trait, symbol=method, not interactive
            for (s in unique(plot_dt$series)) {
              sd_ns <- nonsig[series == s]
              if (nrow(sd_ns) > 0) {
                tr <- as.character(sd_ns$trait[1])
                mt <- sd_ns$method[1]
                p <- add_trace(p, data = sd_ns,
                  x = ~pos_cum, y = ~neglog10p,
                  type = "scatter", mode = "markers",
                  marker = list(size = 4, color = t_colors[tr], opacity = 0.4,
                                symbol = m_symbols[mt]),
                  name = s, legendgroup = s, showlegend = TRUE,
                  hoverinfo = "skip")
              }
            }

            # Significant points: color=trait, symbol=method, interactive
            for (s in unique(sig$series)) {
              sd_s <- sig[series == s]
              if (nrow(sd_s) > 0) {
                tr <- as.character(sd_s$trait[1])
                mt <- sd_s$method[1]
                p <- add_trace(p, data = sd_s,
                  x = ~pos_cum, y = ~neglog10p,
                  type = "scatter", mode = "markers",
                  marker = list(size = 10, color = t_colors[tr], opacity = 0.9,
                                symbol = m_symbols[mt],
                                line = list(color = 'black', width = 1.5)),
                  name = paste0(s, ' *'), legendgroup = s, showlegend = FALSE,
                  text = ~paste0('* ', SNPID, '\n', series,
                                 '\np=', formatC(pvalue, format = "e", digits = 2),
                                 '\n-log10(p)=', round(neglog10p, 2)),
                  hoverinfo = "text", customdata = ~SNPID)
              }
            }

            # Bonferroni threshold line
            n_snps <- nrow(plot_dt) / max(1, length(unique(plot_dt$series)))
            bonf_line <- -log10(0.05 / max(1, n_snps))
            p <- add_trace(p,
              x = c(min(plot_dt$pos_cum), max(plot_dt$pos_cum)),
              y = c(bonf_line, bonf_line),
              type = "scatter", mode = "lines",
              line = list(color = "red", dash = "dash", width = 1),
              name = "Bonferroni", showlegend = TRUE, hoverinfo = "skip")

            # Region rectangles as shapes
            shapes <- list()
            if (isTRUE(show_reg)) {
              rt <- input[[paste0("assoc_region_type_", local_proj)]]
              regions <- if (rt == "combined") ad$regions_combined else ad$regions_trait
              if (!is.null(regions) && nrow(regions) > 0) {
                regions <- regions[order(min_pvalue)]
                top_reg <- head(regions, top_n)
                n_reg <- nrow(top_reg)
                reg_pal <- if (n_reg <= 12) {
                  RColorBrewer::brewer.pal(max(3, n_reg), 'Set3')[1:n_reg]
                } else { viridis::turbo(n_reg) }
                for (i in 1:n_reg) {
                  r <- top_reg[i]
                  offset <- chr_info[chr == r$chr, cum_offset]
                  if (length(offset) == 0) next
                  # Region boundaries (already extended by REGION_DISTANCE in pipeline)
                  shapes[[length(shapes) + 1]] <- list(
                    type = "rect",
                    x0 = r$start + offset, x1 = r$end + offset,
                    y0 = 0, y1 = y_max,
                    fillcolor = reg_pal[i], opacity = 0.35,
                    line = list(color = reg_pal[i], width = 1.5))
                }
              }
            }

            # Layout
            layout(p,
              xaxis = list(title = "Chromosome", tickmode = "array",
                           tickvals = chr_info$midpoint,
                           ticktext = as.character(chr_info$chr),
                           showgrid = FALSE),
              yaxis = list(title = "-log10(p-value)", showgrid = TRUE,
                           gridcolor = '#eee'),
              shapes = shapes, hovermode = "closest", dragmode = "zoom",
              legend = list(orientation = "v", x = 1.02, y = 1,
                            font = list(size = 10)),
              margin = list(r = 150),
              plot_bgcolor = 'white', paper_bgcolor = 'white') %>%
              event_register('plotly_click')
          })

          # ---- VALUE BOXES ----
          output[[paste0("assoc_vb_snps_unique_", local_proj)]] <- renderValueBox({
            n <- if (!is.null(ad$selected)) nrow(ad$selected) else 0
            valueBox(n, "Unique Sig. SNPs", icon = icon("dot-circle-o"), color = "red")
          })
          output[[paste0("assoc_vb_snps_hits_", local_proj)]] <- renderValueBox({
            n <- if (!is.null(mpd)) sum(mpd$data$is_sig) else 0
            valueBox(n, "Sig. Hits on Plot", icon = icon("exclamation-circle"), color = "orange")
          })
          output[[paste0("assoc_vb_regions_", local_proj)]] <- renderValueBox({
            n <- if (!is.null(ad$regions_combined)) nrow(ad$regions_combined) else 0
            valueBox(n, "Regions", icon = icon("th"), color = "blue")
          })
          output[[paste0("assoc_vb_genes_", local_proj)]] <- renderValueBox({
            n <- if (!is.null(ad$genes_collapsed)) nrow(ad$genes_collapsed) else 0
            valueBox(n, "Genes", icon = icon("list"), color = "green")
          })
          output[[paste0("assoc_vb_enrich_", local_proj)]] <- renderValueBox({
            base <- paste0('/pipeline/', local_proj,
                           '_results/tables/association/enrichment')
            n <- if (dir.exists(base)) {
              length(list.files(base, pattern = '_enrichment\\.tsv$', recursive = TRUE))
            } else 0
            valueBox(n, "GO Terms", icon = icon("pie-chart"), color = "yellow")
          })
          output[[paste0("assoc_vb_methods_", local_proj)]] <- renderValueBox({
            n <- length(ad$methods)
            valueBox(n, "Methods", icon = icon("flask"), color = "purple")
          })

          # ---- REGION INFO BAR ----
          output[[paste0("assoc_region_info_", local_proj)]] <- renderUI({
            rid <- sel_region()
            if (is.null(rid) || rid == "") return(NULL)
            rt <- input[[paste0("assoc_region_type_", local_proj)]]
            regions <- if (rt == "combined") ad$regions_combined else ad$regions_trait
            if (is.null(regions)) return(NULL)
            r <- regions[region_id == rid]
            if (nrow(r) == 0) return(NULL)
            # Region boundaries (already extended by REGION_DISTANCE in pipeline)
            fmt_ext <- if (r$length[1] >= 1e6) paste0(round(r$length[1]/1e6, 2), ' Mb')
                       else if (r$length[1] >= 1e3) paste0(round(r$length[1]/1e3, 1), ' kb')
                       else paste0(r$length[1], ' bp')
            div(style = paste0("background: #d9edf7; padding: 12px; margin: 5px 15px;",
                               " border-left: 4px solid #31708f; border-radius: 3px;"),
              strong(paste0("Region: ", rid)),
              span(paste0(" | Chr ", r$chr[1], ": ",
                           format(r$start[1], big.mark = ","), " - ",
                           format(r$end[1], big.mark = ","),
                           " (", fmt_ext, ")",
                           " | ", r$snp_count[1], " SNPs",
                           " | min p = ", formatC(r$min_pvalue[1],
                                                   format = "e", digits = 2))))
          })

          # ---- GENES TABLE ----
          output[[paste0("assoc_genes_prompt_", local_proj)]] <- renderUI({
            rid <- sel_region()
            if (is.null(rid) || rid == "") {
              div(style = "text-align: center; padding: 40px; color: #999;",
                icon("hand-pointer-o"),
                h5("Click a significant SNP or select a region to view genes"))
            }
          })

          output[[paste0("assoc_genes_dt_", local_proj)]] <- renderDT({
            rid <- sel_region()
            req(rid, rid != "")
            genes <- ad$genes
            if (is.null(genes) || nrow(genes) == 0) return(NULL)

            # Try exact match on region_id
            rg <- genes[region_id == rid]
            # Fallback: match by coordinates (for combined regions)
            if (nrow(rg) == 0) {
              rt <- input[[paste0("assoc_region_type_", local_proj)]]
              regions_tbl <- if (rt == "combined") ad$regions_combined else ad$regions_trait
              if (!is.null(regions_tbl)) {
                r <- regions_tbl[region_id == rid]
                if (nrow(r) > 0) {
                  rg <- genes[as.character(chr) == as.character(r$chr[1]) &
                              gene_end >= r$start[1] & gene_start <= r$end[1]]
                }
              }
            }
            if (nrow(rg) == 0) return(NULL)
            rg <- unique(rg, by = 'gene_id')

            dcols <- intersect(
              c('gene_id', 'Name', 'description', 'chr', 'gene_start', 'gene_end',
                'exon_snp_count', 'promoter_snp_count', 'biotype', 'ontology'),
              names(rg))
            datatable(rg[, ..dcols], rownames = FALSE,
              options = list(pageLength = 5, scrollX = TRUE,
                             scrollY = "220px", scrollCollapse = TRUE),
              selection = 'single')
          })

          # ---- ENRICHMENT TABLE ----
          output[[paste0("assoc_enrich_prompt_", local_proj)]] <- renderUI({
            rid <- sel_region()
            if (is.null(rid) || rid == "") {
              div(style = "text-align: center; padding: 40px; color: #999;",
                icon("hand-pointer-o"),
                h5("Click a significant SNP or select a region to view enrichment"))
            }
          })

          output[[paste0("assoc_enrich_dt_", local_proj)]] <- renderDT({
            rid <- sel_region()
            req(rid, rid != "")
            enrich <- load_region_enrichment(local_proj, rid)
            if (is.null(enrich) || nrow(enrich) == 0) {
              return(datatable(
                data.frame(Message = "No enrichment results for this region"),
                rownames = FALSE, options = list(dom = 't')))
            }
            dcols <- intersect(
              c('GO_id', 'description', 'gene_ratio', 'pvalue', 'p_adjust',
                'gene_count', 'gene_ids'),
              names(enrich))
            datatable(enrich[, ..dcols], rownames = FALSE,
              options = list(pageLength = 5, scrollX = TRUE,
                             scrollY = "220px", scrollCollapse = TRUE),
              selection = 'single') %>%
              formatSignif(columns = intersect(c('pvalue', 'p_adjust'), dcols),
                           digits = 3)
          })

        } # end if has_assoc_data
        
        ###############################################################
        # MALADAPTATION TAB SERVER LOGIC
        ###############################################################
        gf_suffs <- find_gf_suffixes(local_proj)
        if (length(gf_suffs) > 0) {

          observeEvent(input[[paste0("gf_suffix_", local_proj)]], {
            sfx <- input[[paste0("gf_suffix_", local_proj)]]
            if (is.null(sfx) || sfx == "") return()

            # Importance plots
            output[[paste0("gf_importance_", local_proj)]] <-
              render_plot(local_proj, paste0('OverallImportance_', sfx, '\\.qs$'))
            output[[paste0("gf_importance_msg_", local_proj)]] <-
              render_noplot_message(local_proj, paste0('OverallImportance_', sfx, '\\.qs$'))

            output[[paste0("gf_cumimp_", local_proj)]] <-
              render_plot(local_proj, paste0('CumulativeImportance_', sfx, '\\.qs$'))
            output[[paste0("gf_cumimp_msg_", local_proj)]] <-
              render_noplot_message(local_proj, paste0('CumulativeImportance_', sfx, '\\.qs$'))

            # Genetic Offset PieMaps
            output[[paste0("gf_piemap_", local_proj)]] <-
              render_plot(local_proj, paste0('GeneticOffsetPieMap_', sfx, '\\.qs$'))
            output[[paste0("gf_piemap_msg_", local_proj)]] <-
              render_noplot_message(local_proj, paste0('GeneticOffsetPieMap_', sfx, '\\.qs$'))

            output[[paste0("gf_piemap_tajima_", local_proj)]] <-
              render_plot(local_proj, paste0('GeneticOffsetPieMap_', sfx, '_TajimaD\\.qs$'))
            output[[paste0("gf_piemap_tajima_msg_", local_proj)]] <-
              render_noplot_message(local_proj, paste0('GeneticOffsetPieMap_', sfx, '_TajimaD\\.qs$'))

            output[[paste0("gf_piemap_diversity_", local_proj)]] <-
              render_plot(local_proj, paste0('GeneticOffsetPieMap_', sfx, '_PiDiversity\\.qs$'))
            output[[paste0("gf_piemap_diversity_msg_", local_proj)]] <-
              render_noplot_message(local_proj, paste0('GeneticOffsetPieMap_', sfx, '_PiDiversity\\.qs$'))

            # Zoom maps (dynamic UI)
            zooms <- find_gf_zooms(local_proj, sfx)
            output[[paste0("gf_zoom_ui_", local_proj)]] <- renderUI({
              if (length(zooms) == 0) return(NULL)
              zoom_boxes <- lapply(zooms, function(ztag) {
                # Format zoom label: 35p45_35p55_32p95_33p10 -> 35.45-35.55, 32.95-33.10
                coords <- gsub('p', '.', strsplit(ztag, '_')[[1]])
                label <- if (length(coords) == 4) {
                  paste0('Lon ', coords[1], '-', coords[2],
                         ', Lat ', coords[3], '-', coords[4])
                } else ztag
                box(width = 4, title = paste("Zoom:", label),
                    status = "info", solidHeader = TRUE,
                  div(style = "position: relative; min-height: 250px;",
                    plotOutput(paste0("gf_zoom_", ztag, "_", local_proj)),
                    uiOutput(paste0("gf_zoom_msg_", ztag, "_", local_proj))
                  )
                )
              })
              fluidRow(zoom_boxes)
            })

            # Render each zoom plot
            for (zt in zooms) {
              local({
                local_zt <- zt
                output[[paste0("gf_zoom_", local_zt, "_", local_proj)]] <-
                  render_plot(local_proj, paste0('GeneticOffsetPieMap_', sfx,
                                                  '_zoom_', local_zt, '\\.qs$'))
                output[[paste0("gf_zoom_msg_", local_zt, "_", local_proj)]] <-
                  render_noplot_message(local_proj, paste0('GeneticOffsetPieMap_', sfx,
                                                            '_zoom_', local_zt, '\\.qs$'))
              })
            }

            # Site values table
            site_file <- paste0('/pipeline/', local_proj,
                                '_results/tables/gradientForest/GeneticOffset_site_',
                                sfx, '.tsv')
            output[[paste0("gf_site_dt_", local_proj)]] <- renderDT({
              if (!file.exists(site_file)) return(NULL)
              dt <- data.table::fread(site_file)
              datatable(dt, rownames = FALSE,
                options = list(pageLength = 10, scrollX = TRUE),
                selection = 'single') %>%
                formatRound(columns = names(dt)[sapply(dt, is.numeric)], digits = 4)
            })

            # Site info summary
            output[[paste0("gf_site_info_", local_proj)]] <- renderUI({
              if (!file.exists(site_file)) return(NULL)
              dt <- data.table::fread(site_file)
              offset_col <- grep('offset|genetic_offset', names(dt),
                                 ignore.case = TRUE, value = TRUE)
              if (length(offset_col) > 0) {
                vals <- dt[[offset_col[1]]]
                div(style = "padding-top: 8px;",
                  strong("Genetic Offset Summary: "),
                  span(paste0("mean = ", round(mean(vals, na.rm = TRUE), 4),
                              ", range = [", round(min(vals, na.rm = TRUE), 4),
                              " - ", round(max(vals, na.rm = TRUE), 4), "]",
                              " | ", nrow(dt), " sites")))
              }
            })
          })
        } # end if gf_suffs

        ###############################################################
        # PHENOTYPE ASSOCIATION TAB SERVER LOGIC
        ###############################################################
        pheno_k_val <- find_pheno_assoc_k(local_proj)
        has_pheno_data <- !is.null(pheno_k_val) && length(pheno_k_val) > 0

        # Initialize to NULL (needed by Overlapping Regions tab)
        pad <- NULL; pmpd <- NULL

        if (has_pheno_data) {
          # Load data eagerly (pipeline outputs don't change during session)
          pad <- load_pheno_assoc_data(local_proj, pheno_k_val[1])
          pmpd <- prep_manhattan(pad$pvals, pad$selected)
          # Use trait name as series (method shown via shape only)
          if (!is.null(pmpd)) pmpd$data[, series := as.character(trait)]
          p_series_names <- if (!is.null(pmpd)) sort(unique(pmpd$data$series)) else character(0)
          p_trait_names <- if (!is.null(pmpd)) sort(unique(as.character(pmpd$data$trait))) else character(0)
          p_method_names <- if (!is.null(pmpd)) sort(unique(pmpd$data$method)) else character(0)
          p_t_colors <- if (length(p_trait_names) > 0) generate_trait_colors(p_trait_names) else character(0)
          p_m_symbols <- if (length(p_method_names) > 0) generate_method_symbols(p_method_names) else character(0)

          # Reactive value for selected region
          pheno_sel_region <- reactiveVal(NULL)

          # Trait selector for PieMap/Distribution
          output[[paste0("pheno_trait_selector_ui_", local_proj)]] <- renderUI({
            if (length(p_trait_names) == 0) return(p("No phenotype traits available"))
            selectInput(paste0("pheno_trait_select_", local_proj),
                        "Select phenotype trait:",
                        choices = p_trait_names,
                        selected = p_trait_names[1])
          })

          # Layer checkboxes
          output[[paste0("pheno_layer_ui_", local_proj)]] <- renderUI({
            if (length(p_series_names) == 0) return(p("No data"))
            checkboxGroupInput(paste0("pheno_layers_", local_proj),
                               label = NULL, choices = p_series_names,
                               selected = p_series_names)
          })

          # Select All / None
          observeEvent(input[[paste0("pheno_all_", local_proj)]], {
            updateCheckboxGroupInput(session, paste0("pheno_layers_", local_proj),
                                     selected = p_series_names)
          })
          observeEvent(input[[paste0("pheno_none_", local_proj)]], {
            updateCheckboxGroupInput(session, paste0("pheno_layers_", local_proj),
                                     selected = character(0))
          })

          # Populate region dropdown based on region type toggle
          observe({
            rt <- input[[paste0("pheno_region_type_", local_proj)]]
            if (is.null(rt)) return()
            regions <- if (rt == "combined") pad$regions_combined else pad$regions_trait
            if (is.null(regions) || nrow(regions) == 0) return()
            regions <- regions[order(min_pvalue)]
            # Region length (already extended by REGION_DISTANCE in pipeline)
            fmt_len <- ifelse(regions$length >= 1e6,
                              paste0(round(regions$length / 1e6, 1), ' Mb'),
                       ifelse(regions$length >= 1e3,
                              paste0(round(regions$length / 1e3, 1), ' kb'),
                              paste0(regions$length, ' bp')))
            ch <- c("(click SNP or select)" = "",
                    setNames(regions$region_id,
                             paste0(regions$region_id, ' (', regions$snp_count,
                                    ' SNPs, ', fmt_len, ')')))
            updateSelectInput(session, paste0("pheno_region_select_", local_proj), choices = ch)
          })

          # Region selection from dropdown
          observeEvent(input[[paste0("pheno_region_select_", local_proj)]], {
            val <- input[[paste0("pheno_region_select_", local_proj)]]
            if (!is.null(val) && val != "") pheno_sel_region(val)
          }, ignoreInit = TRUE)

          # Clear region
          observeEvent(input[[paste0("pheno_clear_region_", local_proj)]], {
            pheno_sel_region(NULL)
            updateSelectInput(session, paste0("pheno_region_select_", local_proj), selected = "")
          })

          # Click on Manhattan plot -> find region for clicked SNP
          observeEvent(event_data("plotly_click", source = paste0("pheno_manh_", local_proj)), {
            click <- event_data("plotly_click", source = paste0("pheno_manh_", local_proj))
            if (is.null(click) || is.null(click$customdata)) return()
            snpid <- click$customdata[1]
            if (is.null(pmpd)) return()
            snp <- pmpd$data[SNPID == snpid][1]
            if (nrow(snp) == 0) return()
            rt <- input[[paste0("pheno_region_type_", local_proj)]]
            regions <- if (rt == "combined") pad$regions_combined else pad$regions_trait
            if (is.null(regions)) return()
            # Region boundaries already extended by REGION_DISTANCE in pipeline
            matching <- regions[as.character(chr) == as.character(snp$chr) &
                                start <= snp$pos &
                                end >= snp$pos]
            if (nrow(matching) > 0) {
              best <- matching[which.min(min_pvalue), region_id]
              pheno_sel_region(best)
              updateSelectInput(session, paste0("pheno_region_select_", local_proj),
                                selected = best)
            }
          })

          # ---- PHENOTYPE TRAIT DISTRIBUTION ----
          output[[paste0("pheno_trait_dist_", local_proj)]] <- renderPlotly({
            sel_trait <- input[[paste0("pheno_trait_select_", local_proj)]]
            if (is.null(sel_trait)) {
              return(plot_ly() %>%
                add_annotations(text = "Select a trait", x = 0.5, y = 0.5,
                                xref = "paper", yref = "paper", showarrow = FALSE))
            }
            meta_file <- paste0('/pipeline/', local_proj,
                                '_results/tables/structure/metadata.tsv')
            if (!file.exists(meta_file)) {
              return(plot_ly() %>%
                add_annotations(text = "Metadata file not found", x = 0.5, y = 0.5,
                                xref = "paper", yref = "paper", showarrow = FALSE))
            }
            meta <- data.table::fread(meta_file)
            if (!sel_trait %in% names(meta)) {
              return(plot_ly() %>%
                add_annotations(text = paste0("Trait '", sel_trait, "' not in metadata"),
                                x = 0.5, y = 0.5,
                                xref = "paper", yref = "paper", showarrow = FALSE))
            }

            site_col <- grep('^(site|pop|population)$', names(meta),
                             ignore.case = TRUE, value = TRUE)

            if (length(site_col) > 0) {
              plot_ly(data = meta, x = ~get(sel_trait), color = ~get(site_col[1]),
                      type = "histogram", opacity = 0.7) %>%
                layout(title = paste0("Distribution of ", sel_trait),
                       xaxis = list(title = sel_trait),
                       yaxis = list(title = "Count"),
                       barmode = "overlay",
                       plot_bgcolor = 'white', paper_bgcolor = 'white')
            } else {
              plot_ly(data = meta, x = ~get(sel_trait), type = "histogram",
                      marker = list(color = '#3c8dbc',
                                    line = list(color = 'white', width = 0.5)),
                      nbinsx = 30) %>%
                layout(title = paste0("Distribution of ", sel_trait),
                       xaxis = list(title = sel_trait),
                       yaxis = list(title = "Count"),
                       plot_bgcolor = 'white', paper_bgcolor = 'white')
            }
          })

          # ---- PHENOTYPE SITE MAP (PIE MAP) ----
          output[[paste0("pheno_site_map_", local_proj)]] <- renderPlot({
            sel_trait <- input[[paste0("pheno_trait_select_", local_proj)]]
            if (is.null(sel_trait)) return(NULL)

            plot_obj <- load_plot(local_proj, paste0('PhenoMap_', sel_trait, '\\.qs$'))
            if (is.null(plot_obj)) {
              return(ggplot() + theme_void() +
                annotate("text", x = 0.5, y = 0.5,
                         label = paste0("Pie map not available for ", sel_trait),
                         size = 5, color = "#999"))
            }
            plot_obj
          }, res = 96)

          # ---- PHENOTYPE MANHATTAN PLOT ----
          output[[paste0("pheno_manhattan_", local_proj)]] <- renderPlotly({
            if (is.null(pmpd)) {
              return(plot_ly() %>%
                add_annotations(text = "No phenotype association data", x = 0.5, y = 0.5,
                                xref = "paper", yref = "paper", showarrow = FALSE))
            }

            active <- input[[paste0("pheno_layers_", local_proj)]]
            show_reg <- input[[paste0("pheno_show_regions_", local_proj)]]
            top_n <- input[[paste0("pheno_top_n_", local_proj)]]
            if (is.null(active)) active <- p_series_names
            if (is.null(top_n)) top_n <- 10

            plot_dt <- pmpd$data[series %in% active]
            chr_info <- pmpd$chr_info

            if (nrow(plot_dt) == 0) {
              return(plot_ly() %>%
                add_annotations(text = "No series selected", x = 0.5, y = 0.5,
                                xref = "paper", yref = "paper", showarrow = FALSE))
            }

            nonsig <- plot_dt[is_sig == FALSE]
            sig <- plot_dt[is_sig == TRUE]
            y_max <- max(plot_dt$neglog10p) * 1.1

            # Downsample non-sig if very large
            max_nonsig <- 50000
            if (nrow(nonsig) > max_nonsig) {
              nonsig <- nonsig[sample(.N, max_nonsig)]
            }

            p <- plot_ly(source = paste0("pheno_manh_", local_proj))

            # Non-significant points (one trace per trait)
            for (s in unique(plot_dt$series)) {
              sd_ns <- nonsig[series == s]
              if (nrow(sd_ns) > 0) {
                tr <- as.character(sd_ns$trait[1])
                mt <- sd_ns$method[1]
                p <- add_trace(p, data = sd_ns,
                  x = ~pos_cum, y = ~neglog10p,
                  type = "scatter", mode = "markers",
                  marker = list(size = 4, color = p_t_colors[tr], opacity = 0.4,
                                symbol = p_m_symbols[mt]),
                  name = s, legendgroup = s, showlegend = TRUE,
                  hoverinfo = "skip")
              }
            }

            # Significant points (shape = method)
            for (s in unique(sig$series)) {
              sd_s <- sig[series == s]
              if (nrow(sd_s) > 0) {
                tr <- as.character(sd_s$trait[1])
                mt <- sd_s$method[1]
                p <- add_trace(p, data = sd_s,
                  x = ~pos_cum, y = ~neglog10p,
                  type = "scatter", mode = "markers",
                  marker = list(size = 10, color = p_t_colors[tr], opacity = 0.9,
                                symbol = p_m_symbols[mt],
                                line = list(color = 'black', width = 1.5)),
                  name = paste0(s, ' *'), legendgroup = s, showlegend = FALSE,
                  text = ~paste0('* ', SNPID, '\n', trait, ' [', method, ']',
                                 '\np=', formatC(pvalue, format = "e", digits = 2),
                                 '\n-log10(p)=', round(neglog10p, 2)),
                  hoverinfo = "text", customdata = ~SNPID)
              }
            }

            # Bonferroni threshold line
            n_snps <- nrow(plot_dt) / max(1, length(unique(plot_dt$series)))
            bonf_line <- -log10(0.05 / max(1, n_snps))
            p <- add_trace(p,
              x = c(min(plot_dt$pos_cum), max(plot_dt$pos_cum)),
              y = c(bonf_line, bonf_line),
              type = "scatter", mode = "lines",
              line = list(color = "red", dash = "dash", width = 1),
              name = "Bonferroni", showlegend = TRUE, hoverinfo = "skip")

            # Region rectangles
            shapes <- list()
            if (isTRUE(show_reg)) {
              rt <- input[[paste0("pheno_region_type_", local_proj)]]
              regions <- if (rt == "combined") pad$regions_combined else pad$regions_trait
              if (!is.null(regions) && nrow(regions) > 0) {
                regions <- regions[order(min_pvalue)]
                top_reg <- head(regions, top_n)
                n_reg <- nrow(top_reg)
                reg_pal <- if (n_reg <= 12) {
                  RColorBrewer::brewer.pal(max(3, n_reg), 'Set3')[1:n_reg]
                } else { viridis::turbo(n_reg) }
                for (i in 1:n_reg) {
                  r <- top_reg[i]
                  offset <- chr_info[chr == r$chr, cum_offset]
                  if (length(offset) == 0) next
                  # Region boundaries already extended in pipeline
                  shapes[[length(shapes) + 1]] <- list(
                    type = "rect",
                    x0 = r$start + offset, x1 = r$end + offset,
                    y0 = 0, y1 = y_max,
                    fillcolor = reg_pal[i], opacity = 0.35,
                    line = list(color = reg_pal[i], width = 1.5))
                }
              }
            }

            layout(p,
              xaxis = list(title = "Chromosome", tickmode = "array",
                           tickvals = chr_info$midpoint,
                           ticktext = as.character(chr_info$chr),
                           showgrid = FALSE),
              yaxis = list(title = "-log10(p-value)", showgrid = TRUE,
                           gridcolor = '#eee'),
              shapes = shapes, hovermode = "closest", dragmode = "zoom",
              legend = list(orientation = "v", x = 1.02, y = 1,
                            font = list(size = 10)),
              margin = list(r = 150),
              plot_bgcolor = 'white', paper_bgcolor = 'white') %>%
              event_register('plotly_click')
          })

          # ---- PHENOTYPE VALUE BOXES ----
          output[[paste0("pheno_vb_snps_unique_", local_proj)]] <- renderValueBox({
            n <- if (!is.null(pad$selected)) nrow(pad$selected) else 0
            valueBox(n, "Unique Sig. SNPs", icon = icon("dot-circle-o"), color = "red")
          })
          output[[paste0("pheno_vb_snps_hits_", local_proj)]] <- renderValueBox({
            n <- if (!is.null(pmpd)) sum(pmpd$data$is_sig) else 0
            valueBox(n, "Sig. Hits on Plot", icon = icon("exclamation-circle"), color = "orange")
          })
          output[[paste0("pheno_vb_regions_", local_proj)]] <- renderValueBox({
            n <- if (!is.null(pad$regions_combined)) nrow(pad$regions_combined) else 0
            valueBox(n, "Regions", icon = icon("th"), color = "blue")
          })
          output[[paste0("pheno_vb_genes_", local_proj)]] <- renderValueBox({
            n <- if (!is.null(pad$genes_collapsed)) nrow(pad$genes_collapsed) else 0
            valueBox(n, "Genes", icon = icon("list"), color = "green")
          })
          output[[paste0("pheno_vb_enrich_", local_proj)]] <- renderValueBox({
            base <- paste0('/pipeline/', local_proj,
                           '_results/tables/association_phenotypes/enrichment')
            n <- if (dir.exists(base)) {
              length(list.files(base, pattern = '_enrichment\\.tsv$', recursive = TRUE))
            } else 0
            valueBox(n, "GO Terms", icon = icon("pie-chart"), color = "yellow")
          })
          output[[paste0("pheno_vb_methods_", local_proj)]] <- renderValueBox({
            n <- length(pad$methods)
            valueBox(n, "Methods", icon = icon("flask"), color = "purple")
          })

          # ---- PHENOTYPE REGION INFO BAR ----
          output[[paste0("pheno_region_info_", local_proj)]] <- renderUI({
            rid <- pheno_sel_region()
            if (is.null(rid) || rid == "") return(NULL)
            rt <- input[[paste0("pheno_region_type_", local_proj)]]
            regions <- if (rt == "combined") pad$regions_combined else pad$regions_trait
            if (is.null(regions)) return(NULL)
            r <- regions[region_id == rid]
            if (nrow(r) == 0) return(NULL)
            # Region boundaries (already extended by REGION_DISTANCE in pipeline)
            fmt_ext <- if (r$length[1] >= 1e6) paste0(round(r$length[1]/1e6, 2), ' Mb')
                       else if (r$length[1] >= 1e3) paste0(round(r$length[1]/1e3, 1), ' kb')
                       else paste0(r$length[1], ' bp')
            div(style = paste0("background: #fcf8e3; padding: 12px; margin: 5px 15px;",
                               " border-left: 4px solid #8a6d3b; border-radius: 3px;"),
              strong(paste0("Region: ", rid)),
              span(paste0(" | Chr ", r$chr[1], ": ",
                           format(r$start[1], big.mark = ","), " - ",
                           format(r$end[1], big.mark = ","),
                           " (", fmt_ext, ")",
                           " | ", r$snp_count[1], " SNPs",
                           " | min p = ", formatC(r$min_pvalue[1],
                                                   format = "e", digits = 2))))
          })

          # ---- PHENOTYPE GENES TABLE ----
          output[[paste0("pheno_genes_prompt_", local_proj)]] <- renderUI({
            rid <- pheno_sel_region()
            if (is.null(rid) || rid == "") {
              div(style = "text-align: center; padding: 40px; color: #999;",
                icon("hand-pointer-o"),
                h5("Click a significant SNP or select a region to view genes"))
            }
          })

          output[[paste0("pheno_genes_dt_", local_proj)]] <- renderDT({
            rid <- pheno_sel_region()
            req(rid, rid != "")
            rt <- input[[paste0("pheno_region_type_", local_proj)]]
            genes <- if (!is.null(rt) && rt == "combined") pad$genes_combined else pad$genes
            if (is.null(genes) || nrow(genes) == 0) return(NULL)

            rg <- genes[region_id == rid]
            if (nrow(rg) == 0) {
              regions_tbl <- if (!is.null(rt) && rt == "combined") pad$regions_combined else pad$regions_trait
              if (!is.null(regions_tbl)) {
                r <- regions_tbl[region_id == rid]
                if (nrow(r) > 0) {
                  # Region boundaries already extended in pipeline
                  rg <- genes[as.character(chr) == as.character(r$chr[1]) &
                              gene_end >= r$start[1] &
                              gene_start <= r$end[1]]
                }
              }
            }
            if (nrow(rg) == 0) return(NULL)
            rg <- unique(rg, by = 'gene_id')

            dcols <- intersect(
              c('gene_id', 'Name', 'description', 'chr', 'gene_start', 'gene_end',
                'exon_snp_count', 'promoter_snp_count', 'biotype', 'ontology'),
              names(rg))
            datatable(rg[, ..dcols], rownames = FALSE,
              options = list(pageLength = 5, scrollX = TRUE,
                             scrollY = "220px", scrollCollapse = TRUE),
              selection = 'single')
          })

          # ---- PHENOTYPE ENRICHMENT TABLE ----
          output[[paste0("pheno_enrich_prompt_", local_proj)]] <- renderUI({
            rid <- pheno_sel_region()
            if (is.null(rid) || rid == "") {
              div(style = "text-align: center; padding: 40px; color: #999;",
                icon("hand-pointer-o"),
                h5("Click a significant SNP or select a region to view enrichment"))
            }
          })

          output[[paste0("pheno_enrich_dt_", local_proj)]] <- renderDT({
            rid <- pheno_sel_region()
            req(rid, rid != "")
            enrich <- load_region_enrichment_pheno(local_proj, rid)
            if (is.null(enrich) || nrow(enrich) == 0) {
              return(datatable(
                data.frame(Message = "No enrichment results for this region"),
                rownames = FALSE, options = list(dom = 't')))
            }
            dcols <- intersect(
              c('GO_id', 'description', 'gene_ratio', 'pvalue', 'p_adjust',
                'gene_count', 'gene_ids'),
              names(enrich))
            datatable(enrich[, ..dcols], rownames = FALSE,
              options = list(pageLength = 5, scrollX = TRUE,
                             scrollY = "220px", scrollCollapse = TRUE),
              selection = 'single') %>%
              formatSignif(columns = intersect(c('pvalue', 'p_adjust'), dcols),
                           digits = 3)
          })

        } # end if has_pheno_data

        ###############################################################
        # OVERLAPPING REGIONS TAB SERVER LOGIC
        ###############################################################
        # Load pipeline-precomputed overlap data
        ovd <- NULL
        has_overlap_file <- file.exists(paste0('/pipeline/', local_proj,
                                               '_results/tables/overlapping/Overlap_summary.tsv'))

        if (has_overlap_file && has_assoc_data && has_pheno_data) {
          ovd <- load_overlap_data(local_proj)

          # Parse overlap summary stats from SUMMARY rows
          ov_summary_rows <- if (!is.null(ovd$overlap_summary))
            ovd$overlap_summary[gea_region_id == 'SUMMARY'] else NULL
          ov_pairs <- if (!is.null(ovd$overlap_summary))
            ovd$overlap_summary[gea_region_id != 'SUMMARY'] else NULL

          ov_n_gea <- if (!is.null(ov_summary_rows))
            as.integer(ov_summary_rows[gwas_region_id == 'gea_total', overlap_length]) else 0L
          ov_n_gwas <- if (!is.null(ov_summary_rows))
            as.integer(ov_summary_rows[gwas_region_id == 'gwas_total', overlap_length]) else 0L
          ov_n_pairs <- if (!is.null(ov_summary_rows))
            as.integer(ov_summary_rows[gwas_region_id == 'overlap_pairs', overlap_length]) else 0L
          ov_n_merged <- if (!is.null(ovd$regions_combined)) nrow(ovd$regions_combined) else 0L
          ov_n_genes <- if (!is.null(ovd$genes) && nrow(ovd$genes) > 0 && 'gene_id' %in% names(ovd$genes))
            length(unique(ovd$genes$gene_id)) else 0L

          # Count GO terms across all overlap enrichment
          ov_n_go <- 0L
          ov_enrich_base <- paste0('/pipeline/', local_proj, '_results/tables/overlapping/enrichment')
          if (dir.exists(ov_enrich_base)) {
            ov_enrich_files <- list.files(ov_enrich_base, pattern = '_enrichment\\.tsv$',
                                          recursive = TRUE, full.names = TRUE)
            if (length(ov_enrich_files) > 0) {
              ov_all_go <- data.table::rbindlist(lapply(ov_enrich_files, function(f)
                tryCatch(data.table::fread(f), error = function(e) NULL)), fill = TRUE)
              if (nrow(ov_all_go) > 0 && 'GO_id' %in% names(ov_all_go))
                ov_n_go <- length(unique(ov_all_go$GO_id))
            }
          }

          # Build unified color palette across ALL traits (GEA + GWAS)
          all_ov_traits <- sort(unique(c(
            if (!is.null(mpd)) as.character(unique(mpd$data$trait)) else character(0),
            if (!is.null(pmpd)) as.character(unique(pmpd$data$trait)) else character(0)
          )))
          ov_colors <- if (length(all_ov_traits) > 0) generate_trait_colors(all_ov_traits) else character(0)

          # Build unified method symbols across ALL methods (GEA + GWAS)
          all_ov_methods <- sort(unique(c(
            if (!is.null(mpd)) unique(mpd$data$method) else character(0),
            if (!is.null(pmpd)) unique(pmpd$data$method) else character(0)
          )))
          ov_symbols <- if (length(all_ov_methods) > 0) generate_method_symbols(all_ov_methods) else character(0)

          # Build series names for GEA and GWAS (for layer checkboxes)
          gea_series <- if (!is.null(mpd)) sort(unique(mpd$data$series)) else character(0)
          gwas_series <- if (!is.null(pmpd)) sort(unique(pmpd$data$series)) else character(0)
          ov_series_gea <- paste0('GEA: ', gea_series)
          ov_series_gwas <- paste0('GWAS: ', gwas_series)
          ov_all_series <- c(ov_series_gea, ov_series_gwas)

          # --- Summary value boxes ---
          output[[paste0("ov_box_gea_", local_proj)]] <- renderValueBox({
            valueBox(ov_n_gea, "GEA Regions", icon = icon("cloud-sun"), color = "blue")
          })
          output[[paste0("ov_box_gwas_", local_proj)]] <- renderValueBox({
            valueBox(ov_n_gwas, "GWAS Regions", icon = icon("dna"), color = "purple")
          })
          output[[paste0("ov_box_pairs_", local_proj)]] <- renderValueBox({
            valueBox(ov_n_pairs, "Overlap Pairs", icon = icon("link"), color = "green")
          })
          output[[paste0("ov_box_merged_", local_proj)]] <- renderValueBox({
            valueBox(ov_n_merged, "New Combined", icon = icon("object-group"), color = "yellow")
          })
          output[[paste0("ov_box_genes_", local_proj)]] <- renderValueBox({
            valueBox(ov_n_genes, "Genes", icon = icon("bars"), color = "orange")
          })
          output[[paste0("ov_box_go_", local_proj)]] <- renderValueBox({
            valueBox(ov_n_go, "GO Terms", icon = icon("sitemap"), color = "red")
          })

          # --- Layer checkboxes ---
          output[[paste0("ov_layer_ui_", local_proj)]] <- renderUI({
            if (length(ov_all_series) == 0) return(p("No data"))
            checkboxGroupInput(paste0("ov_layers_", local_proj),
                               label = NULL, choices = ov_all_series,
                               selected = ov_all_series)
          })

          # Select All / None
          observeEvent(input[[paste0("ov_all_", local_proj)]], {
            updateCheckboxGroupInput(session, paste0("ov_layers_", local_proj),
                                     selected = ov_all_series)
          })
          observeEvent(input[[paste0("ov_none_", local_proj)]], {
            updateCheckboxGroupInput(session, paste0("ov_layers_", local_proj),
                                     selected = character(0))
          })

          # --- Region selection ---
          ov_sel_region <- reactiveVal(NULL)

          # Populate region dropdown based on region type toggle
          observe({
            rt <- input[[paste0("ov_region_type_", local_proj)]]
            if (is.null(rt)) return()
            regions <- if (rt == "combined") ovd$regions_combined else ovd$regions_trait
            if (is.null(regions) || nrow(regions) == 0) return()
            regions <- regions[order(min_pvalue)]
            # Region length (already extended by REGION_DISTANCE in pipeline)
            fmt_len <- ifelse(regions$length >= 1e6,
                              paste0(round(regions$length / 1e6, 1), ' Mb'),
                       ifelse(regions$length >= 1e3,
                              paste0(round(regions$length / 1e3, 1), ' kb'),
                              paste0(regions$length, ' bp')))
            src_label <- if ('source' %in% names(regions)) paste0(' [', regions$source, ']') else ''
            ch <- c("(click SNP or select)" = "",
                    setNames(regions$region_id,
                             paste0(regions$region_id, ' (', regions$snp_count,
                                    ' SNPs, ', fmt_len, src_label, ')')))
            updateSelectInput(session, paste0("ov_region_select_", local_proj), choices = ch)
          })

          # Region selection from dropdown
          observeEvent(input[[paste0("ov_region_select_", local_proj)]], {
            val <- input[[paste0("ov_region_select_", local_proj)]]
            if (!is.null(val) && val != "") ov_sel_region(val)
          }, ignoreInit = TRUE)

          # Clear region
          observeEvent(input[[paste0("ov_clear_region_", local_proj)]], {
            ov_sel_region(NULL)
            updateSelectInput(session, paste0("ov_region_select_", local_proj), selected = "")
          })

          # Click on Miami plot -> find region for clicked SNP
          observeEvent(event_data("plotly_click", source = paste0("ov_miami_", local_proj)), {
            click <- event_data("plotly_click", source = paste0("ov_miami_", local_proj))
            if (is.null(click) || is.null(click$customdata)) return()
            snpid <- click$customdata[1]
            # Search in both GEA and GWAS data
            snp <- NULL
            if (!is.null(mpd)) snp <- mpd$data[SNPID == snpid][1]
            if ((is.null(snp) || nrow(snp) == 0) && !is.null(pmpd))
              snp <- pmpd$data[SNPID == snpid][1]
            if (is.null(snp) || nrow(snp) == 0) return()
            rt <- input[[paste0("ov_region_type_", local_proj)]]
            regions <- if (rt == "combined") ovd$regions_combined else ovd$regions_trait
            if (is.null(regions)) return()
            # Region boundaries already extended by REGION_DISTANCE in pipeline
            matching <- regions[as.character(chr) == as.character(snp$chr) &
                                start <= snp$pos &
                                end >= snp$pos]
            if (nrow(matching) > 0) {
              best <- matching[which.min(min_pvalue), region_id]
              ov_sel_region(best)
              updateSelectInput(session, paste0("ov_region_select_", local_proj),
                                selected = best)
            }
          })

          # --- MIAMI MANHATTAN PLOT ---
          output[[paste0("ov_miami_", local_proj)]] <- renderPlotly({
            req(mpd, pmpd)

            active <- input[[paste0("ov_layers_", local_proj)]]
            show_reg <- input[[paste0("ov_show_regions_", local_proj)]]
            top_n <- input[[paste0("ov_top_n_", local_proj)]]
            if (is.null(active)) active <- ov_all_series
            if (is.null(top_n)) top_n <- 10

            # Filter active GEA and GWAS series
            active_gea <- sub('^GEA: ', '', active[grepl('^GEA: ', active)])
            active_gwas <- sub('^GWAS: ', '', active[grepl('^GWAS: ', active)])

            # GEA data (positive y)
            gea_d <- copy(mpd$data)
            gea_d <- gea_d[series %in% active_gea]
            gea_d[, y_val := neglog10p]

            # GWAS data (negative y)
            gwas_d <- copy(pmpd$data)
            gwas_d <- gwas_d[series %in% active_gwas]
            gwas_d[, y_val := -neglog10p]

            # Shared chr_info for consistent x-axis
            chr_info <- mpd$chr_info

            # Recalculate GWAS cumulative positions using GEA's chromosome layout
            gwas_d <- merge(gwas_d[, !c('cum_offset', 'pos_cum'), with = FALSE],
                            chr_info[, .(chr, cum_offset)], by = 'chr', all.x = TRUE)
            gwas_d[, pos_cum := pos + cum_offset]
            gwas_d <- gwas_d[!is.na(pos_cum)]

            if (nrow(gea_d) == 0 && nrow(gwas_d) == 0) {
              return(plot_ly() %>%
                add_annotations(text = "No series selected", x = 0.5, y = 0.5,
                                xref = "paper", yref = "paper", showarrow = FALSE))
            }

            # Determine y-axis range
            y_max <- max(c(
              if (nrow(gea_d) > 0) max(gea_d$neglog10p) else 1,
              if (nrow(gwas_d) > 0) max(gwas_d$neglog10p) else 1
            ), na.rm = TRUE) * 1.1

            p <- plot_ly(source = paste0("ov_miami_", local_proj))

            # --- GEA non-sig (above) ---
            gea_nonsig <- gea_d[is_sig == FALSE]
            if (nrow(gea_nonsig) > 50000) {
              set.seed(42); gea_nonsig <- gea_nonsig[sample(.N, 50000)]
            }
            for (s in unique(gea_d$series)) {
              sd_ns <- gea_nonsig[series == s]
              if (nrow(sd_ns) > 0) {
                tr <- as.character(sd_ns$trait[1])
                mt <- sd_ns$method[1]
                col <- if (tr %in% names(ov_colors)) ov_colors[[tr]] else '#999999'
                sym <- if (mt %in% names(ov_symbols)) ov_symbols[[mt]] else 'circle'
                p <- add_trace(p, data = sd_ns, x = ~pos_cum, y = ~y_val,
                  type = 'scatter', mode = 'markers',
                  marker = list(size = 4, color = col, opacity = 0.4, symbol = sym),
                  hoverinfo = 'skip', showlegend = TRUE,
                  name = paste0('GEA: ', s), legendgroup = paste0('gea_', s))
              }
            }

            # --- GWAS non-sig (below) ---
            gwas_nonsig <- gwas_d[is_sig == FALSE]
            if (nrow(gwas_nonsig) > 50000) {
              set.seed(42); gwas_nonsig <- gwas_nonsig[sample(.N, 50000)]
            }
            for (s in unique(gwas_d$series)) {
              sd_ns <- gwas_nonsig[series == s]
              if (nrow(sd_ns) > 0) {
                tr <- as.character(sd_ns$trait[1])
                mt <- sd_ns$method[1]
                col <- if (tr %in% names(ov_colors)) ov_colors[[tr]] else '#999999'
                sym <- if (mt %in% names(ov_symbols)) ov_symbols[[mt]] else 'circle'
                p <- add_trace(p, data = sd_ns, x = ~pos_cum, y = ~y_val,
                  type = 'scatter', mode = 'markers',
                  marker = list(size = 4, color = col, opacity = 0.4, symbol = sym),
                  hoverinfo = 'skip', showlegend = TRUE,
                  name = paste0('GWAS: ', s), legendgroup = paste0('gwas_', s))
              }
            }

            # --- GEA sig (above) ---
            gea_sig <- gea_d[is_sig == TRUE]
            for (s in unique(gea_sig$series)) {
              sd_s <- gea_sig[series == s]
              if (nrow(sd_s) > 0) {
                tr <- as.character(sd_s$trait[1])
                mt <- sd_s$method[1]
                col <- if (tr %in% names(ov_colors)) ov_colors[[tr]] else '#333333'
                sym <- if (mt %in% names(ov_symbols)) ov_symbols[[mt]] else 'circle'
                p <- add_trace(p, data = sd_s, x = ~pos_cum, y = ~y_val,
                  type = 'scatter', mode = 'markers',
                  marker = list(size = 10, color = col, opacity = 0.9, symbol = sym,
                                line = list(color = 'black', width = 1.5)),
                  text = ~paste0('* ', SNPID, '\n', series,
                                 '\np=', formatC(pvalue, format = "e", digits = 2),
                                 '\n-log10(p)=', round(neglog10p, 2)),
                  hoverinfo = 'text', customdata = ~SNPID,
                  name = paste0('GEA: ', s, ' *'), legendgroup = paste0('gea_', s),
                  showlegend = FALSE)
              }
            }

            # --- GWAS sig (below) ---
            gwas_sig <- gwas_d[is_sig == TRUE]
            for (s in unique(gwas_sig$series)) {
              sd_s <- gwas_sig[series == s]
              if (nrow(sd_s) > 0) {
                tr <- as.character(sd_s$trait[1])
                mt <- sd_s$method[1]
                col <- if (tr %in% names(ov_colors)) ov_colors[[tr]] else '#333333'
                sym <- if (mt %in% names(ov_symbols)) ov_symbols[[mt]] else 'circle'
                p <- add_trace(p, data = sd_s, x = ~pos_cum, y = ~y_val,
                  type = 'scatter', mode = 'markers',
                  marker = list(size = 10, color = col, opacity = 0.9, symbol = sym,
                                line = list(color = 'black', width = 1.5)),
                  text = ~paste0('* ', SNPID, '\n', series,
                                 '\np=', formatC(pvalue, format = "e", digits = 2),
                                 '\n-log10(p)=', round(neglog10p, 2)),
                  hoverinfo = 'text', customdata = ~SNPID,
                  name = paste0('GWAS: ', s, ' *'), legendgroup = paste0('gwas_', s),
                  showlegend = FALSE)
              }
            }

            # Region rectangles as shapes (from pipeline tables, extended by REGION_DISTANCE)
            shapes <- list()
            if (isTRUE(show_reg)) {
              rt <- input[[paste0("ov_region_type_", local_proj)]]
              regions <- if (rt == "combined") ovd$regions_combined else ovd$regions_trait
              if (!is.null(regions) && nrow(regions) > 0) {
                regions <- regions[order(min_pvalue)]
                top_reg <- head(regions, top_n)
                n_reg <- nrow(top_reg)
                reg_pal <- if (n_reg <= 12) {
                  RColorBrewer::brewer.pal(max(3, n_reg), 'Set3')[1:n_reg]
                } else { viridis::turbo(n_reg) }
                for (i in 1:n_reg) {
                  r <- top_reg[i]
                  offset <- chr_info[chr == r$chr, cum_offset]
                  if (length(offset) == 0) next
                  # Region boundaries already extended in pipeline
                  shapes[[length(shapes) + 1]] <- list(
                    type = "rect",
                    x0 = r$start + offset, x1 = r$end + offset,
                    y0 = -y_max, y1 = y_max,
                    fillcolor = reg_pal[i], opacity = 0.35,
                    line = list(color = reg_pal[i], width = 1.5))
                }
              }
            }

            # Layout with zero line
            layout(p,
              xaxis = list(title = "Chromosome", tickmode = "array",
                           tickvals = chr_info$midpoint,
                           ticktext = as.character(chr_info$chr),
                           showgrid = FALSE),
              yaxis = list(title = "-log10(p-value)",
                           zeroline = TRUE, zerolinecolor = '#333', zerolinewidth = 2,
                           showgrid = TRUE, gridcolor = '#eee'),
              shapes = shapes, hovermode = "closest", dragmode = "zoom",
              legend = list(orientation = "v", x = 1.02, y = 1,
                            font = list(size = 10)),
              margin = list(r = 150),
              plot_bgcolor = 'white', paper_bgcolor = 'white') %>%
              event_register('plotly_click')
          })

          # ---- REGION INFO BAR ----
          output[[paste0("ov_region_info_", local_proj)]] <- renderUI({
            rid <- ov_sel_region()
            if (is.null(rid) || rid == "") return(NULL)
            rt <- input[[paste0("ov_region_type_", local_proj)]]
            regions <- if (rt == "combined") ovd$regions_combined else ovd$regions_trait
            if (is.null(regions)) return(NULL)
            r <- regions[region_id == rid]
            if (nrow(r) == 0) return(NULL)
            # Region boundaries (already extended by REGION_DISTANCE in pipeline)
            fmt_ext <- if (r$length[1] >= 1e6) paste0(round(r$length[1]/1e6, 2), ' Mb')
                       else if (r$length[1] >= 1e3) paste0(round(r$length[1]/1e3, 1), ' kb')
                       else paste0(r$length[1], ' bp')
            src <- if ('source' %in% names(r)) as.character(r$source[1]) else ''
            traits <- if ('traits' %in% names(r)) as.character(r$traits[1]) else
                      if ('trait' %in% names(r)) as.character(r$trait[1]) else ''
            div(style = paste0("background: #d9edf7; padding: 12px; margin: 5px 15px;",
                               " border-left: 4px solid #31708f; border-radius: 3px;"),
              strong(paste0("Region: ", rid)),
              span(paste0(" | Chr ", r$chr[1], ": ",
                           format(r$start[1], big.mark = ","), " - ",
                           format(r$end[1], big.mark = ","),
                           " (", fmt_ext, ")",
                           " | Source: ", src,
                           " | Traits: ", traits,
                           " | ", r$snp_count[1], " SNPs",
                           " | min p = ", formatC(r$min_pvalue[1],
                                                   format = "e", digits = 2))))
          })

          # ---- GENES TABLE ----
          output[[paste0("ov_genes_prompt_", local_proj)]] <- renderUI({
            rid <- ov_sel_region()
            if (is.null(rid) || rid == "") {
              div(style = "text-align: center; padding: 40px; color: #999;",
                icon("hand-pointer-o"),
                h5("Click a significant SNP or select a region to view genes"))
            }
          })

          output[[paste0("ov_genes_dt_", local_proj)]] <- renderDT({
            rid <- ov_sel_region()
            req(rid, rid != "")
            rt <- input[[paste0("ov_region_type_", local_proj)]]
            genes <- if (!is.null(rt) && rt == "combined") ovd$genes_combined else ovd$genes
            if (is.null(genes) || nrow(genes) == 0) return(NULL)
            rg <- genes[region_id == rid]
            # Fallback: coordinate-based matching
            if (nrow(rg) == 0) {
              regions_tbl <- if (!is.null(rt) && rt == "combined") ovd$regions_combined else ovd$regions_trait
              if (!is.null(regions_tbl)) {
                r <- regions_tbl[region_id == rid]
                if (nrow(r) > 0) {
                  rg <- genes[as.character(chr) == as.character(r$chr[1]) &
                              gene_end >= r$start[1] & gene_start <= r$end[1]]
                }
              }
            }
            if (nrow(rg) == 0) return(NULL)
            rg <- unique(rg, by = 'gene_id')
            dcols <- intersect(
              c('gene_id', 'Name', 'description', 'chr', 'gene_start', 'gene_end',
                'exon_snp_count', 'promoter_snp_count', 'biotype', 'ontology'),
              names(rg))
            datatable(rg[, ..dcols], rownames = FALSE,
              options = list(pageLength = 5, scrollX = TRUE,
                             scrollY = "220px", scrollCollapse = TRUE),
              selection = 'single')
          })

          # ---- ENRICHMENT TABLE ----
          output[[paste0("ov_enrich_prompt_", local_proj)]] <- renderUI({
            rid <- ov_sel_region()
            if (is.null(rid) || rid == "") {
              div(style = "text-align: center; padding: 40px; color: #999;",
                icon("hand-pointer-o"),
                h5("Click a significant SNP or select a region to view enrichment"))
            }
          })

          output[[paste0("ov_enrich_dt_", local_proj)]] <- renderDT({
            rid <- ov_sel_region()
            req(rid, rid != "")
            enrich <- load_region_enrichment_overlap(local_proj, rid)
            if (is.null(enrich) || nrow(enrich) == 0) {
              return(datatable(
                data.frame(Message = "No enrichment results for this region"),
                rownames = FALSE, options = list(dom = 't')))
            }
            dcols <- intersect(
              c('GO_id', 'description', 'gene_ratio', 'pvalue', 'p_adjust',
                'gene_count', 'gene_ids'),
              names(enrich))
            datatable(enrich[, ..dcols], rownames = FALSE,
              options = list(pageLength = 5, scrollX = TRUE,
                             scrollY = "220px", scrollCollapse = TRUE),
              selection = 'single') %>%
              formatSignif(columns = intersect(c('pvalue', 'p_adjust'), dcols),
                           digits = 3)
          })

          # ---- OVERLAP PAIRS TABLE ----
          output[[paste0("ov_overlap_dt_", local_proj)]] <- renderDT({
            if (is.null(ov_pairs) || nrow(ov_pairs) == 0) {
              return(datatable(
                data.frame(Message = "No overlapping regions found between GEA and GWAS"),
                rownames = FALSE, options = list(dom = 't')))
            }
            dcols <- intersect(
              c('gea_region_id', 'gwas_region_id', 'chr', 'overlap_length', 'overlap_pct',
                'gea_snp_count', 'gwas_snp_count', 'shared_snps',
                'gea_traits', 'gwas_traits', 'gea_min_pvalue', 'gwas_min_pvalue'),
              names(ov_pairs))
            datatable(ov_pairs[, ..dcols], rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE),
              selection = 'single') %>%
              formatSignif(columns = intersect(c('gea_min_pvalue', 'gwas_min_pvalue'), dcols), digits = 3)
          })

        } # end if has_overlap_file

      })
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)
