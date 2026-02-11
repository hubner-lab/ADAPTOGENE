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
			          })
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

      # Read top N regions default from project config
      assoc_top_n_default <- 10L
      if (has_assoc) {
        cfg_val <- read_project_config(proj, 'ASSOC_TOP_REGIONS')
        if (!is.null(cfg_val)) {
          assoc_top_n_default <- suppressWarnings(as.integer(cfg_val))
          if (is.na(assoc_top_n_default)) assoc_top_n_default <- 10L
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
                    valueBoxOutput(paste0("assoc_vb_snps_", proj), width = 3),
                    valueBoxOutput(paste0("assoc_vb_regions_", proj), width = 3),
                    valueBoxOutput(paste0("assoc_vb_genes_", proj), width = 3),
                    valueBoxOutput(paste0("assoc_vb_enrich_", proj), width = 3)
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
                                      "Top N regions:", min = 1, max = 30,
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
              tabPanel(
                "Maladaptation",
                h3("Maladaptation Analysis"),
                plotOutput(paste0("maladaptation_plot_", proj))
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

        if (has_assoc_data) {
          # Load data eagerly (pipeline outputs don't change during session)
          ad <- load_assoc_data(local_proj, assoc_k_val[1])
          mpd <- prep_manhattan(ad$pvals, ad$selected)
          series_names <- if (!is.null(mpd)) sort(unique(mpd$data$series)) else character(0)
          trait_names <- if (!is.null(mpd)) sort(unique(as.character(mpd$data$trait))) else character(0)
          method_names <- if (!is.null(mpd)) sort(unique(mpd$data$method)) else character(0)
          t_colors <- if (length(trait_names) > 0) generate_trait_colors(trait_names) else character(0)
          m_symbols <- if (length(method_names) > 0) generate_method_symbols(method_names) else character(0)

          # Region extension distance from config (used for rectangles and display)
          region_dist <- 2000000L  # default
          cfg_rd <- read_project_config(local_proj, 'ASSOC_REGION_DISTANCE')
          if (!is.null(cfg_rd)) {
            region_dist <- suppressWarnings(as.integer(cfg_rd))
            if (is.na(region_dist)) region_dist <- 2000000L
          }

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
            # Extended region length: SNP boundaries + REGION_DISTANCE on each side
            ext_len <- regions$length + 2 * region_dist
            fmt_len <- ifelse(ext_len >= 1e6,
                              paste0(round(ext_len / 1e6, 1), ' Mb'),
                       ifelse(ext_len >= 1e3,
                              paste0(round(ext_len / 1e3, 1), ' kb'),
                              paste0(ext_len, ' bp')))
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
            # Match using extended region boundaries
            matching <- regions[as.character(chr) == as.character(snp$chr) &
                                (start - region_dist) <= snp$pos &
                                (end + region_dist) >= snp$pos]
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
                  # Extended region: SNP boundaries + REGION_DISTANCE on each side
                  ext_start <- max(0, r$start - region_dist)
                  ext_end <- r$end + region_dist
                  shapes[[length(shapes) + 1]] <- list(
                    type = "rect",
                    x0 = ext_start + offset, x1 = ext_end + offset,
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
          output[[paste0("assoc_vb_snps_", local_proj)]] <- renderValueBox({
            n <- if (!is.null(ad$selected)) nrow(ad$selected) else 0
            valueBox(n, "Sig. SNPs", icon = icon("exclamation-circle"), color = "red")
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

          # ---- REGION INFO BAR ----
          output[[paste0("assoc_region_info_", local_proj)]] <- renderUI({
            rid <- sel_region()
            if (is.null(rid) || rid == "") return(NULL)
            rt <- input[[paste0("assoc_region_type_", local_proj)]]
            regions <- if (rt == "combined") ad$regions_combined else ad$regions_trait
            if (is.null(regions)) return(NULL)
            r <- regions[region_id == rid]
            if (nrow(r) == 0) return(NULL)
            # Extended region boundaries (SNP positions + REGION_DISTANCE)
            ext_start <- max(0, r$start[1] - region_dist)
            ext_end <- r$end[1] + region_dist
            ext_len <- ext_end - ext_start
            fmt_ext <- if (ext_len >= 1e6) paste0(round(ext_len/1e6, 2), ' Mb')
                       else if (ext_len >= 1e3) paste0(round(ext_len/1e3, 1), ' kb')
                       else paste0(ext_len, ' bp')
            div(style = paste0("background: #d9edf7; padding: 12px; margin: 5px 15px;",
                               " border-left: 4px solid #31708f; border-radius: 3px;"),
              strong(paste0("Region: ", rid)),
              span(paste0(" | Chr ", r$chr[1], ": ",
                           format(ext_start, big.mark = ","), " - ",
                           format(ext_end, big.mark = ","),
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
        
        output[[paste0("maladaptation_plot_", local_proj)]] <- renderPlot({
          plot.new()
          text(0.5, 0.5, "Maladaptation plot placeholder")
        })
      })
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)
