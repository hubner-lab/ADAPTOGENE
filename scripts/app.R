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
# Function to scan directory and find projects

# Function to find available K values from filenames
find_k_values <- function(project, pipeline_path = '/pipeline/') {
  # Get all files in the intermediate directory
  files <- list.files(paste0(pipeline_path, project, '_results/intermediate'))

  # Extract K values from different plot types
  k_values <- unique(as.numeric(stringr::str_extract(
    grep("(PCA_structurePie_K|SNMF_PopDiffTest_K|SNMF_structure_K)[0-9]+\\.qs",
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
  bio_values <- bio_values %>% sort

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

# Function to load PCA plot
load_plot <- function(project, plotname, pipeline_path = '/pipeline/') {
  tryCatch({
	
 #   if (plotname_pattern == T) {
     plot_path = list.files(paste0(pipeline_path, project, '_results/intermediate'),
						     pattern = plotname,
						     full.names = T)
#    } else {
#    plot_path <- file.path(pipeline_path, 
#                          paste0(project, "_results"),
 #                         "intermediate",
#                          plotname)
#    }
   								message(plot_path)
    if (file.exists(plot_path)) {
      pca_plot <- qread(plot_path)
      return(pca_plot)
    } else {
      stop('NO FILE')
      return(NULL)
    }
  }, error = function(e) {
    warning(paste("Error loading", plotname, "plot for project", project, ":", e$message))
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
                plotOutput(paste0("association_plot_", proj))
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
	output[[paste0("PCA_plot_", local_proj)]] <- render_plot(local_proj, 'PCA.qs')
        output[[paste0("PCA_message_", local_proj)]] <- render_noplot_message(local_proj, 'PCA.qs')

        # pcaTW
        output[[paste0("pcaTW_plot_", local_proj)]] <- render_plot(local_proj, 'PCA_TracyWidow.qs')
        output[[paste0("pcaTW_message_", local_proj)]] <- render_noplot_message(local_proj, 'PCA_TracyWidow.qs')
   	
	#   	SNMF_Cross_entropy_K	
	output[[paste0("CrossEntropy_plot_", local_proj)]] <- render_plot(local_proj, 'SNMF_Cross_entropy_K[0-9]*-[0-9]*.qs')
        output[[paste0("CrossEntropy_message_", local_proj)]] <- render_noplot_message(local_proj, 'SNMF_Cross_entropy_K[0-9]*-[0-9]*.qs')
	
	 # New reactive plot outputs for K-dependent plots
        observeEvent(input[[paste0("k_slider_", local_proj)]], {
          k_value <- input[[paste0("k_slider_", local_proj)]]

          # Structure Pie plot
          output[[paste0("structurePie_plot_", local_proj)]] <-
            render_plot(local_proj, paste0('PCA_structurePie_K', k_value, '.qs'))
          output[[paste0("structurePie_message_", local_proj)]] <-
            render_noplot_message(local_proj, paste0('PCA_structurePie_K', k_value, '.qs'))

          # Population Difference Test plot
          output[[paste0("popDiffTest_plot_", local_proj)]] <-
            render_plot(local_proj, paste0('SNMF_PopDiffTest_K', k_value, '.qs'))
          output[[paste0("popDiffTest_message_", local_proj)]] <-
            render_noplot_message(local_proj, paste0('SNMF_PopDiffTest_K', k_value, '.qs'))

          # Structure plot
          output[[paste0("structure_plot_", local_proj)]] <-
            render_plot(local_proj, paste0('SNMF_structure_K', k_value, '.qs'))
          output[[paste0("structure_message_", local_proj)]] <-
            render_noplot_message(local_proj, paste0('SNMF_structure_K', k_value, '.qs'))
        })

################################################################################################
	# STRUCTURE K

        # CorHM
	output[[paste0("CorHM_plot_", local_proj)]] <- render_plot(local_proj, 'CorrelationHeatmap.qs')
        output[[paste0("CorHM_message_", local_proj)]] <- render_noplot_message(local_proj, 'CorrelationHeatmap.qs')

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

        
        output[[paste0("association_plot_", local_proj)]] <- renderPlot({
          plot.new()
          text(0.5, 0.5, "Association plot placeholder")
        })
        
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
