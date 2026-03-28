#' Piemap viewer module UI
#'
#' Displays a piemap image card (bio variable, with optional metric/zoom
#' variants). Selectors are provided externally (e.g., in the parent sidebar).
#'
#' @param id module namespace id
#' @noRd
mod_piemap_viewer_ui <- function(id) {
    htmltools::div(
        class = "piemap-container",
        mod_image_card_ui(shiny::NS(id)("piemap"))
    )
}

#' Piemap viewer module server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param bio reactive integer: bio variable number (1 for bio1, etc.)
#' @param metric reactive character: "none", "tajima_d", or "pi_diversity"
#' @param zoom reactive character: zoom region tag, or NULL/""/"none" for global
#' @noRd
mod_piemap_viewer_server <- function(id, project_data,
                                      bio    = shiny::reactive(1L),
                                      metric = shiny::reactive("none"),
                                      zoom   = shiny::reactive(NULL)) {
    shiny::moduleServer(id, function(input, output, session) {

        effective_metric <- shiny::reactive({
            m <- metric()
            if (is.null(m) || !nzchar(m) || m == "none") NULL else m
        })

        effective_zoom <- shiny::reactive({
            z <- zoom()
            if (is.null(z) || !nzchar(z) || z == "none") NULL else z
        })

        raw_path <- shiny::reactive({
            pd <- project_data()
            piemap_path(pd$name, bio(), effective_metric(), effective_zoom())
        })

        # Suppress the image when the no-variance flag is present (even if file exists)
        path <- shiny::reactive({
            p    <- raw_path()
            flag <- base_flag()
            if (file_ok(p) && file.exists(flag)) NULL else p
        })

        title <- shiny::reactive({
            b <- bio()
            m <- effective_metric()
            z <- effective_zoom()
            base <- paste0("Piemap bio", b)
            if (!is.null(z))      paste0(base, " (zoom: ", z, ")")
            else if (!is.null(m)) paste0(base, " (", gsub("_", " ", m), ")")
            else                  base
        })

        dl_name <- shiny::reactive({
            b <- bio()
            m <- effective_metric()
            z <- effective_zoom()
            base <- paste0("piemap_bio", b)
            if (!is.null(z))      paste0(base, "_zoom_", z)
            else if (!is.null(m)) paste0(base, "_", m)
            else                  base
        })

        # Base piemap path (no metric/zoom) used for the no-variance flag lookup.
        # The flag is written once per bio variable by the pipeline.
        base_flag <- shiny::reactive({
            pd <- project_data()
            base_png <- piemap_path(pd$name, bio(), NULL, NULL)
            sub("\\.png$", "_no_spatial_variance.flag", base_png)
        })

        # Determine placeholder message + suggestion based on file state
        placeholder <- shiny::reactive({
            p <- path()
            if (file_ok(p)) return("Plot not available")
            if (file.exists(base_flag())) "No spatial variance" else "Piemap not available"
        })

        suggestion <- shiny::reactive({
            p <- path()
            if (file_ok(p)) return(NULL)
            flag <- base_flag()
            if (file.exists(flag)) {
                lines <- tryCatch(readLines(flag), error = function(e) character(0))
                range_val <- gsub("^range=", "", lines[grepl("^range=", lines)])
                mean_val  <- gsub("^mean=",  "", lines[grepl("^mean=",  lines)])
                if (length(range_val) && length(mean_val)) {
                    paste0("This climate variable has near-zero spatial variation ",
                           "across sampling sites (range=", range_val,
                           ", mean=", mean_val, "). ",
                           "It is unlikely to drive population differentiation.")
                } else {
                    "This climate variable has near-zero spatial variation across sampling sites."
                }
            } else {
                "Run mode=structure_K to generate piemaps for all bioclimatic variables."
            }
        })

        mod_image_card_server("piemap",
            path        = path,
            title       = title,
            dl_name     = dl_name,
            placeholder = placeholder,
            suggestion  = suggestion
        )
    })
}
