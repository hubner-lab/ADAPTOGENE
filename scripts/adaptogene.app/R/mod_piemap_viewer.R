#' Piemap viewer module UI
#'
#' Displays a piemap image card (bio variable, with optional metric/zoom
#' variants). Selectors are provided externally (e.g., in the parent sidebar).
#'
#' @param id module namespace id
#' @noRd
mod_piemap_viewer_ui <- function(id) {
    mod_image_card_ui(shiny::NS(id)("piemap"))
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

        path <- shiny::reactive({
            pd <- project_data()
            piemap_path(pd$name, bio(), effective_metric(), effective_zoom())
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

        mod_image_card_server("piemap",
            path    = path,
            title   = title,
            dl_name = dl_name
        )
    })
}
