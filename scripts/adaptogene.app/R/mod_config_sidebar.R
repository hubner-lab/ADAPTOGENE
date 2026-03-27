#' Config sidebar module
#'
#' A hidden left sidebar that displays pipeline configuration parameters for
#' the active tab. Parameters are split into mandatory (always visible) and
#' optional (in a collapsed accordion). Values are read from the project YAML.
#'
#' Phase 1: read-only display with interactive inputs (no save yet).
#' Phase 2: adds Save button and YAML write-back.
#' Phase 3: adds Run button and progress bar.
#'
#' @param id module namespace id
#' @param tab_title display title shown in the sidebar header
#' @noRd
mod_config_sidebar_ui <- function(id, tab_title = "Config") {
    ns <- shiny::NS(id)

    bslib::sidebar(
        id       = ns("sidebar"),
        title    = htmltools::span(
            class = "d-flex align-items-center gap-2 config-sidebar-title",
            bsicons::bs_icon("sliders2", size = "1em"),
            tab_title
        ),
        open     = "closed",
        width    = "70vw",
        position = "left",
        # Dynamic content rendered by server
        shiny::uiOutput(ns("sidebar_content"))
    )
}

#' @param config_reactive reactive returning project_data list (with $config)
#' @param tab_name tab identifier matching config_schema()$tab values
#' @noRd
mod_config_sidebar_server <- function(id, config_reactive, tab_name) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        output$sidebar_content <- shiny::renderUI({
            pd  <- config_reactive()
            cfg <- if (!is.null(pd)) pd$config else list()

            entries <- schema_for_tab(tab_name)

            if (length(entries) == 0)
                return(htmltools::p(class = "text-muted small p-2",
                                    "No configuration parameters for this tab."))

            split   <- schema_split(entries)
            mand    <- split$mandatory
            opt     <- split$optional

            htmltools::tagList(
                # ── Mandatory (Core Parameters) ──────────────────────────────
                if (length(mand) > 0)
                    htmltools::div(
                        class = "config-section",
                        htmltools::p(
                            class = "config-section-header",
                            bsicons::bs_icon("asterisk", size = "0.75em"),
                            " Core Parameters"
                        ),
                        render_section_groups(ns, mand, cfg)
                    ),

                # ── Optional (Advanced accordion) ─────────────────────────────
                if (length(opt) > 0)
                    bslib::accordion(
                        id   = ns("adv"),
                        open = FALSE,
                        bslib::accordion_panel(
                            title = htmltools::tagList(
                                bsicons::bs_icon("gear", size = "0.9em"),
                                " Advanced"
                            ),
                            value = "advanced",
                            render_section_groups(ns, opt, cfg)
                        )
                    )
            )
        })
    })
}

# ── Internal helpers ───────────────────────────────────────────────────────────

#' Render inputs grouped by their section label
#' @noRd
render_section_groups <- function(ns, entries, config) {
    sections <- unique(sapply(entries, `[[`, "section"))
    multi_section <- length(sections) > 1

    htmltools::tagList(lapply(sections, function(sec) {
        sec_entries <- Filter(function(e) e$section == sec, entries)
        htmltools::tagList(
            if (multi_section)
                htmltools::p(class = "config-subsection-label", sec),
            lapply(sec_entries, function(e) render_config_field(ns, e, config))
        )
    }))
}

#' Render a single config field (label + input + help)
#' @noRd
render_config_field <- function(ns, entry, config) {
    # Read current value from nested config using dot-path key
    value <- config_get_by_path(config, entry$key)

    input_id <- ns(gsub("\\.", "_", entry$key))

    # Help text element
    help_el <- if (!is.null(entry$help))
        htmltools::tags$small(class = "form-text text-muted config-help", entry$help)
    else NULL

    # Build the appropriate input widget
    input_el <- build_config_input(input_id, entry, value)

    # Checkbox already embeds its label; all others get a separate label
    label_el <- if (entry$type != "checkbox")
        htmltools::tags$label(
            class = "form-label config-field-label",
            `for` = input_id,
            entry$label
        )
    else NULL

    htmltools::div(
        class = "config-field",
        label_el,
        input_el,
        help_el
    )
}

#' Build the correct shiny input widget for a schema entry
#' @noRd
build_config_input <- function(input_id, entry, value) {
    type <- entry$type

    # Normalize value for display
    display_val <- normalize_display_value(value, type)

    switch(type,
        "numeric" = shiny::numericInput(
            input_id,
            label = NULL,
            value = display_val,
            min   = entry$min %||% NA,
            max   = entry$max %||% NA,
            step  = entry$step %||% 1,
            width = "100%"
        ),

        "text" = shiny::textInput(
            input_id,
            label       = NULL,
            value       = display_val,
            placeholder = entry$placeholder %||% "",
            width       = "100%"
        ),

        "select" = shiny::selectInput(
            input_id,
            label    = NULL,
            choices  = entry$choices %||% character(0),
            selected = display_val,
            width    = "100%"
        ),

        "checkbox" = shiny::checkboxInput(
            input_id,
            label = entry$label,
            value = isTRUE(display_val)
        ),

        "textarea" = shiny::textAreaInput(
            input_id,
            label       = NULL,
            value       = display_val,
            rows        = 2,
            width       = "100%",
            placeholder = entry$placeholder %||% ""
        ),

        "method_table" = render_method_display(display_val),

        # Fallback
        shiny::textInput(input_id, label = NULL, value = as.character(display_val %||% ""),
                          width = "100%")
    )
}

#' Convert raw config value to a display-ready scalar
#' @noRd
normalize_display_value <- function(value, type) {
    if (is.null(value)) return(switch(type,
        "numeric"  = NA_real_,
        "checkbox" = FALSE,
        "method_table" = list(),
        ""
    ))

    switch(type,
        "checkbox" = {
            # Config may store TRUE/FALSE (logical) or "T"/"F" strings
            if (is.logical(value)) value
            else value %in% c("T", "TRUE", "true", "1", "yes")
        },
        "textarea" = {
            # Lists (e.g. epsilon_range) → comma-separated string
            if (is.list(value) || length(value) > 1)
                paste(value, collapse = ", ")
            else as.character(value)
        },
        "method_table" = value,  # pass list as-is
        # All others: coerce to character or numeric as appropriate
        "numeric" = {
            n <- suppressWarnings(as.numeric(value))
            if (is.na(n)) NA_real_ else n
        },
        {
            # text / select / default
            if (is.logical(value)) as.character(value)
            else if (length(value) > 1) paste(value, collapse = ", ")
            else as.character(value %||% "")
        }
    )
}

#' Read-only display for association.configs method list
#' @noRd
render_method_display <- function(configs_list) {
    if (is.null(configs_list) || length(configs_list) == 0) {
        return(htmltools::div(
            class = "method-display-empty",
            bsicons::bs_icon("exclamation-triangle", size = "0.9em"),
            " No methods configured"
        ))
    }

    # Ensure lowercase keys (config migration may have uppercase)
    rows <- lapply(configs_list, function(m) {
        method    <- m$method    %||% m$METHOD    %||% "?"
        adjust    <- m$adjust    %||% m$ADJUST    %||% "none"
        threshold <- m$threshold %||% m$THRESHOLD %||% "?"
        htmltools::div(
            class = "method-row d-flex align-items-center gap-2 mb-1",
            htmltools::span(class = "badge bg-primary method-badge", method),
            htmltools::span(
                class = "text-muted small",
                paste0(adjust, " \u2264 ", threshold)
            )
        )
    })

    htmltools::div(
        class = "method-display",
        rows,
        htmltools::tags$small(
            class = "text-muted fst-italic d-block mt-1",
            "(editing methods available in next release)"
        )
    )
}
