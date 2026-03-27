#' Config sidebar module
#'
#' A hidden left sidebar that displays and edits pipeline configuration
#' parameters for the active tab. Parameters are split into mandatory
#' (always visible) and optional (in a collapsed accordion).
#'
#' Reads from config_state$working, writes back on input change.
#' Shows a "modified" badge when any field differs from config_state$saved.
#' Provides a Reset button to restore saved values.
#'
#' Phase 2+3: bidirectional sync + run button (runner module is separate).
#'
#' @param id module namespace id
#' @param tab_title display title shown in the sidebar header
#' @noRd
mod_config_sidebar_ui <- function(id, tab_title = "Config", runner_ui = NULL) {
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
        # Dynamic content rendered by server (fires once on project switch)
        shiny::uiOutput(ns("sidebar_content")),
        # Reset button — shown/hidden via shinyjs
        shinyjs::hidden(
            htmltools::div(
                id    = ns("reset_footer"),
                class = "mt-3 pt-2 border-top",
                shiny::actionButton(
                    ns("reset_btn"),
                    label = htmltools::tagList(
                        bsicons::bs_icon("arrow-counterclockwise", size = "0.85em"),
                        " Reset to saved"
                    ),
                    class = "btn btn-outline-secondary btn-sm w-100 config-reset-btn"
                )
            )
        ),
        # Optional runner UI (Run button + progress) injected from app_ui
        if (!is.null(runner_ui))
            htmltools::div(class = "mt-3 pt-2 border-top", runner_ui)
    )
}

#' @param id module namespace id
#' @param config_state reactiveValues with $saved, $working, $project
#' @param tab_name tab identifier matching config_schema()$tab values
#' @return reactive integer: number of modified fields on this tab
#' @noRd
mod_config_sidebar_server <- function(id, config_state, tab_name) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        entries     <- schema_for_tab(tab_name)
        input_ids   <- setNames(
            gsub("\\.", "_", sapply(entries, `[[`, "key")),
            sapply(entries, `[[`, "key")
        )

        # ── Render UI once per project switch ──────────────────────────────────
        # Uses config_state$working at render time. Does NOT re-render when
        # working changes (that would cause circular update loops).
        # Updates are applied via updateXxxInput() instead.

        last_project <- shiny::reactiveVal(NULL)

        output$sidebar_content <- shiny::renderUI({
            # Depend on project identity, not working config values
            proj <- config_state$project
            last_project(proj)

            cfg <- if (!is.null(proj)) config_state$working else list()

            if (length(entries) == 0)
                return(htmltools::p(class = "text-muted small p-2",
                                    "No configuration parameters for this tab."))

            sp    <- schema_split(entries)
            mand  <- sp$mandatory
            opt   <- sp$optional

            htmltools::tagList(
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

        # ── Input observers: write changed values into config_state$working ───
        # Guard: skip if the new value equals the current working value
        # (prevents circular updates from Reset/updateXxxInput calls).

        for (entry in entries) {
            local({
                e        <- entry
                iid      <- gsub("\\.", "_", e$key)
                key_path <- e$key

                shiny::observeEvent(input[[iid]], {
                    if (e$type == "method_table") return()

                    raw     <- input[[iid]]
                    new_val <- input_to_config_value(raw, e$type)
                    cur_val <- config_get_by_path(config_state$working, key_path)

                    # No-op if equal (prevents circular re-triggers)
                    if (config_values_equal(new_val, cur_val)) return()

                    config_state$working <- config_set_by_path(
                        config_state$working, key_path, new_val
                    )
                }, ignoreInit = TRUE)
            })
        }

        # ── Reset button: restore working from saved for this tab's fields ────
        shiny::observeEvent(input$reset_btn, {
            new_working <- config_state$working
            for (e in entries) {
                if (e$type == "method_table") next
                saved_val <- config_get_by_path(config_state$saved, e$key)
                new_working <- config_set_by_path(new_working, e$key, saved_val)
            }
            config_state$working <- new_working

            # Push saved values back to inputs
            update_sidebar_inputs(session, entries, config_state$saved)
        })

        # ── Dirty count: compare working vs saved for this tab's fields ───────
        dirty_count <- shiny::reactive({
            working <- config_state$working
            saved   <- config_state$saved
            sum(vapply(entries, function(e) {
                if (e$type == "method_table") return(FALSE)
                w <- config_get_by_path(working, e$key)
                s <- config_get_by_path(saved,   e$key)
                !config_values_equal(w, s)
            }, logical(1)))
        })

        # ── Badge + Reset visibility via JS/shinyjs ────────────────────────────
        shiny::observe({
            n <- dirty_count()
            # Send custom message to config-dirty.js
            session$sendCustomMessage("config_dirty_update", list(
                sidebar_id  = ns("sidebar"),
                dirty_count = n
            ))
            # Show/hide reset footer
            if (n > 0) shinyjs::show("reset_footer")
            else       shinyjs::hide("reset_footer")
        })

        # ── Return dirty count for use by pipeline runner ─────────────────────
        return(dirty_count)
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
    value    <- config_get_by_path(config, entry$key)
    input_id <- ns(gsub("\\.", "_", entry$key))
    help_el  <- if (!is.null(entry$help))
        htmltools::tags$small(class = "form-text text-muted config-help", entry$help)
    else NULL
    input_el <- build_config_input(input_id, entry, value)
    label_el <- if (entry$type != "checkbox")
        htmltools::tags$label(
            class = "form-label config-field-label",
            `for` = input_id,
            entry$label
        )
    else NULL

    htmltools::div(class = "config-field", label_el, input_el, help_el)
}

#' Build the correct Shiny input widget for a schema entry
#' @noRd
build_config_input <- function(input_id, entry, value) {
    type        <- entry$type
    display_val <- normalize_display_value(value, type)

    switch(type,
        "numeric" = shiny::numericInput(
            input_id, label = NULL, value = display_val,
            min = entry$min %||% NA, max = entry$max %||% NA,
            step = entry$step %||% 1, width = "100%"
        ),
        "text" = shiny::textInput(
            input_id, label = NULL, value = display_val,
            placeholder = entry$placeholder %||% "", width = "100%"
        ),
        "select" = shiny::selectInput(
            input_id, label = NULL,
            choices = entry$choices %||% character(0),
            selected = display_val, width = "100%"
        ),
        "checkbox" = shiny::checkboxInput(
            input_id, label = entry$label, value = isTRUE(display_val)
        ),
        "textarea" = shiny::textAreaInput(
            input_id, label = NULL, value = display_val,
            rows = 2, width = "100%",
            placeholder = entry$placeholder %||% ""
        ),
        "method_table" = render_method_display(display_val),
        shiny::textInput(input_id, label = NULL,
                          value = as.character(display_val %||% ""), width = "100%")
    )
}

#' Push saved values back to already-rendered inputs (called on Reset)
#' @noRd
update_sidebar_inputs <- function(session, entries, saved_config) {
    for (e in entries) {
        if (e$type == "method_table") next
        iid <- gsub("\\.", "_", e$key)
        val <- config_get_by_path(saved_config, e$key)
        dv  <- normalize_display_value(val, e$type)

        switch(e$type,
            "numeric"  = shiny::updateNumericInput(session, iid, value = dv),
            "text"     = shiny::updateTextInput(session, iid, value = dv %||% ""),
            "select"   = shiny::updateSelectInput(session, iid, selected = dv %||% ""),
            "checkbox" = shiny::updateCheckboxInput(session, iid, value = isTRUE(dv)),
            "textarea" = shiny::updateTextAreaInput(session, iid, value = dv %||% "")
        )
    }
}

#' Convert raw config value to a display-ready scalar
#' @noRd
normalize_display_value <- function(value, type) {
    if (is.null(value)) return(switch(type,
        "numeric"      = NA_real_,
        "checkbox"     = FALSE,
        "method_table" = list(),
        ""
    ))
    switch(type,
        "checkbox" = {
            if (is.logical(value)) value
            else value %in% c("T", "TRUE", "true", "1", "yes")
        },
        "textarea" = {
            if (is.list(value) || length(value) > 1)
                paste(value, collapse = ", ")
            else as.character(value)
        },
        "method_table" = value,
        "numeric" = {
            n <- suppressWarnings(as.numeric(value))
            if (is.na(n)) NA_real_ else n
        },
        {
            if (is.logical(value)) as.character(value)
            else if (length(value) > 1) paste(value, collapse = ", ")
            else as.character(value %||% "")
        }
    )
}

#' Compare two config values for equality, tolerating integer/double coercion
#'
#' Sidebar inputs always return doubles for numeric fields, but YAML may
#' store them as integers (e.g. 5L). Using identical() would always report
#' these as dirty. This helper normalizes numeric types before comparing.
#' @noRd
config_values_equal <- function(a, b) {
    if (identical(a, b)) return(TRUE)
    # Numeric type coercion: integer vs double
    if (is.numeric(a) && is.numeric(b)) {
        return(isTRUE(all.equal(as.double(a), as.double(b),
                                tolerance = .Machine$double.eps ^ 0.5)))
    }
    # NULL vs empty string (both mean "not set")
    if ((is.null(a) || identical(a, "")) && (is.null(b) || identical(b, ""))) return(TRUE)
    FALSE
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
    rows <- lapply(configs_list, function(m) {
        method    <- m$method    %||% m$METHOD    %||% "?"
        adjust    <- m$adjust    %||% m$ADJUST    %||% "none"
        threshold <- m$threshold %||% m$THRESHOLD %||% "?"
        htmltools::div(
            class = "method-row d-flex align-items-center gap-2 mb-1",
            htmltools::span(class = "badge bg-primary method-badge", method),
            htmltools::span(class = "text-muted small",
                            paste0(adjust, " \u2264 ", threshold))
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
