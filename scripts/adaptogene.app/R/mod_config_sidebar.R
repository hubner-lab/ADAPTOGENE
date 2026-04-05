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

            cfg <- if (!is.null(proj)) shiny::isolate(config_state$working) else list()

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
                        render_section_groups(ns, mand, cfg, proj)
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
                            render_section_groups(ns, opt, cfg, proj)
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

                if (e$type == "method_table") {
                    # method_table: observe a hidden JSON bridge input
                    json_id <- paste0(iid, "_json")
                    shiny::observeEvent(input[[json_id]], {
                        raw <- input[[json_id]]
                        if (is.null(raw) || nchar(raw) == 0) return()
                        tryCatch({
                            new_val <- jsonlite::fromJSON(raw, simplifyVector = FALSE)
                            cur_val <- config_get_by_path(config_state$working, key_path)
                            if (identical(new_val, cur_val)) return()
                            config_state$working <- config_set_by_path(
                                config_state$working, key_path, new_val
                            )
                        }, error = function(e) NULL)
                    }, ignoreInit = TRUE, ignoreNULL = TRUE)
                } else {
                    shiny::observeEvent(input[[iid]], {
                        raw     <- input[[iid]]
                        new_val <- input_to_config_value(raw, e$type)
                        cur_val <- config_get_by_path(config_state$working, key_path)

                        # No-op if equal (prevents circular re-triggers)
                        if (config_values_equal(new_val, cur_val)) return()

                        config_state$working <- config_set_by_path(
                            config_state$working, key_path, new_val
                        )
                    }, ignoreInit = TRUE)
                }
            })
        }

        # ── Reset button: restore working from saved for this tab's fields ────
        shiny::observeEvent(input$reset_btn, {
            new_working <- config_state$working
            for (e in entries) {
                saved_val <- config_get_by_path(config_state$saved, e$key)
                new_working <- config_set_by_path(new_working, e$key, saved_val)
            }
            config_state$working <- new_working

            # Push saved values back to inputs (including method tables via JS)
            update_sidebar_inputs(session, entries, config_state$saved)
        })

        # ── Dirty count: compare working vs saved for this tab's fields ───────
        dirty_count <- shiny::reactive({
            working <- config_state$working
            saved   <- config_state$saved
            sum(vapply(entries, function(e) {
                w <- config_get_by_path(working, e$key)
                s <- config_get_by_path(saved,   e$key)
                if (e$type == "method_table") {
                    !identical(w, s)
                } else {
                    !config_values_equal(w, s)
                }
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
render_section_groups <- function(ns, entries, config, project = NULL) {
    sections <- unique(sapply(entries, `[[`, "section"))
    multi_section <- length(sections) > 1

    htmltools::tagList(lapply(sections, function(sec) {
        sec_entries <- Filter(function(e) e$section == sec, entries)
        htmltools::tagList(
            if (multi_section)
                htmltools::p(class = "config-subsection-label", sec),
            lapply(sec_entries, function(e) render_config_field(ns, e, config, project))
        )
    }))
}

#' Render a single config field (label + input + help)
#' @noRd
render_config_field <- function(ns, entry, config, project = NULL) {
    value    <- config_get_by_path(config, entry$key)
    input_id <- ns(gsub("\\.", "_", entry$key))
    help_el  <- if (!is.null(entry$help))
        htmltools::tags$small(class = "form-text text-muted config-help", entry$help)
    else NULL
    input_el <- build_config_input(input_id, entry, value, project)
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
build_config_input <- function(input_id, entry, value, project = NULL) {
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
        "method_table" = {
            # input_id is already ns()-wrapped; extract the local id for textInput
            # by stripping the namespace prefix (everything up to and incl. last "-")
            local_id <- sub("^.*-", "", input_id)
            render_method_editor(input_id, local_id, display_val)
        },
        "bio_chips" = {
            inv <- character(0)
            if (!is.null(project)) {
                inv_path <- climate_invariant_path(project)
                if (file.exists(inv_path)) {
                    inv <- tryCatch(
                        data.table::fread(inv_path, sep = "\t", header = TRUE)$predictor,
                        error = function(e) character(0)
                    )
                }
            }
            render_bio_chips(input_id, display_val, invariant = inv)
        },
        shiny::textInput(input_id, label = NULL,
                          value = as.character(display_val %||% ""), width = "100%")
    )
}

#' Push saved values back to already-rendered inputs (called on Reset)
#' @noRd
update_sidebar_inputs <- function(session, entries, saved_config) {
    for (e in entries) {
        iid <- gsub("\\.", "_", e$key)
        val <- config_get_by_path(saved_config, e$key)

        if (e$type == "method_table") {
            # Reset method editor via JS: overwrite with saved JSON
            # input_id must match the full_json_id used in the editor's JS handler
            json_val <- tryCatch(as.character(jsonlite::toJSON(val %||% list(), auto_unbox = FALSE)),
                                 error = function(e) "[]")
            session$sendCustomMessage("method_editor_reset", list(
                input_id = paste0(session$ns(iid), "_json"),
                configs  = json_val
            ))
            next
        }

        if (e$type == "bio_chips") {
            session$sendCustomMessage("bio_chips_reset", list(
                container_id = paste0(session$ns(iid), "_container"),
                input_id     = session$ns(iid),
                value        = as.character(val %||% "")
            ))
            next
        }

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
        "bio_chips"    = "",
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

#' Inline editable rows for association.configs method list
#'
#' Renders a compact editor: one row per method (method select, adjust select,
#' threshold numeric, remove button), an "Add method" button, and a hidden
#' textInput JSON bridge that Shiny observers watch.
#' @noRd
render_method_editor <- function(input_id, local_id, configs_list) {
    if (is.null(configs_list)) configs_list <- list()

    METHOD_CHOICES <- c("EMMAX", "LFMM", "GLM", "MLM", "CMLM",
                        "ECMLM", "SUPER", "MLMM", "FarmCPU", "BLINK")
    ADJUST_CHOICES <- c("bonf", "qval", "top")

    # Normalise each entry
    configs_norm <- lapply(configs_list, function(m) {
        list(
            method    = m$method    %||% m$METHOD    %||% "EMMAX",
            adjust    = m$adjust    %||% m$ADJUST    %||% "bonf",
            threshold = as.character(m$threshold %||% m$THRESHOLD %||% "0.05")
        )
    })

    # Initial JSON for the bridge (used by JS on page load)
    init_json <- tryCatch(jsonlite::toJSON(configs_norm, auto_unbox = TRUE),
                          error = function(e) "[]")

    # Build one HTML row per config entry (pure HTML, no Shiny widgets)
    make_row <- function(idx, m) {
        row_id <- paste0(input_id, "_row_", idx)

        method_opts <- paste(sapply(METHOD_CHOICES, function(ch) {
            sel <- if (ch == m$method) " selected" else ""
            sprintf('<option value="%s"%s>%s</option>', ch, sel, ch)
        }), collapse = "")

        adjust_opts <- paste(sapply(ADJUST_CHOICES, function(ch) {
            sel <- if (ch == m$adjust) " selected" else ""
            sprintf('<option value="%s"%s>%s</option>', ch, sel, ch)
        }), collapse = "")

        htmltools::tags$div(
            class = "method-row",
            id    = row_id,
            `data-idx` = idx,
            htmltools::tags$select(
                class    = "form-select form-select-sm method-select",
                `data-role` = "method",
                htmltools::HTML(method_opts)
            ),
            htmltools::tags$select(
                class    = "form-select form-select-sm adjust-select",
                `data-role` = "adjust",
                htmltools::HTML(adjust_opts)
            ),
            htmltools::tags$input(
                type      = "number",
                class     = "form-control form-control-sm threshold-input",
                `data-role` = "threshold",
                value     = m$threshold,
                min       = "0", max = "1", step = "0.001"
            ),
            htmltools::tags$button(
                class    = "method-remove-btn",
                type     = "button",
                `data-role` = "remove",
                bsicons::bs_icon("trash3", size = "0.85em")
            )
        )
    }

    rows <- if (length(configs_norm) > 0)
        mapply(make_row, seq_along(configs_norm), configs_norm, SIMPLIFY = FALSE)
    else
        list()

    # local_json_id: non-namespaced, used for shiny::textInput (Shiny adds ns internally)
    # full_json_id: namespaced DOM id used by JS to locate the input element
    local_json_id <- paste0(local_id, "_json")
    full_json_id  <- paste0(input_id, "_json")   # ns already applied to input_id

    htmltools::tagList(
        # Hidden JSON bridge — use full namespaced id so Shiny module routing works
        shiny::textInput(full_json_id, label = NULL, value = init_json) |>
            shinyjs::hidden(),

        # Visible editor container
        htmltools::div(
            class           = "method-editor",
            id              = paste0(input_id, "_editor"),
            `data-json-id`  = full_json_id,
            rows,
            htmltools::tags$button(
                class    = "btn btn-outline-secondary btn-sm method-add-btn",
                type     = "button",
                `data-role` = "add",
                bsicons::bs_icon("plus-circle", size = "0.85em"),
                " Add method"
            )
        ),

        # JS: wire up delegated event handlers for this editor
        htmltools::tags$script(htmltools::HTML(sprintf(
'(function() {
  var editorId  = "%s";
  var jsonInputId = "%s";

  function collectRows(editor) {
    var rows = [];
    editor.querySelectorAll(".method-row").forEach(function(row) {
      rows.push({
        method:    row.querySelector("[data-role=method]").value,
        adjust:    row.querySelector("[data-role=adjust]").value,
        threshold: row.querySelector("[data-role=threshold]").value
      });
    });
    return rows;
  }

  function pushToShiny(editor) {
    var json = JSON.stringify(collectRows(editor));
    var el = document.getElementById(jsonInputId);
    if (el) {
      el.value = json;
      el.dispatchEvent(new Event("change", {bubbles: true}));
      Shiny.setInputValue(jsonInputId, json, {priority: "event"});
    }
  }

  var METHOD_CHOICES = %s;
  var ADJUST_CHOICES = %s;

  function makeRow(m) {
    var div = document.createElement("div");
    div.className = "method-row";
    div.setAttribute("data-role-container", "row");
    ["method", "adjust"].forEach(function(role) {
      var sel = document.createElement("select");
      sel.className = "form-select form-select-sm " + role + "-select";
      sel.setAttribute("data-role", role);
      var choices = role === "method" ? METHOD_CHOICES : ADJUST_CHOICES;
      var def = role === "method" ? (m ? m.method : "EMMAX") : (m ? m.adjust : "bonf");
      choices.forEach(function(ch) {
        var opt = document.createElement("option");
        opt.value = ch; opt.text = ch;
        if (ch === def) opt.selected = true;
        sel.appendChild(opt);
      });
      div.appendChild(sel);
    });
    var num = document.createElement("input");
    num.type = "number"; num.className = "form-control form-control-sm threshold-input";
    num.setAttribute("data-role", "threshold");
    num.value = m ? m.threshold : "0.05";
    num.min = "0"; num.max = "1"; num.step = "0.001";
    div.appendChild(num);
    var btn = document.createElement("button");
    btn.type = "button"; btn.className = "method-remove-btn";
    btn.setAttribute("data-role", "remove");
    btn.innerHTML = "<svg xmlns=\\"http://www.w3.org/2000/svg\\" width=\\"12\\" height=\\"12\\" fill=\\"currentColor\\" viewBox=\\"0 0 16 16\\"><path d=\\"M5.5 5.5A.5.5 0 0 1 6 6v6a.5.5 0 0 1-1 0V6a.5.5 0 0 1 .5-.5zm2.5 0a.5.5 0 0 1 .5.5v6a.5.5 0 0 1-1 0V6a.5.5 0 0 1 .5-.5zm3 .5a.5.5 0 0 0-1 0v6a.5.5 0 0 0 1 0V6z\\"/><path fill-rule=\\"evenodd\\" d=\\"M14.5 3a1 1 0 0 1-1 1H13v9a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V4h-.5a1 1 0 0 1-1-1V2a1 1 0 0 1 1-1H6a1 1 0 0 1 1-1h2a1 1 0 0 1 1 1h3.5a1 1 0 0 1 1 1v1zM4.118 4 4 4.059V13a1 1 0 0 0 1 1h6a1 1 0 0 0 1-1V4.059L11.882 4H4.118zM2.5 3V2h11v1h-11z\\"/></svg>";
    div.appendChild(btn);
    return div;
  }

  document.addEventListener("click", function(e) {
    var editor = document.getElementById(editorId);
    if (!editor) return;
    var role = e.target.closest("[data-role]");
    if (!role) return;
    if (role.getAttribute("data-role") === "remove") {
      role.closest(".method-row").remove();
      pushToShiny(editor);
    } else if (role.getAttribute("data-role") === "add") {
      var addBtn = editor.querySelector("[data-role=add]");
      editor.insertBefore(makeRow(null), addBtn);
      pushToShiny(editor);
    }
  });

  document.addEventListener("change", function(e) {
    var editor = document.getElementById(editorId);
    if (!editor || !editor.contains(e.target)) return;
    if (e.target.getAttribute("data-role") && e.target.getAttribute("data-role") !== "add" && e.target.getAttribute("data-role") !== "remove") {
      pushToShiny(editor);
    }
  });

  // Reset handler: overwrite editor contents from new JSON
  Shiny.addCustomMessageHandler("method_editor_reset", function(msg) {
    if (msg.input_id !== "%s") return;
    var editor = document.getElementById(editorId);
    if (!editor) return;
    editor.querySelectorAll(".method-row").forEach(function(r) { r.remove(); });
    var configs = JSON.parse(msg.configs);
    var addBtn = editor.querySelector("[data-role=add]");
    configs.forEach(function(m) { editor.insertBefore(makeRow(m), addBtn); });
    // Update hidden bridge silently (no Shiny event — avoid round-trip)
    var el = document.getElementById(jsonInputId);
    if (el) el.value = msg.configs;
  });
})();
',
            paste0(input_id, "_editor"),  # %s 1: editorId
            full_json_id,                 # %s 2: jsonInputId (DOM id, namespaced)
            jsonlite::toJSON(METHOD_CHOICES, auto_unbox = FALSE),  # %s 3
            jsonlite::toJSON(ADJUST_CHOICES, auto_unbox = FALSE),  # %s 4
            full_json_id                  # %s 5: method_editor_reset id check
        )))
    )
}

#' Render a grid of bio_1..bio_19 toggle chips
#'
#' Replaces the plain text input for climate.predictors with clickable pills.
#' A hidden textInput (input_id) carries the comma-separated value for Shiny.
#' @noRd
render_bio_chips <- function(input_id, display_val, invariant = character(0)) {
    all_bios <- paste0("bio_", 1:19)

    # Parse current value; default to all selected when empty/NULL
    selected <- if (!is.null(display_val) && nzchar(trimws(display_val))) {
        trimws(strsplit(as.character(display_val), ",")[[1]])
    } else {
        all_bios
    }
    # Force-remove invariant predictors from selected (they cannot be active)
    selected <- setdiff(selected, invariant)

    chips <- lapply(all_bios, function(b) {
        is_inv <- b %in% invariant
        cls <- if (is_inv) {
            "bc-chip bc-invariant"
        } else if (b %in% selected) {
            "bc-chip bc-active"
        } else {
            "bc-chip"
        }
        htmltools::tags$button(
            type               = "button",
            class              = cls,
            `data-bio`         = b,
            `data-invariant`   = if (is_inv) "true" else NULL,
            title              = if (is_inv) paste0(b, ": no variance in current map extent") else NULL,
            b
        )
    })

    container_id <- paste0(input_id, "_container")

    js_code <- sprintf('
(function() {
    var container = document.getElementById("%s");
    var inputEl   = document.getElementById("%s");
    if (!container || !inputEl) return;

    function syncValue() {
        var active = container.querySelectorAll(".bc-chip.bc-active");
        var vals = Array.from(active).map(function(el) { return el.dataset.bio; });
        var csv = vals.join(",");
        inputEl.value = csv;
        Shiny.setInputValue("%s", csv, {priority: "event"});
    }

    container.addEventListener("click", function(e) {
        var chip = e.target.closest(".bc-chip");
        if (!chip) return;
        if (chip.dataset.invariant === "true") return;
        chip.classList.toggle("bc-active");
        syncValue();
    });
})();
', container_id, input_id, input_id)

    # Wrap chips + hidden bridge
    htmltools::tagList(
        htmltools::div(
            id    = container_id,
            class = "bio-chips"
        ) |> htmltools::tagAppendChildren(.list = chips),
        # Hidden bridge: display:none via CSS class
        htmltools::div(
            class = "bio-chips-bridge",
            shiny::textInput(input_id, label = NULL,
                             value = paste(selected, collapse = ","), width = "100%")
        ),
        htmltools::tags$script(htmltools::HTML(js_code))
    )
}

#' JS message handler for bio_chips_reset (registered in global JS)
#'
#' Called by update_sidebar_inputs() when Reset is clicked.
#' Sent via session$sendCustomMessage("bio_chips_reset", list(...)).
#' The matching handler lives in custom.js (or inline in run_app.R).
#' @noRd
NULL  # handler registered in www/bio-chips-reset.js
