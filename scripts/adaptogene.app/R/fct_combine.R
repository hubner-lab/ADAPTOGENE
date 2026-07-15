# On-the-fly combine strategies for interactive Shiny exploration.
# Ports combine logic from scripts/combine_selected_snps.R using data.table
# (avoids GenomicRanges dependency in Shiny).

#' Okabe-Ito colour palette (8 colours, colour-blind safe)
#' @noRd
OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' Map traits to Okabe-Ito colours (sorted order, consistent with build_manhattan_plotly)
#'
#' @param traits character vector of trait names
#' @return named character vector: trait → hex colour
#' @noRd
trait_color_map <- function(traits) {
    traits_sorted <- sort(unique(traits))
    n <- length(traits_sorted)
    colors <- rep(OKABE_ITO, ceiling(n / length(OKABE_ITO)))[seq_len(n)]
    stats::setNames(colors, traits_sorted)
}

#' Plotly symbol names for methods (8 symbols, reused via recycling)
#' @noRd
METHOD_SHAPES <- c("circle", "triangle-up", "square", "diamond",
                   "cross", "x", "star", "triangle-down")

#' Map methods to plotly marker symbols (sorted order, stable across filtering)
#'
#' @param methods character vector of method names
#' @return named character vector: method → plotly symbol name
#' @noRd
method_shape_map <- function(methods) {
    methods_sorted <- sort(unique(methods))
    n <- length(methods_sorted)
    shapes <- rep(METHOD_SHAPES, ceiling(n / length(METHOD_SHAPES)))[seq_len(n)]
    stats::setNames(shapes, methods_sorted)
}

#' Normalize legacy and alias strategy names to current canonical names
#' @noRd
.normalize_strategy <- function(s) {
    switch(s,
        # Legacy pipeline names
        All          = "Union",
        Sum          = "Union",
        Overlap      = "Cross-method",
        MethodOverlap = "Cross-method per-trait",
        PairOverlap  = "Cross-method per-trait",
        s  # pass through "Union", "Cross-method", "Cross-method per-trait" unchanged
    )
}

#' Combine per-method sig SNPs using a specified strategy (on-the-fly in Shiny)
#'
#' @param sigsnps_list Named list of data.tables, one per method.
#'   Each must have columns: SNPID, chr, pos, pvalue, method, trait.
#' @param strategy One of "All", "Overlap", "MethodOverlap", or a method name.
#'   Legacy values "Sum" and "PairOverlap" are accepted and normalized.
#' @param gap Integer. Max bp distance for spatial overlap (Overlap/MethodOverlap).
#' @return Long-format data.table with SNPID, chr, pos, pvalue, method, trait.
#'   Empty data.table if nothing passes.
#' @noRd
combine_sigsnps <- function(sigsnps_list, strategy = "All", gap = 200000L) {
    if (length(sigsnps_list) == 0) return(.empty_sigsnps())

    strategy <- .normalize_strategy(strategy)

    # Drop empty methods
    sigsnps_list <- sigsnps_list[vapply(sigsnps_list, nrow, integer(1)) > 0]
    if (length(sigsnps_list) == 0) return(.empty_sigsnps())

    methods_vec <- names(sigsnps_list)
    gap <- as.integer(gap)

    result <- if (strategy %in% methods_vec) {
        # ── Single-method ───────────────────────────────────────────────────────
        data.table::copy(sigsnps_list[[strategy]])

    } else if (strategy == "Union") {
        # ── Union: all sig SNPs from all methods ────────────────────────────────
        data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)

    } else if (strategy == "Cross-method") {
        # ── Cross-method spatial overlap (trait-agnostic) ───────────────────────
        # A SNP passes if any partner exists in another method within gap bp.
        selected_ids <- character(0)
        for (n1 in methods_vec) {
            for (n2 in methods_vec) {
                if (n1 == n2) next
                ids <- .spatial_overlap_ids(sigsnps_list[[n1]], sigsnps_list[[n2]], gap)
                selected_ids <- c(selected_ids, ids)
            }
        }
        selected_ids <- unique(selected_ids)
        if (length(selected_ids) == 0) return(.empty_sigsnps())
        all_snps <- data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)
        all_snps[SNPID %in% selected_ids]

    } else if (strategy == "Cross-method per-trait") {
        # ── Cross-method spatial overlap, same-trait only ───────────────────────
        all_snps <- data.table::rbindlist(sigsnps_list, use.names = TRUE, fill = TRUE)
        all_traits <- unique(all_snps$trait)
        selected_rows <- list()
        for (n1 in methods_vec) {
            for (n2 in methods_vec) {
                if (n1 == n2) next
                for (bio in all_traits) {
                    m1_bio <- sigsnps_list[[n1]][trait == bio]
                    m2_bio <- sigsnps_list[[n2]][trait == bio]
                    if (nrow(m1_bio) == 0 || nrow(m2_bio) == 0) next
                    ids1 <- .spatial_overlap_ids(m1_bio, m2_bio, gap)
                    ids2 <- .spatial_overlap_ids(m2_bio, m1_bio, gap)
                    if (length(ids1) > 0)
                        selected_rows <- c(selected_rows, list(m1_bio[SNPID %in% ids1]))
                    if (length(ids2) > 0)
                        selected_rows <- c(selected_rows, list(m2_bio[SNPID %in% ids2]))
                }
            }
        }
        if (length(selected_rows) == 0) return(.empty_sigsnps())
        unique(data.table::rbindlist(selected_rows, use.names = TRUE, fill = TRUE))

    } else {
        stop(paste0("Unknown combine strategy: '", strategy,
                    "'. Valid values: Union, Cross-method, Cross-method per-trait, or a method name."))
    }

    if (is.null(result) || nrow(result) == 0) return(.empty_sigsnps())

    # Deduplicate — same SNP can appear in multiple methods/traits
    result <- unique(result, by = c("SNPID", "method", "trait"))

    # Add min_pvalue column for combined Manhattan y-placement
    # (background PNG y-axis spans the minimum p-value per SNP across all methods)
    result[, min_pvalue := min(pvalue, na.rm = TRUE), by = "SNPID"]

    result
}

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Find SNP IDs from dt1 that have a partner in dt2 within gap bp on same chr
#' @noRd
.spatial_overlap_ids <- function(dt1, dt2, gap) {
    if (nrow(dt1) == 0 || nrow(dt2) == 0) return(character(0))

    # Expand dt1 to windows [pos-gap, pos+gap]; dt2 stays as points [pos, pos]
    win <- data.table::data.table(
        chr   = as.character(dt1$chr),
        start = pmax(1L, as.integer(dt1$pos) - gap),
        end   = as.integer(dt1$pos) + gap,
        SNPID = dt1$SNPID
    )
    pts <- data.table::data.table(
        chr   = as.character(dt2$chr),
        start = as.integer(dt2$pos),
        end   = as.integer(dt2$pos),
        SNPID = dt2$SNPID
    )
    data.table::setkey(win, chr, start, end)
    data.table::setkey(pts, chr, start, end)
    ov <- data.table::foverlaps(pts, win, type = "within", nomatch = NULL)
    if (nrow(ov) == 0) return(character(0))
    # Return IDs from dt1 (window side)
    unique(ov$SNPID)
}

#' Empty sig SNPs skeleton
#' @noRd
.empty_sigsnps <- function() {
    data.table::data.table(
        SNPID  = character(),
        chr    = character(),
        pos    = integer(),
        pvalue = numeric(),
        method = character(),
        trait  = character()
    )
}

#' Empty sig SNPs skeleton with region_id column (for interactive combine output)
#' @noRd
.empty_sigsnps_assoc <- function() {
    data.table::data.table(
        SNPID     = character(),
        chr       = character(),
        pos       = integer(),
        pvalue    = numeric(),
        method    = character(),
        trait     = character(),
        region_id = character()
    )
}

# ── Shared UI / server helpers ────────────────────────────────────────────────

#' Extract default threshold type and value from project config.
#'
#' Reads the first entry in GEA.configs or GWAS.configs.
#' Returns list(type="bonf", value=0.05) as safe fallback.
#' @noRd
default_threshold <- function(config, module = MOD_GEA) {
    config_key <- if (module == MOD_GWAS) "GWAS" else "GEA"
    configs <- config_get(config, config_key, "configs", default = list())
    if (length(configs) == 0) return(list(type = "bonf", value = 0.05))

    first <- configs[[1]]
    type  <- if (is.list(first)) (first$adjust  %||% "bonf") else "bonf"
    value <- if (is.list(first)) suppressWarnings(as.numeric(first$threshold %||% 0.05))
             else 0.05
    if (is.na(value) || value <= 0) value <- 0.05
    list(type = type, value = value)
}

#' Sensible default threshold value for a given type.
#' @noRd
threshold_value_default_for_type <- function(type) {
    switch(type,
        bonf   = 0.05,
        qval   = 0.05,
        top    = 100,
        custom = 1e-5,
        0.05
    )
}

#' Is the threshold value valid for the given type?
#' @noRd
threshold_value_valid_for_type <- function(type, value) {
    if (is.null(value) || is.na(value)) return(FALSE)
    v <- suppressWarnings(as.numeric(value))
    if (is.na(v)) return(FALSE)
    switch(type,
        bonf   = v > 0 && v <= 1,
        qval   = v > 0 && v <= 1,
        top    = v >= 1,
        custom = v > 0,
        FALSE
    )
}

#' Build the always-present threshold bar (regime switch + threshold type/value + hint).
#'
#' Separate from the trait matrix so it is never gated on traits being present.
#'
#' @param ns           Shiny namespace function
#' @param input_prefix Character prefix for input IDs (e.g. "" or "gea_")
#' @param regime_value Logical: current WZA regime state
#' @param threshold_type_value  Character: "bonf"/"qval"/"top"/"custom"
#' @param threshold_value_value Numeric: current threshold value
#' @return tagList
#' @noRd
build_threshold_bar_ui <- function(ns, input_prefix = "",
                                   regime_value          = FALSE,
                                   threshold_type_value  = "bonf",
                                   threshold_value_value = 0.05,
                                   regime_context        = "gea") {
    pid <- function(name) ns(paste0(input_prefix, name))

    htmltools::div(
        class = "threshold-bar mb-2",
        # Regime switch row
        htmltools::div(
            class = "filter-row align-items-center mb-2",
            bslib::input_switch(
                pid("regime"), "WZA regime",
                value = isTRUE(regime_value)
            ),
            htmltools::span(
                class = "text-muted small ms-2",
                "Toggle to use Weighted-Z Analysis windows instead of per-SNP p-values"
            ),
            wza_window_note(regime_context)
        ),
        # Threshold type + value
        htmltools::div(
            class = "d-flex align-items-end gap-2",
            htmltools::div(
                class = "d-flex flex-column",
                htmltools::span("Significance threshold", class = "filter-label mb-1"),
                shiny::selectInput(
                    pid("threshold_type"),
                    label = NULL,
                    choices = c(
                        "Bonferroni"     = "bonf",
                        "FDR (qval)"     = "qval",
                        "Top N SNPs"     = "top",
                        "Custom (raw p)" = "custom"
                    ),
                    selected = threshold_type_value,
                    width = "150px"
                )
            ),
            htmltools::div(
                class = "d-flex flex-column",
                htmltools::span(" ", class = "filter-label mb-1"),
                shiny::numericInput(
                    pid("threshold_value"),
                    label = NULL,
                    value = threshold_value_value,
                    min = 0, step = NA,
                    width = "110px"
                )
            ),
            shiny::uiOutput(pid("threshold_hint"))
        )
    )
}

#' Build the trait×method matrix + strategy + clumping filter bar.
#'
#' Pure function — call inside renderUI. Handles its own namespacing via `ns`.
#' The threshold/regime controls live in a separate always-present build_threshold_bar_ui().
#'
#' @param ns                        Shiny namespace function from session$ns
#' @param traits                    Character vector of ALL trait names (full list, not just sig ones)
#' @param methods                   Character vector of method names
#' @param trait_colors              Named character vector: trait -> hex colour (from all_trait_names for stability)
#' @param combo_counts              Named list: "trait::method" -> integer SNP/window count
#' @param combo_thresholds          Named list: "trait::method" -> threshold p-value (from compute_method_thresholds)
#' @param default_strategy_value    Character scalar: "Union", "Cross-method", or "Cross-method per-trait"
#' @param snp_clumping_distance_value Integer scalar: current clumping distance (bp)
#' @param input_prefix              Character scalar: prefix for all Shiny input IDs inside this bar.
#'   Use "" (default) for a single filter bar; use e.g. "gea_" or "gwas_" when two bars
#'   coexist in the same module to avoid input ID collisions.
#' @return tagList (or a muted placeholder div when no data available)
#' @noRd
build_filter_bar_ui <- function(ns, traits, methods, trait_colors,
                                combo_counts, combo_thresholds = list(),
                                default_strategy_value,
                                snp_clumping_distance_value = 100000L,
                                input_prefix = "") {
    # Helper: produce an input id with optional prefix
    pid <- function(name) ns(paste0(input_prefix, name))
    # Guard only for truly empty project (no data at all)
    if (length(methods) == 0 || length(traits) == 0) {
        return(htmltools::div(
            class = "text-muted small fst-italic py-2",
            "No association results available."
        ))
    }

    # Build table header row: blank + method column headers
    header_cells <- lapply(methods, function(m) {
        htmltools::tags$th(
            class = "tm-col-header",
            `data-method` = m,
            m
        )
    })
    header_row <- htmltools::tags$tr(
        htmltools::tags$th(),  # blank corner
        header_cells
    )

    # Build body rows: one per trait
    body_rows <- lapply(traits, function(t) {
        dot_html <- paste0(
            '<span class="filter-trait-dot" style="background:', trait_colors[t],
            ';display:inline-block;width:8px;height:8px;border-radius:50%;',
            'margin-right:4px;"></span>'
        )
        row_header <- htmltools::tags$th(
            class = "tm-row-header",
            `data-trait` = t,
            htmltools::HTML(paste0(dot_html, htmltools::htmlEscape(t)))
        )
        cells <- lapply(methods, function(m) {
            key     <- paste0(t, "::", m)
            n       <- combo_counts[[key]] %||% 0L
            thr     <- combo_thresholds[[key]]
            thr_txt <- if (is.null(thr) || is.na(thr)) "n/a" else
                           formatC(thr, format = "e", digits = 1)
            if (n > 0) {
                htmltools::tags$td(
                    htmltools::tags$button(
                        class = "tm-cell tm-active",
                        `data-trait` = t,
                        `data-method` = m,
                        htmltools::tags$span(class = "tm-count", as.character(n)),
                        htmltools::tags$span(class = "tm-thr", thr_txt)
                    )
                )
            } else {
                htmltools::tags$td(
                    htmltools::tags$button(
                        class = "tm-cell tm-empty",
                        `data-trait` = t,
                        `data-method` = m,
                        htmltools::tags$span(class = "tm-count", "\u2014"),
                        htmltools::tags$span(class = "tm-thr", thr_txt)
                    )
                )
            }
        })
        htmltools::tags$tr(row_header, cells)
    })

    matrix_table <- htmltools::tags$table(
        class = "tm-matrix",
        htmltools::tags$thead(header_row),
        htmltools::tags$tbody(body_rows)
    )

    # Hidden input bridge updated by JS
    hidden_input <- shiny::textInput(
        pid("tm_selection"), label = NULL, value = ""
    )

    # JS: delegated click on matrix container
    container_id <- pid("tm_container")
    input_id     <- pid("tm_selection")
    js_code <- sprintf('
(function() {
    var container = document.getElementById("%s");
    if (!container) return;
    // window-keyed last-sent cache so syncSelection() is idempotent across renderUI
    // re-creations. Each re-run of the renderUI builds a fresh IIFE (new closure), so a
    // closure-local lastSent would reset to undefined — window persists across re-renders.
    if (!window.__tmSelLast) window.__tmSelLast = {};
    function syncSelection() {
        var active = container.querySelectorAll(".tm-cell.tm-active");
        var pairs = Array.from(active).map(function(el) {
            return el.dataset.trait + "::" + el.dataset.method;
        });
        var json = JSON.stringify(pairs);
        // Skip setInputValue when selection unchanged — prevents reactive storms on
        // renderUI re-runs caused by cold-load param restoration (no {priority:"event"};
        // default priority adds a second dedup layer for genuine user interactions).
        if (window.__tmSelLast["%s"] === json) return;
        window.__tmSelLast["%s"] = json;
        Shiny.setInputValue("%s", json);
    }
    container.addEventListener("click", function(e) {
        var cell = e.target.closest(".tm-cell");
        if (cell && !cell.classList.contains("tm-empty")) {
            cell.classList.toggle("tm-active");
            syncSelection();
            return;
        }
        var rh = e.target.closest(".tm-row-header");
        if (rh) {
            var trait = rh.dataset.trait;
            var cells = container.querySelectorAll(".tm-cell[data-trait=\\"" + trait + "\\"]:not(.tm-empty)");
            var allActive = Array.from(cells).every(function(c) { return c.classList.contains("tm-active"); });
            cells.forEach(function(c) { c.classList.toggle("tm-active", !allActive); });
            syncSelection();
            return;
        }
        var ch = e.target.closest(".tm-col-header");
        if (ch) {
            var method = ch.dataset.method;
            var cells = container.querySelectorAll(".tm-cell[data-method=\\"" + method + "\\"]:not(.tm-empty)");
            var allActive = Array.from(cells).every(function(c) { return c.classList.contains("tm-active"); });
            cells.forEach(function(c) { c.classList.toggle("tm-active", !allActive); });
            syncSelection();
        }
    });
    syncSelection();
})();
', container_id, input_id, input_id, input_id)

    strategy_choices <- c("Union", "Cross-method", "Cross-method per-trait")
    strategy_labels  <- setNames(
        strategy_choices,
        c(
            "Union — all sig SNPs from all methods",
            "Cross-method — SNPs significant in ≥2 methods within clumping distance",
            "Cross-method per-trait — same as above, same trait only"
        )
    )

    default_strat_norm <- .normalize_strategy(default_strategy_value)
    if (!default_strat_norm %in% strategy_choices) default_strat_norm <- "Union"

    htmltools::div(
        class = "manhattan-filter-bar",
        htmltools::div(
            class = "filter-row align-items-start",
            # Matrix
            htmltools::div(
                id    = container_id,
                class = "tm-container me-4",
                matrix_table,
                hidden_input
            ),
            # Strategy
            htmltools::div(
                class = "d-flex flex-column me-4",
                htmltools::span("Strategy", class = "filter-label mb-1"),
                shiny::radioButtons(
                    pid("combine_strategy"), label = NULL,
                    choices  = strategy_labels,
                    selected = default_strat_norm,
                    inline   = FALSE
                ),
                htmltools::span(
                    class = "text-muted small mt-1",
                    "Cross-method shows only SNPs/windows significant in ≥2 methods ",
                    "within the clumping distance. A hit significant in a single method ",
                    "appears only under Union."
                )
            ),
            # Clumping distance
            htmltools::div(
                class = "d-flex flex-column me-4",
                htmltools::span("Clumping distance (bp)", class = "filter-label mb-1"),
                shiny::numericInput(
                    pid("snp_clumping_distance"), label = NULL,
                    value = snp_clumping_distance_value, min = 1000L, step = 100000L,
                    width = "160px"
                ),
                htmltools::span(
                    class = "text-muted small mt-1",
                    "Merge distance for regions; also used for Cross-method overlap"
                )
            )
        ),
        htmltools::tags$script(htmltools::HTML(js_code))
    )
}

#' Compute interactive sig SNPs from per-method data + matrix selection
#'
#' Pure function — call inside reactive(). Filters, combines, and stamps region IDs.
#'
#' @param all_method_sigsnps Named list of data.tables (method -> dt with SNPID/chr/pos/pvalue/method/trait)
#' @param tm_selection_json  JSON string from tm_selection hidden input (list of "trait::method" pairs)
#' @param combo_counts       Named list: "trait::method" -> integer count (for fallback)
#' @param known_traits       Character vector of traits to restrict to (avoids cross-tab bleed)
#' @param strategy           Character: "Union", "Cross-method", "Cross-method per-trait", or method name
#' @param clumping_distance  Integer: SNP clumping distance in bp
#' @param project_name       Character: project name (for assign_region_ids)
#' @param module             Character: MOD_GEA, MOD_GWAS, or MOD_GEAXGWAS
#' @return data.table (SNPID/chr/pos/pvalue/method/trait/region_id) or .empty_sigsnps_assoc()
#' @noRd
compute_interactive_sigsnps <- function(all_method_sigsnps, tm_selection_json,
                                        combo_counts, known_traits,
                                        strategy, clumping_distance,
                                        project_name, module) {
    if (length(all_method_sigsnps) == 0) return(NULL)
    strategy <- .normalize_strategy(strategy)

    # All valid trait::method pairs (positive count, restricted to known_traits)
    all_valid_pairs <- names(combo_counts)[
        vapply(combo_counts, function(n) n > 0, logical(1))
    ]
    all_valid_pairs <- all_valid_pairs[
        sub("::.*", "", all_valid_pairs) %in% known_traits
    ]

    # Parse selected pairs from JSON; fall back to all valid.
    # An empty JSON array ("[]") from JS means the matrix rendered with no
    # tm-active cells (un-initialized state), not a deliberate user deselect.
    # Treat length-0 parsed result the same as NULL / "" to show all by default.
    selected_pairs <- if (is.null(tm_selection_json) || !nzchar(tm_selection_json)) {
        all_valid_pairs
    } else {
        parsed <- tryCatch(
            jsonlite::fromJSON(tm_selection_json),
            error = function(e) NULL
        )
        if (is.null(parsed) || length(parsed) == 0) all_valid_pairs else parsed
    }

    # Filter each method to selected traits for that method
    filtered <- lapply(names(all_method_sigsnps), function(m) {
        paired_traits <- sub("::.*", "",
            selected_pairs[grepl(paste0("::", m, "$"), selected_pairs, fixed = FALSE)])
        if (length(paired_traits) == 0) return(NULL)
        dt <- all_method_sigsnps[[m]]
        if (nrow(dt) == 0) return(NULL)
        dt[trait %in% paired_traits]
    })
    names(filtered) <- names(all_method_sigsnps)
    filtered <- Filter(function(x) !is.null(x) && nrow(x) > 0, filtered)

    if (length(filtered) == 0) return(.empty_sigsnps_assoc())

    combined_dt <- combine_sigsnps(filtered, strategy = strategy,
                                   gap = as.integer(clumping_distance))

    if (nrow(combined_dt) == 0) return(.empty_sigsnps_assoc())

    # region_id is stamped downstream from the live-computed regions (see
    # assign_region_ids_from_regions() / plotted_sigsnps in mod_gea.R / mod_gwas.R),
    # not from the static pipeline regions_combined.tsv — keeps the clicked SNP's
    # region in sync with the interactive region table/rectangles.
    combined_dt[, region_id := NA_character_]
    combined_dt
}
