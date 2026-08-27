# lab_90_dispatch.R — masks module UIs with variant dispatchers.
#
# Sourced LAST (name sorts after the variant files). Captures each app original
# BEFORE masking, so ?v=orig remains a genuine untouched control to compare
# against rather than a reimplementation of it.

# ── Processing ───────────────────────────────────────────────────────────────
.lab_orig_mod_processing_ui <- mod_processing_ui

mod_processing_ui <- function(id) {
    switch(lab_variant(),
        a    = mod_processing_ui_a(id),
        b    = mod_processing_ui_b(id),
        c    = mod_processing_ui_c(id),
        .lab_orig_mod_processing_ui(id)
    )
}

# ── PreStructure ─────────────────────────────────────────────────────────────
.lab_orig_mod_prestructure_ui     <- mod_prestructure_ui
.lab_orig_mod_prestructure_server <- mod_prestructure_server

mod_prestructure_ui <- function(id) {
    switch(lab_variant(),
        a = mod_prestructure_ui_a(id),
        .lab_orig_mod_prestructure_ui(id)
    )
}

# The real server is called UNCHANGED; the lab only ADDS the outputs it does not
# define. Same id -> same namespace, so extras can read the real server's inputs.
#
# `...` is forwarded on every wrapper: mod_gea_server() takes run_trigger/module
# beyond (id, project_data), and app_server passes them BY NAME. A two-arg
# wrapper would drop them silently.
mod_prestructure_server <- function(id, project_data, ...) {
    .lab_orig_mod_prestructure_server(id, project_data, ...)
    if (identical(lab_variant(), "a")) lab_prestructure_extras(id, project_data)
    invisible(NULL)
}

# ── Structure ────────────────────────────────────────────────────────────────
.lab_orig_mod_structure_ui     <- mod_structure_ui
.lab_orig_mod_structure_server <- mod_structure_server

mod_structure_ui <- function(id) {
    switch(lab_variant(),
        a = mod_structure_ui_a(id),
        .lab_orig_mod_structure_ui(id)
    )
}

mod_structure_server <- function(id, project_data, ...) {
    .lab_orig_mod_structure_server(id, project_data, ...)
    if (identical(lab_variant(), "a")) lab_structure_extras(id, project_data)
    invisible(NULL)
}

# ── Environmental (climate) ──────────────────────────────────────────────────
.lab_orig_mod_climate_ui     <- mod_climate_ui
.lab_orig_mod_climate_server <- mod_climate_server

mod_climate_ui <- function(id) {
    switch(lab_variant(),
        a = mod_climate_ui_a(id),
        .lab_orig_mod_climate_ui(id)
    )
}

mod_climate_server <- function(id, project_data, ...) {
    .lab_orig_mod_climate_server(id, project_data, ...)
    if (identical(lab_variant(), "a")) lab_climate_extras(id, project_data)
    invisible(NULL)
}

# ── Phenotypic (traits) ──────────────────────────────────────────────────────
.lab_orig_mod_traits_ui     <- mod_traits_ui
.lab_orig_mod_traits_server <- mod_traits_server

mod_traits_ui <- function(id) {
    switch(lab_variant(),
        a = mod_traits_ui_a(id),
        .lab_orig_mod_traits_ui(id)
    )
}

mod_traits_server <- function(id, project_data, ...) {
    .lab_orig_mod_traits_server(id, project_data, ...)
    if (identical(lab_variant(), "a")) lab_traits_extras(id, project_data)
    invisible(NULL)
}

# ── GEA ──────────────────────────────────────────────────────────────────────
.lab_orig_mod_gea_ui     <- mod_gea_ui
.lab_orig_mod_gea_server <- mod_gea_server

mod_gea_ui <- function(id) {
    switch(lab_variant(),
        a = mod_gea_ui_a(id),
        .lab_orig_mod_gea_ui(id)
    )
}

mod_gea_server <- function(id, project_data, ...) {
    .lab_orig_mod_gea_server(id, project_data, ...)
    if (identical(lab_variant(), "a")) lab_gea_extras(id, project_data)
    invisible(NULL)
}

# ── GWAS ─────────────────────────────────────────────────────────────────────
.lab_orig_mod_gwas_ui     <- mod_gwas_ui
.lab_orig_mod_gwas_server <- mod_gwas_server

mod_gwas_ui <- function(id) {
    switch(lab_variant(),
        a = mod_gwas_ui_a(id),
        .lab_orig_mod_gwas_ui(id)
    )
}

mod_gwas_server <- function(id, project_data, ...) {
    .lab_orig_mod_gwas_server(id, project_data, ...)
    if (identical(lab_variant(), "a")) lab_gwas_extras(id, project_data)
    invisible(NULL)
}

# ── GEAxGWAS ─────────────────────────────────────────────────────────────────
.lab_orig_mod_gea_x_gwas_ui     <- mod_gea_x_gwas_ui
.lab_orig_mod_gea_x_gwas_server <- mod_gea_x_gwas_server

mod_gea_x_gwas_ui <- function(id) {
    switch(lab_variant(),
        a = mod_gea_x_gwas_ui_a(id),
        .lab_orig_mod_gea_x_gwas_ui(id)
    )
}

mod_gea_x_gwas_server <- function(id, project_data, ...) {
    .lab_orig_mod_gea_x_gwas_server(id, project_data, ...)
    if (identical(lab_variant(), "a")) lab_geaxgwas_extras(id, project_data)
    invisible(NULL)
}

# ── PreGEA ───────────────────────────────────────────────────────────────────
.lab_orig_mod_pregea_ui     <- mod_pregea_ui
.lab_orig_mod_pregea_server <- mod_pregea_server

mod_pregea_ui <- function(id) {
    switch(lab_variant(),
        a = mod_pregea_ui_a(id),
        .lab_orig_mod_pregea_ui(id)
    )
}

mod_pregea_server <- function(id, project_data, ...) {
    .lab_orig_mod_pregea_server(id, project_data, ...)
    if (identical(lab_variant(), "a")) lab_pregea_extras(id, project_data)
    invisible(NULL)
}

# ── Maladaptation ────────────────────────────────────────────────────────────
.lab_orig_mod_maladaptation_ui     <- mod_maladaptation_ui
.lab_orig_mod_maladaptation_server <- mod_maladaptation_server

mod_maladaptation_ui <- function(id) {
    switch(lab_variant(),
        a = mod_maladaptation_ui_a(id),
        .lab_orig_mod_maladaptation_ui(id)
    )
}

mod_maladaptation_server <- function(id, project_data, ...) {
    .lab_orig_mod_maladaptation_server(id, project_data, ...)
    if (identical(lab_variant(), "a")) lab_malad_extras(id, project_data)
    invisible(NULL)
}
