# lab_40_variant_c.R — Variant C: NARRATIVE FLOW
#
# QC laid out as the two-rail pipeline it actually is. An always-visible SPINE
# (KPIs + both rails side by side) over a STAGE STEPPER you click to pull one
# stage's evidence.
#
# The two rails cannot overlap by construction: filter_vcf applies --keep and
# --maf/--geno in a single plink call, so sample selection is already final
# before any SNP filter runs. filtering_attrition.png is the SNP rail ONLY --
# the sample drops (51 -> 50 -> 47) live nowhere but filtering_table's
# n_samples column, which is why the table gets equal billing rather than being
# a footnote under the figure.
#
# Stops are THRESHOLD DECISION POINTS, not stages-that-fired: `After het outlier
# removal` is conditional and absent from SIMDATA's table, but the server binds
# all eight cards unconditionally, so a stop cannot be dropped without breaking
# the DOM contract. Whether a stage fired is answered by its note badge going
# NULL and by the row's absence from filtering_table.
#
# navset_card_underline keeps every .tab-pane in the DOM (inactive ones are just
# display:none), so the contract holds while hidden stages never read their PNG.

.lab_c_stop <- function(n, label, icon, rail = c("sample", "snp")) {
    rail <- match.arg(rail)
    htmltools::span(
        class = paste0("lab-c-stop lab-c-rail-", rail),
        htmltools::span(class = "lab-c-stop-num", n),
        bsicons::bs_icon(icon, class = "lab-c-stop-icon"),
        htmltools::span(class = "lab-c-stop-label", label)
    )
}

.lab_c_phase <- function(label, icon) {
    bslib::nav_item(htmltools::span(
        class = "lab-c-phase-sep", bsicons::bs_icon(icon), " ", label
    ))
}

mod_processing_ui_c <- function(id) {
    ns <- shiny::NS(id)

    lab_root("c",

        lab_kpi_row(ns),

        # ── SPINE: the two rails, equal billing ───────────────────────────────
        bslib::layout_columns(
            # 4/8, not 7/5: the attrition image is height-capped in the spine,
            # so a wide card just strands it in empty space. A narrower column
            # lets the image nearly fill its card while the table -- which is
            # the SAMPLE rail and needs the width -- gets the rest.
            col_widths = c(4, 8),
            class = "lab-c-spine",
            fill = FALSE, fillable = FALSE,
            gap = "0.6rem",
            mod_image_card_ui(ns("attrition")),
            lab_filtering_summary_card(ns)
        ),

        # ── STAGE STEPPER ─────────────────────────────────────────────────────
        bslib::navset_card_underline(
            id = ns("stage_stepper"),
            selected = "sample_miss",
            title = htmltools::div(
                class = "lab-c-stepper-title",
                bsicons::bs_icon("diagram-3"), " Filter stages ",
                htmltools::span(class = "text-muted",
                                "— click a stop to inspect its evidence")
            ),

            .lab_c_phase("Sample QC", "people-fill"),
            bslib::nav_panel(
                title = .lab_c_stop(1, "Sample missingness", "person-x", "sample"),
                value = "sample_miss",
                mod_image_card_ui(ns("sample_miss"))
            ),
            bslib::nav_panel(
                title = .lab_c_stop(2, "Het outliers", "activity", "sample"),
                value = "het",
                bslib::layout_columns(
                    col_widths = c(7, 5), fill = FALSE, fillable = FALSE,
                    mod_image_card_ui(ns("het_miss")),
                    lab_het_table_card(ns)
                )
            ),
            bslib::nav_panel(
                title = .lab_c_stop(3, "Relatedness", "diagram-2", "sample"),
                value = "relatedness",
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("relatedness")),
                    mod_image_card_ui(ns("relatedness_mds"))
                )
            ),

            # sample selection freezes here — everything downstream is SNP-only
            .lab_c_phase("SNP QC", "database-fill"),
            bslib::nav_panel(
                title = .lab_c_stop(4, "MAF", "bar-chart-line", "snp"),
                value = "maf",
                mod_image_card_ui(ns("maf"))
            ),
            bslib::nav_panel(
                title = .lab_c_stop(5, "SNP missingness", "grid-3x3-gap", "snp"),
                value = "snp_miss",
                mod_image_card_ui(ns("snp_miss"))
            ),
            bslib::nav_panel(
                # snp_density plots dens_raw AND dens_filt, so it is the net
                # before/after portrait -- it belongs at the survivors end.
                title = .lab_c_stop(6, "LD pruning → survivors", "scissors", "snp"),
                value = "survivors",
                mod_image_card_ui(ns("snp_density"))
            )
        ),

        lab_depth_below_fold(ns, open = FALSE)
    )
}
