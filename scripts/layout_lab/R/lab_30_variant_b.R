# lab_30_variant_b.R — Variant B: ANALYST WORKBENCH
#
# One figure large in a focus canvas, the other seven live as a 2-up filmstrip,
# both tables docked right. Click a thumbnail to swap what is in the canvas.
#
# Why nothing is ever hidden: a navset inactive pane / conditionalPanel /
# collapsed accordion is display:none, which SUSPENDS renderImage -- the Step 0
# smoke test confirmed the image is not even transmitted. A filmstrip you
# navigate BY LOOKING needs live mini-plots, so B renders all eight cards
# always-visible and makes thumbnail-vs-focus a pure CSS size change driven by
# one data-focus attribute. The trade: all eight PNGs decode up front.
#
# Reuses mod_processing_server() unchanged: every id below is byte-identical.

.lab_b_fig <- function(ns, fig_id, domain, label) {
    htmltools::div(
        class      = paste0("lab-b-fig dom-", domain),
        `data-fig` = fig_id,
        # A real <button> gives native focus ring, Tab order and Enter/Space
        # with no roving-tabindex code.
        htmltools::tags$button(
            class = "lab-b-focus-btn", type = "button",
            `aria-label` = paste("Focus", label)
        ),
        mod_image_card_ui(ns(fig_id))
    )
}

mod_processing_ui_b <- function(id) {
    ns <- shiny::NS(id)

    lab_root("b",

        htmltools::tags$script(htmltools::HTML(
            "(function(){if(window.__labB)return;window.__labB=1;
             document.addEventListener('click',function(e){
               var b=e.target.closest('.lab-b-focus-btn'); if(!b)return;
               var f=b.parentElement, s=f.closest('.lab-b-stage'),
                   id=f.getAttribute('data-fig');
               if(s&&id)s.setAttribute('data-focus',id);});})();"
        )),

        lab_kpi_row(ns),

        htmltools::div(
            class = "lab-b-body",

            # STAGE — focus canvas (col 1) + 2-up filmstrip (cols 2-3).
            # Every figure is always display:block; only its grid placement and
            # size change.
            htmltools::div(
                class = "lab-b-stage", `data-focus` = "attrition",
                .lab_b_fig(ns, "attrition",       "attr", "Filtering Attrition"),
                .lab_b_fig(ns, "sample_miss",     "samp", "Sample Missingness"),
                .lab_b_fig(ns, "het_miss",        "samp", "Heterozygosity vs Missingness"),
                .lab_b_fig(ns, "relatedness",     "samp", "Relatedness (IBS)"),
                .lab_b_fig(ns, "relatedness_mds", "samp", "Relatedness MDS (IBS)"),
                .lab_b_fig(ns, "snp_miss",        "snp",  "SNP Missingness"),
                .lab_b_fig(ns, "maf",             "snp",  "MAF Distribution"),
                .lab_b_fig(ns, "snp_density",     "snp",  "SNP Density by Chromosome")
            ),

            # DOCK — both tables always visible (display:none would suspend them)
            htmltools::div(
                class = "lab-b-dock",
                lab_filtering_summary_card(ns),
                lab_het_table_card(ns)
            )
        ),

        # Depth below the fold. OPEN, not collapsed: B's premise is that nothing
        # is hidden, and a closed panel would leave the card blank until opened.
        lab_depth_below_fold(ns, open = TRUE)
    )
}
