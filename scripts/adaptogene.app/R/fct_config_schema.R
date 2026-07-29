#' Config parameter schema for all pipeline tabs
#'
#' Each entry describes one config parameter: where it lives in the YAML,
#' how to display it, which tab it belongs to, and validation rules.
#'
#' Fields:
#'   key        dot-path into config (e.g. "sNMF.k_start")
#'   label      human-readable label
#'   tab        tab name: home | processing | prestructure | structure | gea |
#'                         gwas | gea_x_gwas | haplotype | maladaptation
#'   section    grouping label within the tab (subsection when `group` is set)
#'   group      optional top-level grouping label. When ANY entry on a tab sets
#'              `group`, that tab renders as flat named sections (one per group,
#'              `section` becomes a subsection label within it) instead of the
#'              default mandatory-core / Advanced-accordion split.
#'   type       numeric | text | select | checkbox | checkbox_invert | textarea |
#'                       method_table
#'              checkbox_invert displays/writes the negation of the stored value
#'              (e.g. a "Disable X" box backed by an `X.enabled` key)
#'   mandatory  TRUE = shown in core section; FALSE = shown in Advanced accordion
#'              (ignored for tabs using `group`)
#'   show_if    optional dot-key of a checkbox entry; this field only renders
#'              while that checkbox is checked (client-side conditionalPanel)
#'   help       tooltip / hint text shown below input
#'   placeholder hint text for empty text inputs
#'   min/max/step  for numeric type
#'   choices    for select type (character vector)
#'
#' @noRd
config_schema <- function() {
    # helper to build one entry compactly
    s <- function(key, label, tab, section, type, mandatory, ..., group = NULL, show_if = NULL) {
        args <- list(...)
        list(
            key         = key,
            label       = label,
            tab         = tab,
            section     = section,
            group       = group,
            type        = type,
            mandatory   = mandatory,
            show_if     = show_if,
            help        = args$help,
            placeholder = args$placeholder,
            min         = args$min,
            max         = args$max,
            step        = args$step,
            choices     = args$choices
        )
    }

    list(
        # ── HOME TAB (project input files only) ───────────────────────────────
        s("Input.dir",      "Input directory",    "home", "Input Files",
          "text",    TRUE,
          help = "Directory containing input files (relative to /pipeline/)",
          placeholder = "data/"),
        s("Input.vcf",      "VCF file",           "home", "Input Files",
          "text",    TRUE,
          help = "VCF filename inside input dir (.vcf or .vcf.gz)",
          placeholder = "samples.vcf.gz"),
        s("Input.metadata", "Metadata file",      "home", "Input Files",
          "text",    TRUE,
          help = "Tab-separated: site, sample, lat, lon, [phenotype columns…]",
          placeholder = "metadata.tsv"),
        s("Input.gff",      "GFF3 annotation",    "home", "Input Files",
          "text",    TRUE,
          help = "Optional. Enables gene annotation and GO enrichment.",
          placeholder = "annotation.gff3"),

        # ── PROCESSING TAB ────────────────────────────────────────────────────
        # Two top-level groups (SNP Filtering / Sample Filtering), no Advanced
        # accordion — every processing filter is a first-look decision, not a
        # rarely-touched knob. LD Pruning is a subsection within SNP Filtering.
        s("Filter.maf",         "Minor allele freq",  "processing", NULL,
          "numeric", TRUE, group = "SNP Filtering",
          min = 0, max = 0.5, step = 0.01,
          help = "Minimum MAF threshold (e.g. 0.05)"),
        s("Filter.snp_miss",    "SNP missingness",    "processing", NULL,
          "numeric", TRUE, group = "SNP Filtering",
          min = 0, max = 1, step = 0.01,
          help = "Maximum fraction of missing genotypes per SNP"),
        s("LD.window", "LD window (kb)",  "processing", "LD Pruning",
          "numeric", TRUE, group = "SNP Filtering",
          min = 1, step = 1,
          help = "Sliding window size for LD pruning in kilobases"),
        s("LD.step",   "LD step size",   "processing", "LD Pruning",
          "numeric", TRUE, group = "SNP Filtering",
          min = 1, step = 1,
          help = "Number of SNPs to slide per step"),
        s("LD.r2",     "LD r\u00b2 threshold", "processing", "LD Pruning",
          "numeric", TRUE, group = "SNP Filtering",
          min = 0, max = 1, step = 0.05,
          help = "Prune SNP pairs in LD above this r\u00b2 threshold"),
        s("Filter.sample_miss", "Sample missingness", "processing", NULL,
          "numeric", TRUE, group = "Sample Filtering",
          min = 0, max = 1, step = 0.01,
          help = "Maximum fraction of missing genotypes per sample (default 0.5)"),
        s("Filter.relatedness", "Max relatedness (IBS)", "processing", NULL,
          "numeric", FALSE, group = "Sample Filtering",
          min = 0, max = 1, step = 0.05,
          help = "Optional. Colors/counts related pairs (plink IBS allele-sharing) above this in the Processing relatedness histogram (blank = skip). Does NOT remove samples by itself \u2014 see 'Relatedness action' below. IBS is model-free \u2014 works for selfers and outcrossers alike; duplicates/clones cluster near 1.0 \u2014 set from the histogram, not a fixed outbred default."),
        s("Filter.relatedness_action", "Relatedness action", "processing", NULL,
          "select",  FALSE, group = "Sample Filtering",
          choices = c("keep", "remove"),
          help = "keep (default): only visualize/count related pairs, remove nothing. remove: drop the higher-missingness member of each pair above the threshold. Recommended: set the threshold, inspect the relatedness histogram, then switch to remove and re-run Processing."),

        # ── PRESTRUCTURE TAB ──────────────────────────────────────────────────
        s("sNMF.k_start", "K range start",  "prestructure", "sNMF",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Minimum K (ancestral populations) to test in cross-validation"),
        s("sNMF.k_end",   "K range end",    "prestructure", "sNMF",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Maximum K to test in cross-validation"),
        s("sNMF.repeats", "Repeats per K",  "prestructure", "sNMF",
          "numeric", TRUE,
          min = 1, max = 100, step = 1,
          help = "Number of independent sNMF runs per K for robust cross-entropy estimation"),

        # ── STRUCTURE TAB ─────────────────────────────────────────────────────
        s("sNMF.k_best",           "Best K",            "structure", "Structure",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Selected K for all downstream analysis. Review cross-entropy plot first."),
        s("Climate.enabled",       "Disable climate",   "structure", "Climate",
          "checkbox_invert", TRUE,
          help = "For datasets without coordinates. Check to run GWAS analysis only — skips WorldClim download, GEA, and maladaptation."),
        s("Climate.predictors",    "Climate predictors", "gea", "Climate",
          "bio_chips", TRUE,
          help = "Click to toggle variables for GEA. Review piemaps in Structure first to drop collinear or low-variance predictors."),
        # Optional — map
        s("Map.climate_extent", "Map extent",        "structure", "Map",
          "text",    TRUE,
          help = "Bounding box [min_lon,max_lon,min_lat,max_lat] or 'auto'",
          placeholder = "auto"),
        s("Map.gap",            "Auto extent gap",   "structure", "Map",
          "numeric", FALSE,
          min = 0, max = 10, step = 0.5,
          help = "Degrees of buffer around sample coordinates (auto extent only)"),
        s("Map.resolution",     "Resolution (min)",  "structure", "Map",
          "numeric", TRUE,
          min = 0.5, max = 10, step = 0.5,
          help = "WorldClim resolution in arc-minutes (lower = finer but slower)"),
        s("Map.zoom_extent",    "Zoom extent",       "structure", "Map",
          "text",    FALSE,
          help = "Optional zoom region for piemaps: xmin,xmax,ymin,ymax",
          placeholder = "NULL"),
        # Population stats — core (not Advanced): dependents reveal only once
        # calc_stats is checked, so the section stays compact when unused.
        s("Population.calc_stats",       "Calculate pop stats", "structure", "Pop Stats",
          "checkbox", TRUE,
          help = "Compute Tajima's D, Pi diversity, AMOVA, and IBD analysis"),
        s("Population.window_size",      "Window size (bp)",   "structure", "Pop Stats",
          "numeric",  TRUE,
          min = 1000, step = 1000,
          show_if = "Population.calc_stats",
          help = "Sliding window size in base pairs for Tajima's D and Pi diversity"),
        s("Population.custom_trait_file","Custom trait file",  "structure", "Pop Stats",
          "text",     TRUE,
          show_if = "Population.calc_stats",
          help = "Custom trait file for piemap sizing (NULL = use metadata phenotype columns)",
          placeholder = "NULL"),
        # Optional — piemap. Order: size/appearance first, then the pie->points
        # style switch, then label controls (label size gated on show_labels).
        s("Piemap.pie_scale",   "Pie size scale",    "structure", "Piemap",
          "numeric", FALSE,
          min = 0.1, max = 5, step = 0.1,
          help = "Multiplier for pie chart diameter (1.0 = default)"),
        s("Piemap.alpha",       "Pie transparency",  "structure", "Piemap",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.05,
          help = "Transparency of pie slices (0 = transparent, 1 = opaque)"),
        s("Piemap.use_points",  "Replace pies with points", "structure", "Piemap",
          "checkbox", FALSE,
          help = "Draw sample sites as simple dots instead of pie charts — clearer when dense sampling makes pies overlap. Affects on-demand haplotype-viz maps; the main Structure/Maladaptation/GWAS maps always emit both and switch live via the Points toggle in the app."),
        s("Piemap.show_labels", "Show pop labels",   "structure", "Piemap",
          "checkbox", FALSE,
          help = "Overlay population name labels on piemaps"),
        s("Piemap.label_size",  "Label size",        "structure", "Piemap",
          "numeric",  FALSE,
          min = 4, max = 20, step = 1,
          show_if = "Piemap.show_labels",
          help = "Font size for population labels"),
        # Optional — LD decay
        s("LDdecay.group_by",    "Group by",        "structure", "LD Decay",
          "select", TRUE,
          choices = c("site", "cluster"),
          help = "Variable to group samples by for LD decay curves"),
        s("LDdecay.min_samples", "Min samples",     "structure", "LD Decay",
          "numeric", TRUE,
          min = 2, step = 1,
          help = "Minimum group size; smaller groups are excluded"),
        s("LDdecay.max_distance","Max dist (kb)",   "structure", "LD Decay",
          "numeric", TRUE,
          min = 10, step = 50,
          help = "Maximum pairwise distance in kb for LD calculations"),
        s("LDdecay.scope",       "Scope",           "structure", "LD Decay",
          "select", TRUE,
          choices = c("both", "genome_wide", "per_chromosome"),
          help = "Compute LD decay genome-wide, per chromosome, or both"),

        # ── PREGEA TAB ────────────────────────────────────────────────────────
        # Hyperparameter-exploration mode (LD-pruned only) — sweeps LFMM-K,
        # EMMAX-#PC, and RDA Condition()-PC ladders + a shared spatial
        # variance-partitioning panel before committing to the expensive
        # full-SNP GEA run. Named `group`s mirror the config YAML's own
        # LFMM/EMMAX/RDA/Varpart/TransferGuard sub-blocks 1:1 (processing-tab
        # style: every entry always visible, no mandatory/Advanced split —
        # this tab's whole audience is a power user tuning hyperparameters).
        s("PreGEA.seed", "Random seed", "pregea", NULL,
          "numeric", TRUE, group = "General",
          min = 0, step = 1,
          help = "Seed for all preGEA permutation tests (RDA anova.cca, ordiR2step, varpart) — reproducibility only, not a tuning knob."),

        s("PreGEA.LFMM.enabled", "Run LFMM-K ladder", "pregea", NULL,
          "checkbox", TRUE, group = "LFMM K ladder",
          help = "Sweep LFMM K around sNMF.k_best with genomic_control OFF (below) so lambda is actually informative for K selection."),
        s("PreGEA.LFMM.k_offset_low", "K offset (low)", "pregea", NULL,
          "numeric", TRUE, group = "LFMM K ladder", show_if = "PreGEA.LFMM.enabled",
          min = -5, max = 0, step = 1,
          help = "Sweep starts at sNMF.k_best + this (usually negative)"),
        s("PreGEA.LFMM.k_offset_high", "K offset (high)", "pregea", NULL,
          "numeric", TRUE, group = "LFMM K ladder", show_if = "PreGEA.LFMM.enabled",
          min = 0, max = 10, step = 1,
          help = "Sweep ends at sNMF.k_best + this"),
        s("PreGEA.LFMM.k_min", "K hard floor", "pregea", NULL,
          "numeric", TRUE, group = "LFMM K ladder", show_if = "PreGEA.LFMM.enabled",
          min = 1, step = 1,
          help = "lfmm2() requires K >= 1 — the swept range is clamped to this floor"),
        s("PreGEA.LFMM.predictors", "Predictors", "pregea", NULL,
          "text", TRUE, group = "LFMM K ladder", show_if = "PreGEA.LFMM.enabled",
          help = "'auto' = Climate.predictors, or a comma-separated subset",
          placeholder = "auto"),
        s("PreGEA.LFMM.genomic_control", "Genomic control", "pregea", NULL,
          "checkbox", TRUE, group = "LFMM K ladder", show_if = "PreGEA.LFMM.enabled",
          help = "Keep OFF for K selection: lfmm2.test's genomic.control recalibration forces lambda toward 1 by construction, so a lambda computed WITH it on is not informative about which K is correct. The production GEA run uses genomic_control=TRUE separately — this only affects the ladder's own lambda-vs-K diagnostic."),
        s("PreGEA.LFMM.fdr", "FDR (hit-count series)", "pregea", NULL,
          "numeric", FALSE, group = "LFMM K ladder", show_if = "PreGEA.LFMM.enabled",
          min = 0, max = 1, step = 0.01,
          help = "qvalue FDR threshold for the hits-vs-K plot's first series"),
        s("PreGEA.LFMM.bonf_alpha", "Bonferroni alpha (hit-count series)", "pregea", NULL,
          "numeric", FALSE, group = "LFMM K ladder", show_if = "PreGEA.LFMM.enabled",
          min = 0, max = 1, step = 0.01,
          help = "Bonferroni alpha for the hits-vs-K plot's second series"),

        s("PreGEA.EMMAX.enabled", "Run EMMAX-#PC ladder", "pregea", NULL,
          "checkbox", TRUE, group = "EMMAX #PC ladder",
          help = "Sweep the number of PCA covariates on top of a fixed kinship matrix. n_pcs=0 is the kinship-only reference panel."),
        s("PreGEA.EMMAX.n_pcs_min", "#PCs (min)", "pregea", NULL,
          "numeric", TRUE, group = "EMMAX #PC ladder", show_if = "PreGEA.EMMAX.enabled",
          min = 0, step = 1,
          help = "0 = kinship-only reference panel"),
        s("PreGEA.EMMAX.n_pcs_max", "#PCs (max)", "pregea", NULL,
          "numeric", TRUE, group = "EMMAX #PC ladder", show_if = "PreGEA.EMMAX.enabled",
          min = 1, step = 1,
          help = "Keep small on tiny sample sizes"),
        s("PreGEA.EMMAX.kinship", "Kinship", "pregea", NULL,
          "select", TRUE, group = "EMMAX #PC ladder", show_if = "PreGEA.EMMAX.enabled",
          choices = c("BN", "IBS"),
          help = "Matches GEA.EMMAX.kinship semantics — held fixed across the whole ladder"),
        s("PreGEA.EMMAX.predictors", "Predictors", "pregea", NULL,
          "text", TRUE, group = "EMMAX #PC ladder", show_if = "PreGEA.EMMAX.enabled",
          help = "'auto' = Climate.predictors, or a comma-separated subset",
          placeholder = "auto"),
        s("PreGEA.EMMAX.lambda_tol", "Lambda tolerance", "pregea", NULL,
          "numeric", FALSE, group = "EMMAX #PC ladder", show_if = "PreGEA.EMMAX.enabled",
          min = 0, max = 1, step = 0.01,
          help = "|lambda - 1| tolerance defining the recommended 'knee' rung"),
        s("PreGEA.EMMAX.deflation_floor", "Deflation floor", "pregea", NULL,
          "numeric", FALSE, group = "EMMAX #PC ladder", show_if = "PreGEA.EMMAX.enabled",
          min = 0, max = 1, step = 0.01,
          help = "Reference line — lambda below this flags over-correction"),

        s("PreGEA.RDA.enabled", "Run RDA setup", "pregea", NULL,
          "checkbox", TRUE, group = "RDA setup",
          help = "Predictor collinearity screen + Condition()-PC ladder + ordiR2step forward selection path."),
        s("PreGEA.RDA.condition_pcs_min", "Condition PCs (min)", "pregea", NULL,
          "numeric", TRUE, group = "RDA setup", show_if = "PreGEA.RDA.enabled",
          min = 0, step = 1,
          help = "Population-structure PCs partialled out via Condition() — sweep floor"),
        s("PreGEA.RDA.condition_pcs_max", "Condition PCs (max)", "pregea", NULL,
          "numeric", TRUE, group = "RDA setup", show_if = "PreGEA.RDA.enabled",
          min = 1, step = 1,
          help = "Sweep ceiling"),
        s("PreGEA.RDA.collinearity_r", "Collinearity |r| pre-screen", "pregea", NULL,
          "numeric", FALSE, group = "RDA setup", show_if = "PreGEA.RDA.enabled",
          min = 0, max = 1, step = 0.05,
          help = "Predictor pairs above this |r| are pre-screened before VIF"),
        s("PreGEA.RDA.vif_max", "VIF max (post-fit)", "pregea", NULL,
          "numeric", FALSE, group = "RDA setup", show_if = "PreGEA.RDA.enabled",
          min = 1, step = 0.5,
          help = "vif.cca() cutoff applied after the |r| pre-screen"),
        s("PreGEA.RDA.min_predictors", "Min predictors kept", "pregea", NULL,
          "numeric", TRUE, group = "RDA setup", show_if = "PreGEA.RDA.enabled",
          min = 2, step = 1,
          help = "Never pre-screen below this — RDA and the covRob loadings estimator both need >= 2 predictors"),
        s("PreGEA.RDA.axis_alpha", "Axis significance alpha", "pregea", NULL,
          "numeric", FALSE, group = "RDA setup", show_if = "PreGEA.RDA.enabled",
          min = 0, max = 1, step = 0.01,
          help = "anova.cca per-axis significance threshold"),
        s("PreGEA.RDA.permutations", "Permutations", "pregea", NULL,
          "numeric", FALSE, group = "RDA setup", show_if = "PreGEA.RDA.enabled",
          min = 99, step = 100,
          help = "anova.cca permutation count (production default 999; SIMDATA-scale runs can use fewer)"),

        s("PreGEA.Varpart.enabled", "Run varpart / dbMEM", "pregea", NULL,
          "checkbox", TRUE, group = "Varpart & spatial (dbMEM)",
          help = "Shared spatial-vector artifact (adespatial::dbmem) + climate/structure/geography variance partitioning."),
        s("PreGEA.Varpart.response", "Response matrix", "pregea", NULL,
          "select", TRUE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          choices = c("pcs", "snps"),
          help = "'pcs' (genomic PCA scores, scalable) or 'snps' (raw genotypes)"),
        s("PreGEA.Varpart.response_var_cutoff", "Response variance cutoff", "pregea", NULL,
          "numeric", TRUE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          min = 0, max = 1, step = 0.05,
          help = "Cumulative variance of genomic PCs retained as Y (response = 'pcs' only)"),
        s("PreGEA.Varpart.response_max_pcs", "Response PCs (max)", "pregea", NULL,
          "numeric", FALSE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          min = 1, step = 1,
          help = "Cap on retained PCs regardless of variance cutoff"),
        s("PreGEA.Varpart.response_min_pcs", "Response PCs (min)", "pregea", NULL,
          "numeric", FALSE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          min = 1, step = 1,
          help = "Floor so varpart always has a multivariate Y"),
        s("PreGEA.Varpart.structure_table", "Structure covariate", "pregea", NULL,
          "select", TRUE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          choices = c("qmatrix", "none"),
          help = "'qmatrix' = sNMF Q at sNMF.k_best (last column dropped); 'none' skips the structure fraction"),
        s("PreGEA.Varpart.spatial_level", "Spatial level", "pregea", NULL,
          "select", TRUE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          choices = c("auto", "site", "sample"),
          help = "'auto' uses site-level MEMs when there are enough distinct sites (avoids pseudoreplication), else falls back to sample-level"),
        s("PreGEA.Varpart.min_sites", "Min distinct sites", "pregea", NULL,
          "numeric", TRUE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          min = 2, step = 1,
          help = "Below this, dbMEM writes a skip record instead of crashing"),
        s("PreGEA.Varpart.confounding_flag", "Flag climate/geography confounding", "pregea", NULL,
          "checkbox", FALSE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          help = "Auto-flag when the shared climate/geography fraction exceeds the larger unique fraction"),
        s("PreGEA.Varpart.ordir2step_pin", "ordiR2step p-in", "pregea", NULL,
          "numeric", FALSE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          min = 0, max = 1, step = 0.01,
          help = "Blanchet double stopping rule — permutation p-value threshold to add a term"),
        s("PreGEA.Varpart.ordir2step_permutations", "ordiR2step permutations", "pregea", NULL,
          "numeric", FALSE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          min = 99, step = 100,
          help = "Permutations per forward-selection step"),
        s("PreGEA.Varpart.r2_permutations", "R2permutations", "pregea", NULL,
          "numeric", FALSE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          min = 99, step = 100,
          help = "Blanchet stopping rule's R2 permutation count (spec default 1000; SIMDATA-scale runs can use fewer)"),
        s("PreGEA.Varpart.varpart_permutations", "varpart anova permutations", "pregea", NULL,
          "numeric", FALSE, group = "Varpart & spatial (dbMEM)", show_if = "PreGEA.Varpart.enabled",
          min = 99, step = 100,
          help = "anova.cca permutation count for each testable variance fraction"),

        s("PreGEA.TransferGuard.enabled", "Run transfer guard", "pregea", NULL,
          "checkbox", FALSE, group = "Transfer guard (opt-in)",
          help = "Re-estimates LFMM/EMMAX lambda on the FULL marker set at the recommended hyperparameters, before committing to the full GEA run. OFF by default — pulls the full-set imputation chain into this otherwise LD-pruned-only mode."),
        s("PreGEA.TransferGuard.lfmm_k", "LFMM K override", "pregea", NULL,
          "text", FALSE, group = "Transfer guard (opt-in)", show_if = "PreGEA.TransferGuard.enabled",
          help = "'auto' resolves from pregea_recommendations.tsv, or enter an explicit K",
          placeholder = "auto"),
        s("PreGEA.TransferGuard.emmax_n_pcs", "EMMAX #PCs override", "pregea", NULL,
          "text", FALSE, group = "Transfer guard (opt-in)", show_if = "PreGEA.TransferGuard.enabled",
          help = "'auto' resolves from pregea_recommendations.tsv, or enter an explicit #PCs",
          placeholder = "auto"),

        # ── GEA TAB ───────────────────────────────────────────────────────────
        s("GEA.configs", "GEA methods", "gea", "Methods",
          "method_table", TRUE,
          help = "Methods to run, with optional per-method hyperparameters (expand a row's Params to edit). Method list and hyperparameter defaults come from the pipeline's method registry (workflow/methods/gea.py) — read the preGEA ladders before overriding a default."),
        # NOTE: SNP clumping distance is no longer config-editable \u2014 it is a purely
        # interactive parameter in the GEA/GWAS filter bar (default 100 kb), kept in
        # sync with the region table/rectangles. See mod_gea.R / mod_gwas.R.
        s("GEA.wza.window_size", "WZA window size",    "gea", "WZA",
          "text",    FALSE,
          help = "Window size for Weighted-Z Analysis (shown in the WZA regime view). auto_genome_wide (default) uses genome-wide LD r\u00b2=0.2; auto_per_chromosome uses per-chr; or enter a fixed bp integer. Changing this requires re-running the GEA module.",
          placeholder = "auto_genome_wide"),
        s("GEA.wza.fallback_window_bp", "WZA fallback window (bp)", "gea", "WZA",
          "numeric", FALSE,
          min = 1000, step = 1000,
          help = "Window size (bp) used when LD decay data is unavailable (default: 10000 bp). Applies only in auto mode. Changing this requires re-running the GEA module."),
        # Gene annotation
        s("GFF.feature",   "GFF feature type", "gea", "Gene Annotation",
          "text",    FALSE,
          help = "Feature type to extract from GFF3 for gene overlaps",
          placeholder = "mRNA"),
        s("GFF.gene_name", "Gene name field",  "gea", "Gene Annotation",
          "text",    FALSE,
          help = "GFF attribute used as gene display name in output tables",
          placeholder = "description"),
        s("GFF.biotype",   "Biotype field",    "gea", "Gene Annotation",
          "text",    FALSE,
          help = "GFF attribute for gene biotype (used in filtering)",
          placeholder = "biotype"),
        s("GFF.go_field", "GO field in GFF",  "gea", "Gene Annotation",
          "text",    FALSE,
          help = "GFF attribute containing GO term IDs. NULL disables enrichment analysis.",
          placeholder = "NULL"),
        s("GEA.promoter_length",  "Promoter length (bp)","gea","Gene Annotation",
          "numeric", FALSE,
          min = 0, step = 1000,
          help = "Upstream distance from gene start counted as promoter region"),
        # Enrichment
        s("Enrichment.top_terms",        "Top GO terms",      "gea", "Enrichment",
          "numeric", FALSE,
          min = 5, max = 100, step = 5,
          help = "Number of enriched GO terms shown in dot/emap/cnetplot"),
        s("Enrichment.cnet_label",       "Cnet gene label",   "gea", "Enrichment",
          "select",  FALSE,
          choices = c("gene_id", "Name", "description"),
          help = "Which gene attribute to show as labels in the concept network plot"),

        # ── GWAS TAB ─────────────────────────────────────────────────────────
        s("GWAS.traits",           "Phenotype traits",    "gwas", "Traits",
          "pheno_chips", TRUE,
          help = "Click to toggle which phenotype traits to include in GWAS. Discovered from metadata columns 5+. Empty = all traits."),
        s("GWAS.configs",          "GWAS methods",        "gwas", "Methods",
          "method_table", TRUE,
          help = "Association methods for phenotype GWAS. Same options as GEA methods."),
        s("GWAS.missing_strategy", "Missing data strategy","gwas", "Methods",
          "select",       TRUE,
          choices = c("DROP", "MEAN", "MEDIAN"),
          help = "How to handle samples missing phenotype data. DROP removes them per trait."),
        # NOTE: SNP clumping distance is no longer config-editable \u2014 see GEA note above.
        s("GWAS.wza.window_size", "WZA window size",    "gwas", "WZA",
          "text",    FALSE,
          help = "Window size for Weighted-Z Analysis (shown in the WZA regime view). auto_genome_wide (default) uses genome-wide LD r\u00b2=0.2; auto_per_chromosome uses per-chr; or enter a fixed bp integer. Changing this requires re-running the GWAS module.",
          placeholder = "auto_genome_wide"),
        s("GWAS.wza.fallback_window_bp", "WZA fallback window (bp)", "gwas", "WZA",
          "numeric", FALSE,
          min = 1000, step = 1000,
          help = "Window size (bp) used when LD decay data is unavailable (default: 10000 bp). Applies only in auto mode. Changing this requires re-running the GWAS module."),
        s("GWAS.promoter_length",  "Promoter length (bp)","gwas", "Gene Annotation",
          "numeric", FALSE,
          min = 0, step = 1000,
          help = "Promoter length (bp) for gene annotation around significant regions."),

        # ── GEAxGWAS TAB ──────────────────────────────────────────────────────
        s("GEAxGWAS.pairwise.window_size", "Pairwise window (bp)","gea_x_gwas", "Pairwise Analysis",
          "numeric", FALSE,
          min = 0, step = 100000,
          help = "Window size for scoring GEA vs GWAS pairwise region overlap"),
        s("GEAxGWAS.pairwise.min_snps",    "Min pairwise SNPs",   "gea_x_gwas", "Pairwise Analysis",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Minimum shared significant SNPs to call a pairwise overlap"),

        # ── HAPLOTYPE TAB ──────────────────────────────────────────────────────
        s("haplotype.scan.regions_source", "Regions source",    "haplotype", "Scan Setup",
          "select",  TRUE,
          choices = c("gea", "gwas", "gea_x_gwas", "all", "custom"),
          help = "Which pipeline output to pull candidate regions from"),
        s("haplotype.epsilon_selected",    "Selected epsilon",  "haplotype", "Scan Setup",
          "text",    TRUE,
          help = "Epsilon chosen after reviewing clustree. Required before running haplotype mode.",
          placeholder = "NULL"),
        s("haplotype.scan.top_regions",       "Top regions",       "haplotype", "Scan Parameters",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Number of top regions to scan per source"),
        s("haplotype.scan.min_snps",          "Min SNPs per region","haplotype", "Scan Parameters",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Skip regions with fewer SNPs than this threshold"),
        s("haplotype.scan.min_group_size",    "Min group size (MGmin)","haplotype","Scan Parameters",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "DBSCAN MinPts parameter controlling cluster density"),
        s("haplotype.scan.min_haplotype_size","Min haplotype size","haplotype", "Scan Parameters",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Minimum number of samples to call a valid haplotype group"),
        s("haplotype.scan.epsilon_range",     "Epsilon range",     "haplotype", "Scan Parameters",
          "textarea", FALSE,
          help = "Comma-separated epsilon values to scan (e.g. 0.1,0.2,0.4,0.6,0.8)",
          placeholder = "0.01,0.1,0.2,0.4,0.6,0.8,1.0"),
        s("haplotype.scan.metadata_type",     "Metadata type",     "haplotype", "Scan Parameters",
          "select",  FALSE,
          choices = c("site", "cluster", "both"),
          help = "Metadata variable to use for haplotype group annotation"),
        s("haplotype.scan.regions_file",      "Custom regions file","haplotype","Scan Parameters",
          "text",    FALSE,
          help = "Path to custom regions TSV — only used when source = 'custom'",
          placeholder = "NULL"),

        # ── MALADAPTATION TAB ──────────────────────────────────────────────────
        s("Future.ssp",    "Climate scenario (SSP)","maladaptation","Future Climate",
          "select",  TRUE,
          choices = c("126", "245", "370", "585"),
          help = "Shared Socioeconomic Pathway (126 = low emissions, 585 = high emissions)"),
        s("Future.year",   "Projection period",    "maladaptation","Future Climate",
          "select",  TRUE,
          choices = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"),
          help = "Decade range for future climate projections"),
        s("Future.models", "Climate models",       "maladaptation","Future Climate",
          "text",    FALSE,
          help = "Comma-separated CMIP6 model names to average. Blank = use all available.",
          placeholder = "(all models)"),
        s("Maladaptation.snp_sets",
          "SNP sets to run",    "maladaptation","SNP Sets",
          "snp_set_picker", TRUE,
          help = "Saved sets from the GEA tab. Selected sets define which adaptive-SNP sets GF runs on."),
        s("Maladaptation.methods.gradient_forest.ntree",
          "Number of trees",    "maladaptation","Gradient Forest",
          "numeric", FALSE,
          min = 100, step = 100,
          help = "Number of trees per Gradient Forest run (more = more stable, slower)"),
        s("Maladaptation.methods.gradient_forest.cor_threshold",
          "R\u00b2 threshold",  "maladaptation","Gradient Forest",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.05,
          help = "Minimum R\u00b2 to include a climate predictor in the model"),
        s("Maladaptation.methods.gradient_forest.spatial_correction",
          "Spatial correction", "maladaptation","Gradient Forest",
          "select",  FALSE,
          choices = c("with", "without", "both"),
          help = "Include spatial eigenvectors (PCNM) as covariates. 'both' runs each SNP set as spatial AND nospatial."),
        s("Maladaptation.methods.gradient_forest.random_model",
          "Random model",       "maladaptation","Gradient Forest",
          "checkbox", FALSE,
          help = "Also build a random-SNP neutral model for comparison with adaptive model"),
        s("Maladaptation.methods.geometric_offset.k",
          "Latent factors (K)", "maladaptation","Geometric Offset",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Number of LFMM2 latent factors. Leave blank to use sNMF k_best."),
        s("Maladaptation.methods.geometric_offset.scale",
          "Scale environment",  "maladaptation","Geometric Offset",
          "checkbox", FALSE,
          help = "Standardise climate predictors before fitting LFMM2 (recommended).")
    )
}

#' Get config value by dot-path key
#'
#' @param config config list from read_project_config()
#' @param key_path dot-separated path (e.g. "snmf.k_start")
#' @return the value, or NULL if not found
#' @noRd
config_get_by_path <- function(config, key_path) {
    keys <- strsplit(key_path, "\\.")[[1]]
    val <- config
    for (k in keys) {
        if (is.null(val) || !is.list(val) || !k %in% names(val)) return(NULL)
        val <- val[[k]]
    }
    val
}

#' Set config value by dot-path key (returns modified config list)
#'
#' @param config config list
#' @param key_path dot-separated path (e.g. "snmf.k_start")
#' @param value value to set
#' @return modified config list
#' @noRd
config_set_by_path <- function(config, key_path, value) {
    keys <- strsplit(key_path, "\\.")[[1]]
    if (length(keys) == 1) {
        config[[keys]] <- value
        return(config)
    }
    # Recurse: ensure intermediate list exists
    top <- keys[1]
    rest <- paste(keys[-1], collapse = ".")
    if (is.null(config[[top]]) || !is.list(config[[top]])) config[[top]] <- list()
    config[[top]] <- config_set_by_path(config[[top]], rest, value)
    config
}

#' Coerce a sidebar input value back to the appropriate R type for YAML
#'
#' @param raw_value the raw value from input[[id]]
#' @param type schema type: numeric | text | select | checkbox | textarea | method_table
#' @return R value suitable for yaml::write_yaml()
#' @noRd
input_to_config_value <- function(raw_value, type) {
    switch(type,
        "numeric" = {
            n <- suppressWarnings(as.numeric(raw_value))
            if (is.na(n) || length(n) == 0) NULL else n
        },
        "checkbox" = {
            isTRUE(raw_value)
        },
        "checkbox_invert" = {
            !isTRUE(raw_value)
        },
        "textarea" = {
            # Comma-separated string -> list of trimmed values
            if (is.null(raw_value) || !nzchar(trimws(raw_value))) return(NULL)
            parts <- trimws(strsplit(as.character(raw_value), ",")[[1]])
            parts <- parts[nzchar(parts)]
            if (length(parts) == 0) return(NULL)
            # Try numeric coercion; fall back to character list
            nums <- suppressWarnings(as.numeric(parts))
            if (!anyNA(nums)) as.list(nums) else as.list(parts)
        },
        "method_table" = raw_value,  # read-only, pass through
        {
            # text / select / default
            if (is.null(raw_value)) return(NULL)
            val <- as.character(raw_value)
            if (!nzchar(trimws(val))) NULL else val
        }
    )
}

#' Get schema entries for a specific tab
#' @noRd
schema_for_tab <- function(tab_name) {
    Filter(function(e) e$tab == tab_name, config_schema())
}

#' Split schema entries for a tab into mandatory and optional
#' @noRd
schema_split <- function(entries) {
    list(
        mandatory = Filter(function(e) isTRUE(e$mandatory), entries),
        optional  = Filter(function(e) !isTRUE(e$mandatory), entries)
    )
}
