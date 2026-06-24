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
#'   section    grouping label within the tab
#'   type       numeric | text | select | checkbox | textarea | method_table
#'   mandatory  TRUE = shown in core section; FALSE = shown in Advanced accordion
#'   help       tooltip / hint text shown below input
#'   placeholder hint text for empty text inputs
#'   min/max/step  for numeric type
#'   choices    for select type (character vector)
#'
#' @noRd
config_schema <- function() {
    # helper to build one entry compactly
    s <- function(key, label, tab, section, type, mandatory, ...) {
        args <- list(...)
        list(
            key         = key,
            label       = label,
            tab         = tab,
            section     = section,
            type        = type,
            mandatory   = mandatory,
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
        s("Filter.maf",         "Minor allele freq",  "processing", "Filtering",
          "numeric", TRUE,
          min = 0, max = 0.5, step = 0.01,
          help = "Minimum MAF threshold (e.g. 0.05)"),
        s("Filter.snp_miss",    "SNP missingness",    "processing", "Filtering",
          "numeric", TRUE,
          min = 0, max = 1, step = 0.01,
          help = "Maximum fraction of missing genotypes per SNP"),
        s("Filter.sample_miss", "Sample missingness", "processing", "Filtering",
          "numeric", TRUE,
          min = 0, max = 1, step = 0.01,
          help = "Maximum fraction of missing genotypes per sample (default 0.5)"),
        s("LD.window", "LD window (kb)",  "processing", "LD Pruning",
          "numeric", TRUE,
          min = 1, step = 1,
          help = "Sliding window size for LD pruning in kilobases"),
        s("LD.step",   "LD step size",   "processing", "LD Pruning",
          "numeric", TRUE,
          min = 1, step = 1,
          help = "Number of SNPs to slide per step"),
        s("LD.r2",     "LD r\u00b2 threshold", "processing", "LD Pruning",
          "numeric", TRUE,
          min = 0, max = 1, step = 0.05,
          help = "Prune SNP pairs in LD above this r\u00b2 threshold"),

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
          "numeric", FALSE,
          min = 1, max = 100, step = 1,
          help = "Number of independent sNMF runs per K for robust cross-entropy estimation"),

        # ── STRUCTURE TAB ─────────────────────────────────────────────────────
        s("sNMF.k_best",           "Best K",            "structure", "Structure",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Selected K for all downstream analysis. Review cross-entropy plot first."),
        s("Climate.enabled",       "Enable climate",    "structure", "Climate",
          "checkbox", TRUE,
          help = "Download WorldClim data and enable GEA / maladaptation analysis"),
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
        # Optional — population stats
        s("Population.calc_stats",       "Calculate pop stats", "structure", "Pop Stats",
          "checkbox", FALSE,
          help = "Compute Tajima's D, Pi diversity, AMOVA, and IBD analysis"),
        s("Population.window_size",      "Window size (bp)",   "structure", "Pop Stats",
          "numeric",  FALSE,
          min = 1000, step = 1000,
          help = "Sliding window size in base pairs for Tajima's D and Pi diversity"),
        s("Population.custom_trait_file","Custom trait file",  "structure", "Pop Stats",
          "text",     FALSE,
          help = "Custom trait file for piemap sizing (NULL = use metadata phenotype columns)",
          placeholder = "NULL"),
        # Optional — piemap
        s("Piemap.alpha",       "Pie transparency",  "structure", "Piemap",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.05,
          help = "Transparency of pie slices (0 = transparent, 1 = opaque)"),
        s("Piemap.show_labels", "Show pop labels",   "structure", "Piemap",
          "checkbox", FALSE,
          help = "Overlay population name labels on piemaps"),
        s("Piemap.label_size",  "Label size",        "structure", "Piemap",
          "numeric",  FALSE,
          min = 4, max = 20, step = 1,
          help = "Font size for population labels"),
        s("Piemap.pie_scale",   "Pie size scale",    "structure", "Piemap",
          "numeric", FALSE,
          min = 0.1, max = 5, step = 0.1,
          help = "Multiplier for pie chart diameter (1.0 = default)"),
        s("Piemap.use_points",  "Use points",        "structure", "Piemap",
          "checkbox", FALSE,
          help = "Replace pie charts with single colored points (for large datasets)"),
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

        # ── GEA TAB ───────────────────────────────────────────────────────────
        s("GEA.configs", "GEA methods", "gea", "Methods",
          "method_table", TRUE,
          help = "Methods to run: EMMAX, LFMM, GLM, MLM, CMLM, ECMLM, SUPER, MLMM, FarmCPU, BLINK"),
        s("GFF.go_field", "GO field in GFF",  "gea", "Methods",
          "text",    FALSE,
          help = "GFF attribute containing GO term IDs. NULL disables enrichment analysis.",
          placeholder = "NULL"),
        # Results processing
        s("GEA.snp_clumping_distance", "SNP clumping distance", "gea", "Results",
          "text",    FALSE,
          help = "Single-linkage merge distance: SNPs within this distance form one region; region bounds are padded by this distance. auto_per_chromosome (default) uses per-chr LD r\u00b2=0.2; auto_genome_wide uses genome-wide; or enter a fixed bp value.",
          placeholder = "auto_per_chromosome"),
        s("GEA.clumping_r2_threshold", "LD r\u00b2 threshold",  "gea", "Results",
          "numeric", FALSE,
          min = 0.05, max = 0.95, step = 0.05,
          help = "r\u00b2 level at which LD is considered background. Hill-Weir curve is inverted at this value to get the clumping distance (default: 0.2)."),
        s("GEA.ld_decay_group",        "LD decay group",        "gea", "Results",
          "text",    FALSE,
          help = "Which sample group from the LD decay table to use for distance computation (default: All).",
          placeholder = "All"),
        s("GEA.promoter_length",  "Promoter length (bp)","gea","Results",
          "numeric", FALSE,
          min = 0, step = 1000,
          help = "Upstream distance from gene start counted as promoter region"),
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
        s("GWAS.snp_clumping_distance", "SNP clumping distance", "gwas", "Results",
          "text",    FALSE,
          help = "Single-linkage merge distance: SNPs within this distance form one region; region bounds are padded by this distance. auto_per_chromosome (default) uses per-chr LD r\u00b2=0.2; auto_genome_wide uses genome-wide; or enter a fixed bp value.",
          placeholder = "auto_per_chromosome"),
        s("GWAS.clumping_r2_threshold", "LD r\u00b2 threshold",  "gwas", "Results",
          "numeric", FALSE,
          min = 0.05, max = 0.95, step = 0.05,
          help = "r\u00b2 level at which LD is considered background. Hill-Weir curve is inverted at this value (default: 0.2)."),
        s("GWAS.ld_decay_group",        "LD decay group",        "gwas", "Results",
          "text",    FALSE,
          help = "Which sample group from the LD decay table to use for distance computation (default: All).",
          placeholder = "All"),
        s("GWAS.promoter_length",  "Promoter length (bp)","gwas", "Results",
          "numeric", FALSE,
          min = 0, step = 1000,
          help = "Promoter length (bp) for gene annotation around significant regions."),
        s("GWAS.wza.window_size", "WZA window size",    "gwas", "WZA",
          "text",    FALSE,
          help = "Window size for Weighted-Z Analysis (shown in the WZA regime view). auto_genome_wide (default) uses genome-wide LD r\u00b2=0.2; auto_per_chromosome uses per-chr; or enter a fixed bp integer. Changing this requires re-running the GWAS module.",
          placeholder = "auto_genome_wide"),
        s("GWAS.wza.fallback_window_bp", "WZA fallback window (bp)", "gwas", "WZA",
          "numeric", FALSE,
          min = 1000, step = 1000,
          help = "Window size (bp) used when LD decay data is unavailable (default: 10000 bp). Applies only in auto mode. Changing this requires re-running the GWAS module."),

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
          help = "Also build a random-SNP neutral model for comparison with adaptive model")
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
