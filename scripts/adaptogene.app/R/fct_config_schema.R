#' Config parameter schema for all pipeline tabs
#'
#' Each entry describes one config parameter: where it lives in the YAML,
#' how to display it, which tab it belongs to, and validation rules.
#'
#' Fields:
#'   key        dot-path into config (e.g. "snmf.k_start")
#'   label      human-readable label
#'   tab        tab name: home | structure | structure_k | association |
#'                         phenotype | overlapping | haplotype | maladaptation
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
        # ── HOME TAB (Processing mode) ─────────────────────────────────────────
        # Mandatory — project identity
        s("project_name",   "Project name",       "home", "Project",
          "text",    TRUE,
          help = "Output directory prefix. Results go to {name}_results/",
          placeholder = "MYPROJECT"),
        s("cpu",            "CPU cores",           "home", "Project",
          "numeric", TRUE,
          min = 1, max = 64, step = 1,
          help = "Number of parallel Snakemake jobs"),
        # Mandatory — input files
        s("input.dir",      "Input directory",    "home", "Input Files",
          "text",    TRUE,
          help = "Directory containing input files (relative to /pipeline/)",
          placeholder = "data/"),
        s("input.vcf",      "VCF file",           "home", "Input Files",
          "text",    TRUE,
          help = "VCF filename inside input dir (.vcf or .vcf.gz)",
          placeholder = "samples.vcf.gz"),
        s("input.metadata", "Metadata file",      "home", "Input Files",
          "text",    TRUE,
          help = "Tab-separated: site, sample, lat, lon, [phenotype columns…]",
          placeholder = "metadata.tsv"),
        s("input.gff",      "GFF3 annotation",    "home", "Input Files",
          "text",    FALSE,
          help = "Optional. Enables gene annotation and GO enrichment.",
          placeholder = "annotation.gff3"),
        # Optional — filtering
        s("filter.maf",         "Minor allele freq",  "home", "Filtering",
          "numeric", FALSE,
          min = 0, max = 0.5, step = 0.01,
          help = "Minimum MAF threshold (e.g. 0.05)"),
        s("filter.snp_miss",    "SNP missingness",    "home", "Filtering",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.01,
          help = "Maximum fraction of missing genotypes per SNP"),
        s("filter.sample_miss", "Sample missingness", "home", "Filtering",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.01,
          help = "Maximum fraction of missing genotypes per sample (default 0.5)"),
        # Optional — LD pruning
        s("ld.window", "LD window (kb)",  "home", "LD Pruning",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Sliding window size for LD pruning in kilobases"),
        s("ld.step",   "LD step size",   "home", "LD Pruning",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Number of SNPs to slide per step"),
        s("ld.r2",     "LD r\u00b2 threshold", "home", "LD Pruning",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.05,
          help = "Prune SNP pairs in LD above this r\u00b2 threshold"),

        # ── STRUCTURE TAB ──────────────────────────────────────────────────────
        s("snmf.k_start", "K range start",  "structure", "sNMF",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Minimum K (ancestral populations) to test in cross-validation"),
        s("snmf.k_end",   "K range end",    "structure", "sNMF",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Maximum K to test in cross-validation"),
        s("snmf.repeats", "Repeats per K",  "structure", "sNMF",
          "numeric", FALSE,
          min = 1, max = 100, step = 1,
          help = "Number of independent sNMF runs per K for robust cross-entropy estimation"),

        # ── STRUCTURE K TAB ────────────────────────────────────────────────────
        s("snmf.k_best",        "Best K",            "structure_k", "Structure",
          "numeric", TRUE,
          min = 1, max = 20, step = 1,
          help = "Selected K for all downstream analysis. Review cross-entropy plot first."),
        s("climate.enabled",    "Enable climate",    "structure_k", "Climate",
          "checkbox", TRUE,
          help = "Download WorldClim data and enable GEA / maladaptation analysis"),
        s("climate.predictors", "Climate predictors", "association", "Climate",
          "text",    TRUE,
          help = "Comma-separated bioclim variable names for GEA (review piemaps in Structure K first, e.g. bio_1,bio_12)",
          placeholder = "bio_1,bio_12"),
        # Optional — map
        s("map.climate_extent", "Map extent",        "structure_k", "Map",
          "text",    FALSE,
          help = "Bounding box [min_lon,max_lon,min_lat,max_lat] or 'auto'",
          placeholder = "auto"),
        s("map.gap",            "Auto extent gap",   "structure_k", "Map",
          "numeric", FALSE,
          min = 0, max = 10, step = 0.5,
          help = "Degrees of buffer around sample coordinates (auto extent only)"),
        s("map.resolution",     "Resolution (min)",  "structure_k", "Map",
          "numeric", FALSE,
          min = 0.5, max = 10, step = 0.5,
          help = "WorldClim resolution in arc-minutes (lower = finer but slower)"),
        s("map.zoom_extent",    "Zoom extent",       "structure_k", "Map",
          "text",    FALSE,
          help = "Optional zoom region for piemaps: xmin,xmax,ymin,ymax",
          placeholder = "NULL"),
        # Optional — population stats
        s("pop.calc_stats",       "Calculate pop stats", "structure_k", "Pop Stats",
          "checkbox", FALSE,
          help = "Compute Tajima's D, Pi diversity, AMOVA, and IBD analysis"),
        s("pop.window_size",      "Window size (bp)",   "structure_k", "Pop Stats",
          "numeric",  FALSE,
          min = 1000, step = 1000,
          help = "Sliding window size in base pairs for Tajima's D and Pi diversity"),
        s("pop.custom_trait_file","Custom trait file",  "structure_k", "Pop Stats",
          "text",     FALSE,
          help = "Custom trait file for piemap sizing (NULL = use metadata phenotype columns)",
          placeholder = "NULL"),
        # Optional — piemap
        s("piemap.alpha",       "Pie transparency",  "structure_k", "Piemap",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.05,
          help = "Transparency of pie slices (0 = transparent, 1 = opaque)"),
        s("piemap.show_labels", "Show pop labels",   "structure_k", "Piemap",
          "checkbox", FALSE,
          help = "Overlay population name labels on piemaps"),
        s("piemap.label_size",  "Label size",        "structure_k", "Piemap",
          "numeric",  FALSE,
          min = 4, max = 20, step = 1,
          help = "Font size for population labels"),
        s("piemap.pie_scale",   "Pie size scale",    "structure_k", "Piemap",
          "numeric", FALSE,
          min = 0.1, max = 5, step = 0.1,
          help = "Multiplier for pie chart diameter (1.0 = default)"),
        s("piemap.use_points",  "Use points",        "structure_k", "Piemap",
          "checkbox", FALSE,
          help = "Replace pie charts with single colored points (for large datasets)"),
        # Optional — LD decay
        s("ld_decay.group_by",    "Group by",        "structure_k", "LD Decay",
          "select", FALSE,
          choices = c("site", "cluster"),
          help = "Variable to group samples by for LD decay curves"),
        s("ld_decay.min_samples", "Min samples",     "structure_k", "LD Decay",
          "numeric", FALSE,
          min = 2, step = 1,
          help = "Minimum group size; smaller groups are excluded"),
        s("ld_decay.max_distance","Max dist (kb)",   "structure_k", "LD Decay",
          "numeric", FALSE,
          min = 10, step = 50,
          help = "Maximum pairwise distance in kb for LD calculations"),
        s("ld_decay.scope",       "Scope",           "structure_k", "LD Decay",
          "select", FALSE,
          choices = c("both", "genome_wide", "per_chromosome"),
          help = "Compute LD decay genome-wide, per chromosome, or both"),

        # ── ASSOCIATION TAB ────────────────────────────────────────────────────
        s("association.configs", "Association methods", "association", "Methods",
          "method_table", TRUE,
          help = "Methods to run: EMMAX, LFMM, GLM, MLM, CMLM, ECMLM, SUPER, MLMM, FarmCPU, BLINK"),
        s("association.go_field", "GO field in GFF",  "association", "Methods",
          "text",    FALSE,
          help = "GFF attribute containing GO term IDs. NULL disables enrichment analysis.",
          placeholder = "NULL"),
        # Results processing (combine_gap and region_distance removed — interactive-only in filter bar)
        s("association.promoter_length",  "Promoter length (bp)","association","Results",
          "numeric", FALSE,
          min = 0, step = 1000,
          help = "Upstream distance from gene start counted as promoter region"),
        # Gene annotation
        s("gff.feature",   "GFF feature type", "association", "Gene Annotation",
          "text",    FALSE,
          help = "Feature type to extract from GFF3 for gene overlaps",
          placeholder = "mRNA"),
        s("gff.gene_name", "Gene name field",  "association", "Gene Annotation",
          "text",    FALSE,
          help = "GFF attribute used as gene display name in output tables",
          placeholder = "description"),
        s("gff.biotype",   "Biotype field",    "association", "Gene Annotation",
          "text",    FALSE,
          help = "GFF attribute for gene biotype (used in filtering)",
          placeholder = "biotype"),
        # Enrichment
        s("enrichment.top_terms",        "Top GO terms",      "association", "Enrichment",
          "numeric", FALSE,
          min = 5, max = 100, step = 5,
          help = "Number of enriched GO terms shown in dot/emap/cnetplot"),
        s("enrichment.cnet_label",       "Cnet gene label",   "association", "Enrichment",
          "select",  FALSE,
          choices = c("gene_id", "Name", "description"),
          help = "Which gene attribute to show as labels in the concept network plot"),
        # Region plot (merged from legacy regionplot tab)
        s("regionplot.genes", "Genes to label", "association", "Region Plot",
          "text",    FALSE,
          help = "Which genes to label in regional Manhattan: 'all', 'nameless', or 'gene1|gene2'",
          placeholder = "all"),
        # ── PHENOTYPE ASSOCIATION TAB ──────────────────────────────────────────
        s("phenotype_association.configs",          "GWAS methods",        "phenotype", "Methods",
          "method_table", TRUE,
          help = "Association methods for phenotype GWAS. Same options as GEA methods."),
        s("phenotype_association.missing_strategy", "Missing data strategy","phenotype", "Methods",
          "select",       TRUE,
          choices = c("DROP", "MEAN", "MEDIAN"),
          help = "How to handle samples missing phenotype data. DROP removes them per trait."),
        # Inherited overrides (optional — blank = inherit from association)
        s("phenotype_association.combine_method",   "Combine method",      "phenotype", "Results (overrides GEA)",
          "select",  FALSE,
          choices = c("", "Sum", "Overlap", "EMMAX", "LFMM", "GLM", "MLM", "FarmCPU", "BLINK"),
          help = "Override GEA combine method for GWAS only. Leave blank to inherit."),
        s("phenotype_association.region_distance",  "Region extent (bp)",  "phenotype", "Results (overrides GEA)",
          "text",    FALSE,
          help = "Override GEA region_distance for GWAS. Use 'auto' or blank to inherit.",
          placeholder = "(inherit from GEA)"),
        s("phenotype_association.promoter_length",  "Promoter length (bp)","phenotype", "Results (overrides GEA)",
          "numeric", FALSE,
          min = 0, step = 1000,
          help = "Override GEA promoter_length for GWAS. Blank = inherit."),
        s("regionplot.genes", "Genes to label", "phenotype", "Region Plot",
          "text",    FALSE,
          help = "Which genes to label in regional Manhattan: 'all', 'nameless', or 'gene1|gene2'",
          placeholder = "all"),

        # ── OVERLAPPING REGIONS TAB ────────────────────────────────────────────
        s("overlap.region_distance",      "Region extent (bp)",  "overlapping", "Region Merging",
          "text",    FALSE,
          help = "Region distance for combined GEA+GWAS analysis. Defaults to max(GEA, GWAS).",
          placeholder = "(max of GEA and GWAS)"),
        s("overlap.promoter_length",      "Promoter length (bp)","overlapping", "Region Merging",
          "numeric", FALSE,
          min = 0, step = 1000,
          help = "Promoter length for combined overlap analysis"),
        s("overlap.pairwise.window_size", "Pairwise window (bp)","overlapping", "Pairwise Analysis",
          "numeric", FALSE,
          min = 0, step = 100000,
          help = "Window size for scoring GEA vs GWAS pairwise region overlap"),
        s("overlap.pairwise.min_snps",    "Min pairwise SNPs",   "overlapping", "Pairwise Analysis",
          "numeric", FALSE,
          min = 1, step = 1,
          help = "Minimum shared significant SNPs to call a pairwise overlap"),
        s("regionplot.genes", "Genes to label", "overlapping", "Region Plot",
          "text",    FALSE,
          help = "Which genes to label in regional Manhattan: 'all', 'nameless', or 'gene1|gene2'",
          placeholder = "all"),

        # ── HAPLOTYPE TAB ──────────────────────────────────────────────────────
        s("haplotype.scan.regions_source", "Regions source",    "haplotype", "Scan Setup",
          "select",  TRUE,
          choices = c("association", "association_phenotypes", "overlapping", "all", "custom"),
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
        s("future.ssp",    "Climate scenario (SSP)","maladaptation","Future Climate",
          "select",  TRUE,
          choices = c("126", "245", "370", "585"),
          help = "Shared Socioeconomic Pathway (126 = low emissions, 585 = high emissions)"),
        s("future.year",   "Projection period",    "maladaptation","Future Climate",
          "select",  TRUE,
          choices = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"),
          help = "Decade range for future climate projections"),
        s("future.models", "Climate models",       "maladaptation","Future Climate",
          "text",    FALSE,
          help = "Comma-separated CMIP6 model names to average. Blank = use all available.",
          placeholder = "(all models)"),
        s("gradient_forest.combine_method",    "SNP selection",      "maladaptation","Gradient Forest",
          "select",  TRUE,
          choices = c("All", "Overlap", "MethodOverlap", "EMMAX", "LFMM", "GLM", "MLM", "FarmCPU", "BLINK"),
          help = "Which sig SNPs feed into GF model. Explore strategies in the Association tab first, then set here."),
        s("gradient_forest.combine_gap",       "SNP combine gap (bp)","maladaptation","Gradient Forest",
          "numeric", TRUE,
          min = 0, step = 10000,
          help = "Max distance for cross-method overlap (only used for Overlap/MethodOverlap strategies)"),
        s("gradient_forest.ntree",             "Number of trees",    "maladaptation","Gradient Forest",
          "numeric", FALSE,
          min = 100, step = 100,
          help = "Number of trees per Gradient Forest run (more = more stable, slower)"),
        s("gradient_forest.cor_threshold",     "R\u00b2 threshold",  "maladaptation","Gradient Forest",
          "numeric", FALSE,
          min = 0, max = 1, step = 0.05,
          help = "Minimum R\u00b2 to include a climate predictor in the model"),
        s("gradient_forest.spatial_correction","Spatial correction", "maladaptation","Gradient Forest",
          "select",  FALSE,
          choices = c("with", "without"),
          help = "Include spatial eigenvectors (PCNM) as covariates to account for IBD"),
        s("gradient_forest.run_label",         "Run label",          "maladaptation","Gradient Forest",
          "text",    FALSE,
          help = "Unique suffix to distinguish multiple GF runs in the same project",
          placeholder = ""),
        s("gradient_forest.random_model",      "Random model",       "maladaptation","Gradient Forest",
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
