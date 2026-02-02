#!/usr/bin/env Rscript
# Run per-region GO enrichment analysis using clusterProfiler::enricher()
# Processes each region separately and creates per-trait summary tables

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(clusterProfiler)
library(qs)

args = commandArgs(trailingOnly=TRUE)
################
GENES_PER_REGION = args[1]  # Genes_per_region.tsv (redundant: one row per region-gene pair)
GFF = args[2]
GO_FIELD = args[3]          # Field name containing GO terms in GFF (e.g., "ontology")
GFF_FEATURE = args[4]       # Feature type (e.g., "mRNA")
TABLES_DIR = args[5]        # Base directory: tables/association/enrichment/
INTERMEDIATE_DIR = args[6]  # Base directory: intermediate/enrichment/
################

message('INFO: Running per-region GO enrichment analysis')
message(paste0('INFO: GO field: ', GO_FIELD))
message(paste0('INFO: GFF feature: ', GFF_FEATURE))
message(paste0('INFO: TSV output directory: ', TABLES_DIR))
message(paste0('INFO: QS output directory: ', INTERMEDIATE_DIR))

# Check if GO_FIELD is NULL or empty
if (is.null(GO_FIELD) || GO_FIELD == 'NULL' || GO_FIELD == '') {
    message('WARNING: GO_FIELD not specified, skipping enrichment')
    quit(save = "no", status = 0)
}

######################## Functions

# Extract gene ID from Parent or ID field
extract_gene_id <- function(attr) {
    id <- str_extract(attr, "(?<=Parent=)[^;]+")
    id <- ifelse(is.na(id), str_extract(attr, "(?<=ID=)[^;]+"), id)
    id <- str_remove(id, "\\.[0-9]+(_exon_[0-9]+)?$")
    return(id)
}

# Extract specific field value from GFF description
extract_field <- function(description, field_name) {
    pattern <- paste0("(?<=", field_name, "=)[^;]+")
    str_extract(description, pattern)
}

# Get GO term descriptions from GO.db (if available)
get_go_descriptions <- function(go_ids) {
    if (requireNamespace("GO.db", quietly = TRUE)) {
        message('INFO: Using GO.db for term descriptions')
        library(GO.db)

        go_terms <- AnnotationDbi::select(GO.db,
                                          keys = go_ids,
                                          columns = c("GOID", "TERM"),
                                          keytype = "GOID")

        term2name <- data.frame(
            term = go_terms$GOID,
            name = go_terms$TERM,
            stringsAsFactors = FALSE
        )

        term2name$name[is.na(term2name$name)] <- term2name$term[is.na(term2name$name)]
        return(term2name)
    } else {
        message('WARNING: GO.db not available, descriptions will be GO IDs')
        return(NULL)
    }
}

# Run enrichment for a single region
run_enrichment_for_region <- function(region_id, trait, genes_in_region, term2gene, term2name, all_genes_with_go) {
    message(paste0('INFO: Running enrichment for region: ', region_id))

    # Filter to genes with GO annotations
    genes_with_go <- intersect(genes_in_region, all_genes_with_go)

    if (length(genes_with_go) < 2) {
        message(paste0('WARNING: Region ', region_id, ' has <2 genes with GO annotations (', length(genes_with_go), ')'))
        return(NULL)
    }

    message(paste0('INFO: Region ', region_id, ': ', length(genes_in_region), ' total genes, ', length(genes_with_go), ' with GO'))

    tryCatch({
        enrich_result <- enricher(
            gene = genes_with_go,
            TERM2GENE = term2gene,
            TERM2NAME = term2name,
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            pAdjustMethod = "BH",
            minGSSize = 2,
            maxGSSize = 500
        )

        if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
            message(paste0('WARNING: No significant enrichment for region ', region_id))
            return(NULL)
        }

        # Fix qvalue if NA (common with small datasets where qvalue package fails)
        if (nrow(enrich_result@result) > 0) {
            if (all(is.na(enrich_result@result$qvalue))) {
                message('INFO: qvalue is NA, using p.adjust as qvalue')
                enrich_result@result$qvalue <- enrich_result@result$p.adjust
            }
        }

        message(paste0('INFO: Region ', region_id, ': Found ', nrow(enrich_result@result), ' enriched terms'))

        # Return the enrichResult object directly (will be saved with qs)
        return(list(
            region_id = region_id,
            trait = trait,
            enrich_obj = enrich_result
        ))

    }, error = function(e) {
        message(paste0('ERROR in enrichment for region ', region_id, ': ', e$message))
        return(NULL)
    })
}

#################### Main

# Load genes from regions (redundant table: one row per region-gene pair)
genes <- fread(GENES_PER_REGION)

if (nrow(genes) == 0) {
    message('WARNING: No genes to process')
    quit(save = "no", status = 0)
}

message(paste0('INFO: Loaded ', nrow(genes), ' region-gene pairs'))

# Extract trait from region_id (format: chr:start-end_trait)
# region_id examples: "1:10238549-10683918_bio_12", "2:25222568-25222568_bio_1"
# Extract everything after the first underscore (which comes after genomic coordinates)
genes <- genes %>%
    dplyr::mutate(trait = sub("^[^_]+_", "", region_id))

# Check if trait extraction worked
if (all(is.na(genes$trait))) {
    message('ERROR: Could not extract trait from region_id. Expected format: chr:start-end_TRAIT')
    message('Example region_id from data: ', genes$region_id[1])
    quit(save = "no", status = 1)
}

message(paste0('INFO: Found traits: ', paste(unique(genes$trait), collapse = ', ')))
message(paste0('INFO: Found regions: ', n_distinct(genes$region_id)))

# Load GFF and extract GO terms for all genes (background)
message('INFO: Loading GFF for background genes')
gff_raw <- fread(cmd = paste("grep -v '#'", GFF), header = F) %>%
    dplyr::filter(V3 == !!GFF_FEATURE) %>%
    dplyr::select(V9) %>%
    setNames('description')

# Extract gene IDs and GO terms
gff_data <- gff_raw %>%
    dplyr::mutate(
        gene_id = extract_gene_id(description),
        go_terms = extract_field(description, GO_FIELD)
    ) %>%
    dplyr::filter(!is.na(go_terms) & go_terms != '' & go_terms != 'NA') %>%
    dplyr::select(gene_id, go_terms)

if (nrow(gff_data) == 0) {
    message(paste0('WARNING: No GO terms found in field "', GO_FIELD, '"'))
    quit(save = "no", status = 0)
}

message(paste0('INFO: Found GO annotations for ', n_distinct(gff_data$gene_id), ' background genes'))

# Create TERM2GENE mapping (GO term -> gene)
term2gene <- gff_data %>%
    tidyr::separate_rows(go_terms, sep = ',') %>%
    dplyr::rename(term = go_terms, gene = gene_id) %>%
    dplyr::filter(term != '' & !is.na(term)) %>%
    dplyr::select(term, gene) %>%
    as.data.frame()

message(paste0('INFO: Created TERM2GENE with ', nrow(term2gene), ' term-gene pairs'))
message(paste0('INFO: Unique GO terms: ', n_distinct(term2gene$term)))

# Get GO term descriptions
unique_go_ids <- unique(term2gene$term)
term2name <- get_go_descriptions(unique_go_ids)

# Get all genes with GO annotations (for filtering)
all_genes_with_go <- unique(term2gene$gene)

# Create output directories
traits <- unique(genes$trait)
for (trait in traits) {
    # Create tables directory for TSV files
    trait_dir_tables <- file.path(TABLES_DIR, trait)
    if (!dir.exists(trait_dir_tables)) {
        dir.create(trait_dir_tables, recursive = TRUE)
        message(paste0('INFO: Created directory: ', trait_dir_tables))
    }

    # Create intermediate directory for .qs files
    trait_dir_inter <- file.path(INTERMEDIATE_DIR, trait)
    if (!dir.exists(trait_dir_inter)) {
        dir.create(trait_dir_inter, recursive = TRUE)
        message(paste0('INFO: Created directory: ', trait_dir_inter))
    }
}

# Group genes by region
genes_by_region <- genes %>%
    group_by(region_id, trait) %>%
    summarise(gene_ids = list(unique(gene_id)), .groups = 'drop')

message(paste0('INFO: Processing ', nrow(genes_by_region), ' unique regions'))

# Process each region separately
all_results <- list()
result_counter <- 0

for (i in seq_len(nrow(genes_by_region))) {
    region_id <- genes_by_region$region_id[i]
    trait <- genes_by_region$trait[i]
    genes_in_region <- genes_by_region$gene_ids[[i]]

    # Skip if no genes
    if (length(genes_in_region) == 0) {
        message(paste0('WARNING: Region ', region_id, ' has no genes'))
        next
    }

    # Run enrichment
    result <- run_enrichment_for_region(
        region_id = region_id,
        trait = trait,
        genes_in_region = genes_in_region,
        term2gene = term2gene,
        term2name = term2name,
        all_genes_with_go = all_genes_with_go
    )

    if (!is.null(result)) {
        result_counter <- result_counter + 1

        region_id <- result$region_id
        trait <- result$trait
        enrich_obj <- result$enrich_obj

        # Save enrichResult object as .qs file to intermediate directory
        qs_file <- file.path(INTERMEDIATE_DIR, trait, paste0('Region_', region_id, '_enrichment.qs'))
        qs::qsave(enrich_obj, qs_file)
        message(paste0('INFO: Saved enrichment object to ', qs_file))

        # Also save TSV summary for easy viewing
        result_dt <- enrich_obj@result %>%
            as.data.table() %>%
            dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
            setNames(c('GO_id', 'description', 'gene_ratio', 'bg_ratio', 'pvalue', 'p_adjust', 'qvalue', 'gene_ids', 'gene_count'))
        result_dt$region_id <- region_id
        result_dt$trait <- trait
        result_dt <- result_dt %>%
            dplyr::select(region_id, trait, GO_id, description, gene_ratio, bg_ratio, pvalue, p_adjust, qvalue, gene_ids, gene_count)

        tsv_file <- file.path(TABLES_DIR, trait, paste0('Region_', region_id, '_enrichment.tsv'))
        result_dt %>% fwrite(tsv_file, sep = '\t')
        message(paste0('INFO: Saved TSV summary to ', tsv_file))

        # Store for summary
        all_results[[result_counter]] <- result_dt
    }
}

# Create per-trait summary files
if (length(all_results) > 0) {
    combined_results <- rbindlist(all_results)

    # Create per-trait summaries
    for (trait in traits) {
        trait_results <- combined_results %>% dplyr::filter(trait == !!trait)

        if (nrow(trait_results) > 0) {
            summary_file <- file.path(TABLES_DIR, trait, paste0('Enrichment_', trait, '_summary.tsv'))
            trait_results %>% fwrite(summary_file, sep = '\t')
            message(paste0('INFO: Created trait summary: ', summary_file))
            message(paste0('INFO: Trait ', trait, ': ', n_distinct(trait_results$region_id), ' regions, ', nrow(trait_results), ' total enriched terms'))
        }
    }

    message(paste0('INFO: Enrichment complete: ', result_counter, ' regions with significant terms'))
} else {
    message('WARNING: No regions had significant enrichment')
}

message('INFO: Complete')
