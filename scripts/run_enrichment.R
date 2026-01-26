#!/usr/bin/env Rscript
# Run GO enrichment analysis using clusterProfiler::enricher()
# Uses genes from collapsed gene table and GO terms from GFF

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(clusterProfiler)

args = commandArgs(trailingOnly=TRUE)
################
GENES_COLLAPSED = args[1]  # Genes_per_region_collapsed.tsv
GFF = args[2]
GO_FIELD = args[3]         # Field name containing GO terms in GFF (e.g., "ontology")
GFF_FEATURE = args[4]      # Feature type (e.g., "mRNA")
TABLES_DIR = args[5]       # Output directory for enrichment tables
OUTPUT = args[6]           # Main combined enrichment output
################

message('INFO: Running GO enrichment analysis')
message(paste0('INFO: GO field: ', GO_FIELD))
message(paste0('INFO: GFF feature: ', GFF_FEATURE))

# Check if GO_FIELD is NULL or empty
if (is.null(GO_FIELD) || GO_FIELD == 'NULL' || GO_FIELD == '') {
    message('WARNING: GO_FIELD not specified, skipping enrichment')
    empty_dt <- data.table(
        GO_id = character(),
        description = character(),
        gene_ratio = character(),
        bg_ratio = character(),
        pvalue = numeric(),
        p_adjust = numeric(),
        qvalue = numeric(),
        gene_ids = character(),
        gene_count = integer()
    )
    empty_dt %>% fwrite(OUTPUT, sep = '\t')
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
    # Try to load GO.db
    if (requireNamespace("GO.db", quietly = TRUE)) {
        message('INFO: Using GO.db for term descriptions')
        library(GO.db)

        # Get term descriptions
        go_terms <- AnnotationDbi::select(GO.db,
                                          keys = go_ids,
                                          columns = c("GOID", "TERM"),
                                          keytype = "GOID")

        # Create mapping
        term2name <- data.frame(
            term = go_terms$GOID,
            name = go_terms$TERM,
            stringsAsFactors = FALSE
        )

        # Fill in missing descriptions with GO ID
        term2name$name[is.na(term2name$name)] <- term2name$term[is.na(term2name$name)]

        return(term2name)
    } else {
        message('WARNING: GO.db not available, descriptions will be GO IDs')
        return(NULL)
    }
}

#################### Main

# Load genes from regions
genes <- fread(GENES_COLLAPSED)

if (nrow(genes) == 0) {
    message('WARNING: No genes to process')
    empty_dt <- data.table(
        GO_id = character(),
        description = character(),
        gene_ratio = character(),
        bg_ratio = character(),
        pvalue = numeric(),
        p_adjust = numeric(),
        qvalue = numeric(),
        gene_ids = character(),
        gene_count = integer()
    )
    empty_dt %>% fwrite(OUTPUT, sep = '\t')
    quit(save = "no", status = 0)
}

message(paste0('INFO: Processing ', nrow(genes), ' genes'))

# Get gene IDs of interest
genes_of_interest <- unique(genes$gene_id)
message(paste0('INFO: Genes of interest: ', length(genes_of_interest)))

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
    empty_dt <- data.table(
        GO_id = character(),
        description = character(),
        gene_ratio = character(),
        bg_ratio = character(),
        pvalue = numeric(),
        p_adjust = numeric(),
        qvalue = numeric(),
        gene_ids = character(),
        gene_count = integer()
    )
    empty_dt %>% fwrite(OUTPUT, sep = '\t')
    quit(save = "no", status = 0)
}

message(paste0('INFO: Found GO annotations for ', n_distinct(gff_data$gene_id), ' genes'))

# Create TERM2GENE mapping (GO term -> gene)
# GO terms may be comma-separated
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

if (!is.null(term2name)) {
    message(paste0('INFO: Got descriptions for ', sum(!is.na(term2name$name) & term2name$name != term2name$term), ' GO terms'))
}

# Filter genes of interest to those with GO annotations
genes_with_go <- intersect(genes_of_interest, unique(term2gene$gene))
message(paste0('INFO: Genes of interest with GO annotations: ', length(genes_with_go)))

if (length(genes_with_go) == 0) {
    message('WARNING: No genes of interest have GO annotations')
    empty_dt <- data.table(
        GO_id = character(),
        description = character(),
        gene_ratio = character(),
        bg_ratio = character(),
        pvalue = numeric(),
        p_adjust = numeric(),
        qvalue = numeric(),
        gene_ids = character(),
        gene_count = integer()
    )
    empty_dt %>% fwrite(OUTPUT, sep = '\t')
    quit(save = "no", status = 0)
}

# Run enrichment analysis
message('INFO: Running enrichment analysis')

tryCatch({
    enrich_result <- enricher(
        gene = genes_with_go,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,  # Now includes GO term descriptions
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        pAdjustMethod = "BH",
        minGSSize = 2,
        maxGSSize = 500
    )

    if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
        message('WARNING: No significant enrichment found')
        empty_dt <- data.table(
            GO_id = character(),
            description = character(),
            gene_ratio = character(),
            bg_ratio = character(),
            pvalue = numeric(),
            p_adjust = numeric(),
            qvalue = numeric(),
            gene_ids = character(),
            gene_count = integer()
        )
        empty_dt %>% fwrite(OUTPUT, sep = '\t')
    } else {
        # Extract results
        result_dt <- enrich_result@result %>%
            as.data.table() %>%
            dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count) %>%
            setNames(c('GO_id', 'description', 'gene_ratio', 'bg_ratio', 'pvalue', 'p_adjust', 'qvalue', 'gene_ids', 'gene_count'))

        # If description is still NA or same as GO_id, try to look it up again
        if (!is.null(term2name)) {
            for (i in seq_len(nrow(result_dt))) {
                if (is.na(result_dt$description[i]) || result_dt$description[i] == result_dt$GO_id[i]) {
                    match_idx <- which(term2name$term == result_dt$GO_id[i])
                    if (length(match_idx) > 0 && !is.na(term2name$name[match_idx[1]])) {
                        result_dt$description[i] <- term2name$name[match_idx[1]]
                    }
                }
            }
        }

        message(paste0('INFO: Found ', nrow(result_dt), ' enriched terms'))

        # Save combined result
        result_dt %>% fwrite(OUTPUT, sep = '\t')

        message(paste0('INFO: Saved enrichment results to ', OUTPUT))
    }
}, error = function(e) {
    message(paste0('ERROR in enrichment: ', e$message))
    empty_dt <- data.table(
        GO_id = character(),
        description = character(),
        gene_ratio = character(),
        bg_ratio = character(),
        pvalue = numeric(),
        p_adjust = numeric(),
        qvalue = numeric(),
        gene_ids = character(),
        gene_count = integer()
    )
    empty_dt %>% fwrite(OUTPUT, sep = '\t')
})

message('INFO: Complete')
