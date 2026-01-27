#!/usr/bin/env Rscript
# Find genes around significant regions
# Outputs: redundant (one row per gene-region) and collapsed (one row per gene)

library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
################
GFF = args[1]
REGIONS = args[2]
GFF_FEATURE = args[3]
DISTANCE = args[4] %>% as.numeric
PROMOTER_LENGTH = args[5] %>% as.numeric
ALL_SNPS = args[6]  # vcfsnp file
CPU = args[7] %>% as.numeric
OUTPUT_GENES = args[8]
OUTPUT_COLLAPSED = args[9]
################

message('INFO: Finding genes around significant regions')
message(paste0('INFO: Distance: ', DISTANCE))
message(paste0('INFO: GFF feature: ', GFF_FEATURE))

######################## Functions

# Harmonize chromosome names between two character vectors
# Strips common prefixes (chr, Chr, CHR) if one side has them and the other doesn't
# Handles formats: chr1, Chr1, CHR1, 1, chr1H, chr1H_1
harmonize_chr_names <- function(chr_a, chr_b) {
    has_prefix_a <- any(grepl("^[Cc][Hh][Rr]", chr_a))
    has_prefix_b <- any(grepl("^[Cc][Hh][Rr]", chr_b))

    if (has_prefix_a && !has_prefix_b) {
        # GFF has prefix, regions don't — strip prefix from GFF side
        message('INFO: Harmonizing chr names: stripping chr prefix from source A to match source B')
        return(list(a = sub("^[Cc][Hh][Rr]", "", chr_a), b = chr_b))
    } else if (!has_prefix_a && has_prefix_b) {
        # Regions have prefix, GFF doesn't — strip prefix from regions side
        message('INFO: Harmonizing chr names: stripping chr prefix from source B to match source A')
        return(list(a = chr_a, b = sub("^[Cc][Hh][Rr]", "", chr_b)))
    }
    # Both have same convention — no change needed
    return(list(a = chr_a, b = chr_b))
}

# Extract gene ID from Parent or ID field (only these are assumed standard)
extract_gene_id <- function(attr) {
    id <- str_extract(attr, "(?<=Parent=)[^;]+")
    id <- ifelse(is.na(id), str_extract(attr, "(?<=ID=)[^;]+"), id)
    id <- str_remove(id, "\\.[0-9]+(_exon_[0-9]+)?$")
    return(id)
}

# Clean attribute values by removing key= prefix
clean_attr_value <- function(x) {
    # Remove patterns like "ID=", "Parent=", "description=", etc.
    str_remove(x, "^[^=]+=")
}

#################### Main

# Load regions
regions <- fread(REGIONS)
if ('chr' %in% colnames(regions)) regions$chr <- as.character(regions$chr)

if (nrow(regions) == 0) {
    message('WARNING: No regions to process')
    empty_genes <- data.table(
        region_id = character(), gene_id = character(), chr = character(),
        gene_start = integer(), gene_end = integer()
    )
    empty_genes %>% fwrite(OUTPUT_GENES, sep = '\t')
    empty_genes %>% fwrite(OUTPUT_COLLAPSED, sep = '\t')
    quit(save = "no", status = 0)
}

message(paste0('INFO: Processing ', nrow(regions), ' regions'))

# Load all SNPs for exon/promoter counting
message('INFO: Loading all SNPs')
allsnps <- fread(ALL_SNPS) %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chr = as.character(chr))

# Load GFF
message('INFO: Loading GFF')
genes_gff_raw <- fread(cmd = paste("grep -v '#'", GFF), header = F) %>%
    dplyr::filter(V3 == !!GFF_FEATURE) %>%
    dplyr::select(V1, V4, V5, V9) %>%
    setNames(c('chr', 'start', 'end', 'description')) %>%
    dplyr::mutate(chr = as.character(chr))

# Extract gene ID
genes_gff_raw <- genes_gff_raw %>%
    dplyr::mutate(gene_id = description %>% extract_gene_id)

# Parse description fields dynamically
# Extract all keys from first description
sample_desc <- genes_gff_raw$description[1]
fields <- sample_desc %>%
    str_extract_all("(?<=^|;)[^=;]+(?==)", simplify = TRUE) %>%
    as.character() %>%
    unique() %>%
    .[. != ""]

message(paste0('INFO: Found GFF fields: ', paste(fields, collapse = ', ')))

# Split description into columns
genes_gff <- genes_gff_raw %>%
    tidyr::separate(col = "description", into = fields, sep = ";", fill = "right", extra = "drop")

# Clean attribute values (remove "key=" prefixes)
for (field in fields) {
    if (field %in% colnames(genes_gff)) {
        genes_gff[[field]] <- clean_attr_value(genes_gff[[field]])
        genes_gff[[field]] <- na_if(genes_gff[[field]], "NA")
    }
}

genes_gff <- genes_gff %>%
    dplyr::mutate(start = as.integer(start), end = as.integer(end))

# Harmonize chromosome names between regions and GFF
# (e.g., VCF-derived regions may have "1" while GFF has "chr1")
chr_harmonized <- harmonize_chr_names(genes_gff$chr, regions$chr)
genes_gff$chr <- chr_harmonized$a
regions$chr <- chr_harmonized$b
# Also harmonize allsnps
allsnps$chr <- harmonize_chr_names(allsnps$chr, regions$chr)$a

# Create GRanges for genes
genes_gr <- genes_gff %>%
    dplyr::select(chr, start, end) %>%
    GRanges()

# Extend genes by DISTANCE/2 on each side
genes_extended_gr <- GRanges(
    seqnames = seqnames(genes_gr),
    ranges = IRanges(
        start = pmax(1, start(genes_gr) - DISTANCE / 2),
        end = end(genes_gr) + DISTANCE / 2
    )
)

# Create GRanges for regions
regions_gr <- regions %>%
    dplyr::select(chr, start, end, region_id) %>%
    GRanges()

# Find overlaps between regions and extended genes
message('INFO: Finding gene-region overlaps')
overlaps <- findOverlaps(regions_gr, genes_extended_gr)

if (length(overlaps) == 0) {
    message('WARNING: No genes found around regions')
    empty_genes <- data.table(
        region_id = character(), gene_id = character(), chr = character(),
        gene_start = integer(), gene_end = integer()
    )
    empty_genes %>% fwrite(OUTPUT_GENES, sep = '\t')
    empty_genes %>% fwrite(OUTPUT_COLLAPSED, sep = '\t')
    quit(save = "no", status = 0)
}

# Build redundant table (one row per region-gene pair)
genes_per_region <- data.table(
    region_id = regions$region_id[queryHits(overlaps)],
    genes_gff[subjectHits(overlaps), ]
)

# Rename columns for clarity
colnames(genes_per_region)[colnames(genes_per_region) == 'start'] <- 'gene_start'
colnames(genes_per_region)[colnames(genes_per_region) == 'end'] <- 'gene_end'

message(paste0('INFO: Found ', nrow(genes_per_region), ' gene-region pairs'))

# Count SNPs in exons and promoters
message('INFO: Counting exon/promoter SNPs')

# Load exon data for SNP counting (apply same chr harmonization as genes_gff)
exon_gff <- fread(cmd = paste("grep -v '#'", GFF), header = F) %>%
    dplyr::filter(V3 == 'exon') %>%
    dplyr::select(V1, V4, V5, V9) %>%
    setNames(c('chr', 'start', 'end', 'description')) %>%
    dplyr::mutate(chr = as.character(chr),
                  gene_id = description %>% extract_gene_id) %>%
    dplyr::select(-description)
exon_gff$chr <- harmonize_chr_names(exon_gff$chr, regions$chr)$a

# Get unique gene IDs
unique_genes <- unique(genes_per_region$gene_id)

# Count SNPs in exons for each gene
exon_snp_counts <- if (nrow(exon_gff) > 0) {
    exons_gr <- exon_gff %>%
        dplyr::filter(gene_id %in% unique_genes) %>%
        GRanges()

    snps_gr <- allsnps %>%
        dplyr::mutate(start = pos, end = pos) %>%
        dplyr::select(chr, start, end) %>%
        GRanges()

    if (length(exons_gr) > 0 && length(snps_gr) > 0) {
        exon_overlaps <- findOverlapPairs(snps_gr, exons_gr)
        if (length(exon_overlaps) > 0) {
            exon_overlaps %>%
                as.data.table %>%
                dplyr::select(first.seqnames, first.start, second.gene_id) %>%
                setNames(c('chr', 'pos', 'gene_id')) %>%
                dplyr::mutate(SNPID = paste0(chr, ':', pos)) %>%
                group_by(gene_id) %>%
                summarise(exon_snps = paste(unique(SNPID), collapse = ','), .groups = 'drop')
        } else {
            data.table(gene_id = character(), exon_snps = character())
        }
    } else {
        data.table(gene_id = character(), exon_snps = character())
    }
} else {
    data.table(gene_id = character(), exon_snps = character())
}

# Count SNPs in promoters for each gene
promoter_gff <- genes_gff %>%
    dplyr::filter(gene_id %in% unique_genes) %>%
    dplyr::mutate(
        promoter_start = pmax(1, start - PROMOTER_LENGTH),
        promoter_end = start
    ) %>%
    dplyr::select(chr, promoter_start, promoter_end, gene_id)

promoter_snp_counts <- if (nrow(promoter_gff) > 0) {
    promoters_gr <- promoter_gff %>%
        dplyr::rename(start = promoter_start, end = promoter_end) %>%
        GRanges()

    snps_gr <- allsnps %>%
        dplyr::mutate(start = pos, end = pos) %>%
        dplyr::select(chr, start, end) %>%
        GRanges()

    if (length(promoters_gr) > 0 && length(snps_gr) > 0) {
        prom_overlaps <- findOverlapPairs(snps_gr, promoters_gr)
        if (length(prom_overlaps) > 0) {
            prom_overlaps %>%
                as.data.table %>%
                dplyr::select(first.seqnames, first.start, second.gene_id) %>%
                setNames(c('chr', 'pos', 'gene_id')) %>%
                dplyr::mutate(SNPID = paste0(chr, ':', pos)) %>%
                group_by(gene_id) %>%
                summarise(promoter_snps = paste(unique(SNPID), collapse = ','), .groups = 'drop')
        } else {
            data.table(gene_id = character(), promoter_snps = character())
        }
    } else {
        data.table(gene_id = character(), promoter_snps = character())
    }
} else {
    data.table(gene_id = character(), promoter_snps = character())
}

# Join SNP counts to genes
genes_per_region <- genes_per_region %>%
    left_join(exon_snp_counts, by = 'gene_id') %>%
    left_join(promoter_snp_counts, by = 'gene_id')

# Replace NA with empty string for SNP columns
genes_per_region$exon_snps[is.na(genes_per_region$exon_snps)] <- ''
genes_per_region$promoter_snps[is.na(genes_per_region$promoter_snps)] <- ''

# Reorder columns: region_id first, then gene_id, then rest
col_order <- c('region_id', 'gene_id', 'chr', 'gene_start', 'gene_end',
               setdiff(colnames(genes_per_region),
                       c('region_id', 'gene_id', 'chr', 'gene_start', 'gene_end')))
genes_per_region <- genes_per_region %>% dplyr::select(all_of(col_order))

# Save redundant table
message(paste0('INFO: Saving redundant table to ', OUTPUT_GENES))
genes_per_region %>% fwrite(OUTPUT_GENES, sep = '\t')

# Create collapsed table (one row per gene)
message('INFO: Creating collapsed table')

# Define columns to collapse (all except key identifiers)
key_cols <- c('gene_id', 'chr', 'gene_start', 'gene_end')
collapse_cols <- setdiff(colnames(genes_per_region), c(key_cols, 'region_id', 'exon_snps', 'promoter_snps'))

genes_collapsed <- genes_per_region %>%
    group_by(gene_id, chr, gene_start, gene_end) %>%
    summarise(
        region_id = paste(unique(region_id), collapse = ','),
        exon_snps = paste(unique(exon_snps[exon_snps != '']), collapse = ','),
        promoter_snps = paste(unique(promoter_snps[promoter_snps != '']), collapse = ','),
        across(all_of(collapse_cols), ~ paste(unique(.x[!is.na(.x)]), collapse = ',')),
        .groups = 'drop'
    ) %>%
    dplyr::arrange(chr, gene_start)

# Save collapsed table
message(paste0('INFO: Saving collapsed table to ', OUTPUT_COLLAPSED))
genes_collapsed %>% fwrite(OUTPUT_COLLAPSED, sep = '\t')

message(paste0('INFO: Found ', n_distinct(genes_per_region$gene_id), ' unique genes across ', nrow(regions), ' regions'))
message('INFO: Complete')
