#!/usr/bin/env Rscript
# Find genes overlapping regions
# Regions already have extended boundaries from create_regions.R — used as-is here
# Outputs: redundant (one row per gene-region) and collapsed (one row per gene)

suppressPackageStartupMessages({
    library(tidyr)
    library(data.table)
    library(dplyr)
    library(stringr)
})

source("/pipeline/scripts/R/lib/gff_parsing.R")
source("/pipeline/scripts/R/lib/genes_in_regions.R")

args = commandArgs(trailingOnly = TRUE)
################
GFF              = args[1]
REGIONS          = args[2]
GFF_FEATURE      = args[3]
PROMOTER_LENGTH  = as.numeric(args[4])
ALL_SNPS         = args[5]  # vcfsnp file
CPU              = as.numeric(args[6])
OUTPUT_GENES     = args[7]
OUTPUT_COLLAPSED = args[8]
################

message('INFO: Finding genes overlapping regions')
message(paste0('INFO: GFF feature: ', GFF_FEATURE))

regions <- fread(REGIONS, colClasses = c(chr = "character"))

empty_gene_tables <- function(reason) {
    message(paste0('WARNING: ', reason, ' — writing empty gene tables'))
    empty <- data.table(region_id = character(), gene_id = character(),
                        chr = character(), gene_start = integer(), gene_end = integer())
    fwrite(empty, OUTPUT_GENES,     sep = '\t')
    fwrite(empty, OUTPUT_COLLAPSED, sep = '\t')
    quit(save = "no", status = 0)
}

if (nrow(regions) == 0) {
    empty_gene_tables('No regions to process')
}

# Input.gff is optional (common.smk guards every use with `if GFF:`), and normalize_gff
# writes a zero-byte file when it is unset — the normal case for simulated datasets, which
# ship no annotation. Without this guard read_gff() runs `grep -v '#'` on the empty file,
# grep exits 1 for "no lines matched", and fread reports it as a misleading
# "disk is full in the temporary directory" error. Gene annotation is simply inert with no
# GFF, so emit the same empty tables as the no-regions case.
if (!file.exists(GFF) || file.size(GFF) == 0) {
    empty_gene_tables(paste0('No GFF annotation supplied (', GFF, ' is empty)'))
}

allsnps <- fread(ALL_SNPS, header = FALSE) %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    as.data.table()

gff_dt <- read_gff(GFF, GFF_FEATURE)

res <- find_genes_for_regions(regions, gff_dt, gff_path = GFF,
                               promoter_length = PROMOTER_LENGTH,
                               allsnps_dt = allsnps)

message(paste0('INFO: Saving redundant table to ', OUTPUT_GENES))
fwrite(res$genes_per_region, OUTPUT_GENES, sep = '\t', quote = FALSE)

message(paste0('INFO: Saving collapsed table to ', OUTPUT_COLLAPSED))
fwrite(res$genes_collapsed, OUTPUT_COLLAPSED, sep = '\t', quote = FALSE)

message('INFO: Complete')
