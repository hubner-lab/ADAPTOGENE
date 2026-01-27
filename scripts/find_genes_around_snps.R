#!/usr/bin/env Rscript
# Find genes around significant SNPs (refactored)
# Outputs to specified paths

library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
################
GFF = args[1]
SIGSNPS = args[2]
GFF_FEATURE = args[3]
DISTANCE = args[4] %>% as.numeric
PROMOTER_LENGTH = args[5] %>% as.numeric
ALL_SNPS = args[6]  # vcfsnp file
CPU = args[7] %>% as.numeric
OUTPUT_GENES = args[8]
OUTPUT_COLLAPSED = args[9]
################

message('INFO: Finding genes around significant SNPs')
message(paste0('INFO: Distance: ', DISTANCE))
message(paste0('INFO: GFF feature: ', GFF_FEATURE))

######################## Functions

# Extract gene ID from Parent or ID field
extract_gene_id <- function(attr) {
    id <- str_extract(attr, "(?<=Parent=)[^;]+")
    id <- ifelse(is.na(id), str_extract(attr, "(?<=ID=)[^;]+"), id)
    id <- str_remove(id, "\\.[0-9]+(_exon_[0-9]+)?$")
    return(id)
}

FUN_find_overlapped_genes <- function(GFF, gff_feature, sigSNPs, DISTANCE, CPU) {

    snps_sig <- fread(sigSNPs)
    if ('chr' %in% colnames(snps_sig)) snps_sig$chr <- as.character(snps_sig$chr)

    if (nrow(snps_sig) == 0) {
        message('WARNING: No significant SNPs to process')
        return(data.table())
    }

    traits <- snps_sig$trait %>% unique

    # Load GFF
    genes_gff <- fread(cmd = paste("grep -v '#'", GFF), header = F) %>%
        dplyr::filter(V3 == !!gff_feature) %>%
        dplyr::select(V1, V4, V5, V9) %>%
        setNames(c('chr', 'start', 'end', 'description')) %>%
        dplyr::mutate(chr = as.character(chr),
                      id = description %>% extract_gene_id)

    # Extract description fields
    fields <- genes_gff$description[1] %>%
        str_extract_all("(?<=^|;)[^=;]+(?==)", simplify = TRUE) %>%
        as.character() %>%
        unique() %>%
        .[. != ""]

    genes_gff <- genes_gff %>%
        tidyr::separate(col = "description", into = fields, sep = ";", fill = "right", extra = "drop") %>%
        dplyr::mutate(across(all_of(fields), ~ na_if(.x, "NA"))) %>%
        dplyr::mutate(start = as.integer(start), end = as.integer(end))

    genes_gr <- genes_gff %>% GRanges()

    genes_extended_gr <- GRanges(
        seqnames = seqnames(genes_gr),
        ranges = IRanges(
            start = pmax(1, start(genes_gr) - DISTANCE / 2),
            end = end(genes_gr) + DISTANCE / 2
        )
    )

    # Convert SNPs to GRanges
    snps_sig_gr <- snps_sig %>%
        dplyr::mutate(start = pos, end = pos) %>%
        dplyr::select(chr, start, end, trait, method) %>%
        dplyr::mutate(SNPID = paste0(chr, ':', start)) %>%
        GRanges()

    # Find overlaps
    genes_overlapped <- mclapply(1:length(snps_sig_gr), function(i) {
        overlaps <- findOverlaps(snps_sig_gr[i], genes_extended_gr)

        genes_gff[subjectHits(overlaps), ] %>%
            dplyr::mutate(
                SNPID = snps_sig_gr$SNPID[i],
                trait = snps_sig_gr$trait[i],
                method = snps_sig_gr$method[i]
            )
    }, mc.cores = CPU) %>%
        do.call(rbind, .)

    if (is.null(genes_overlapped) || nrow(genes_overlapped) == 0) {
        return(data.table())
    }

    genes_overlapped <- genes_overlapped %>%
        dplyr::select(trait, SNPID, method, everything())

    return(genes_overlapped)
}

# Count SNPs in exons and promoter regions
count_exon_promoter_snps <- function(genesID, snps, GFF, promoter_len, gff_feature) {

    if (length(genesID) == 0) {
        return(data.table(id = character(), promoterSNP_N = integer(),
                          promoterSNP_IDs = character(), exonSNP_N = integer(),
                          exonSNP_IDs = character()))
    }

    gff_like_dt <- fread(cmd = paste("grep -v '#'", GFF)) %>%
        dplyr::select(V1, V3, V4, V5, V9) %>%
        setNames(c('chr', 'feature', 'start', 'end', 'description')) %>%
        dplyr::filter(feature %in% c(!!gff_feature, 'exon')) %>%
        dplyr::mutate(id = description %>% extract_gene_id) %>%
        dplyr::select(-description)

    exons_gr <- gff_like_dt %>%
        dplyr::filter(id %in% !!genesID & feature == 'exon') %>%
        dplyr::select(-feature) %>%
        GRanges()

    promoters_gr <- gff_like_dt %>%
        dplyr::filter(id %in% !!genesID & feature == !!gff_feature) %>%
        dplyr::mutate(end = start, start = start - promoter_len) %>%
        dplyr::select(-feature) %>%
        GRanges()

    snps_gr <- snps %>%
        dplyr::mutate(start = pos, end = pos) %>%
        dplyr::select(chr, start, end) %>%
        GRanges()

    # Exon overlaps
    res_exon <- if (length(exons_gr) > 0 && length(snps_gr) > 0) {
        findOverlapPairs(snps_gr, exons_gr) %>%
            as.data.table %>%
            dplyr::select(first.seqnames, first.start, second.id) %>%
            setNames(c('chr', 'pos', 'id')) %>%
            dplyr::mutate(SNPID = paste0(chr, ':', pos)) %>%
            group_by(id) %>%
            dplyr::summarise(exonSNP_N = n(), exonSNP_IDs = paste(SNPID, collapse = ','))
    } else {
        data.table(id = character(), exonSNP_N = integer(), exonSNP_IDs = character())
    }

    # Promoter overlaps
    res_promoter <- if (length(promoters_gr) > 0 && length(snps_gr) > 0) {
        findOverlapPairs(snps_gr, promoters_gr) %>%
            as.data.table %>%
            dplyr::select(first.seqnames, first.start, second.id) %>%
            setNames(c('chr', 'pos', 'id')) %>%
            dplyr::mutate(SNPID = paste0(chr, ':', pos)) %>%
            group_by(id) %>%
            dplyr::summarise(promoterSNP_N = n(), promoterSNP_IDs = paste(SNPID, collapse = ','))
    } else {
        data.table(id = character(), promoterSNP_N = integer(), promoterSNP_IDs = character())
    }

    return(full_join(res_promoter, res_exon, by = 'id'))
}

#################### Main

message('INFO: Loading all SNPs')
allsnps <- fread(ALL_SNPS) %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chr = as.character(chr))
message(allsnps %>% str)

message('INFO: Finding overlapped genes')
genes_overlapped_dt <- FUN_find_overlapped_genes(
    GFF = GFF,
    gff_feature = GFF_FEATURE,
    sigSNPs = SIGSNPS,
    DISTANCE = DISTANCE,
    CPU = CPU
)

if (nrow(genes_overlapped_dt) == 0) {
    message('WARNING: No genes found around significant SNPs')

    # Create empty output files
    data.table() %>% fwrite(OUTPUT_GENES, sep = '\t')
    data.table() %>% fwrite(OUTPUT_COLLAPSED, sep = '\t')

} else {
    message(genes_overlapped_dt %>% str)

    message('INFO: Counting exon/promoter SNPs')
    exon_promoter_snp_counts <- count_exon_promoter_snps(
        genes_overlapped_dt$id,
        allsnps,
        GFF,
        PROMOTER_LENGTH,
        GFF_FEATURE
    )
    message(exon_promoter_snp_counts %>% str)

    message('INFO: Joining results')
    res_dt <- genes_overlapped_dt %>%
        left_join(exon_promoter_snp_counts, by = 'id') %>%
        dplyr::arrange(desc(exonSNP_N), desc(promoterSNP_N)) %>%
        as.data.table

    # Collapse by gene
    res_dt_collapsed <- res_dt[, lapply(.SD, function(x) {
        if (length(unique(x)) == 1) {
            as.character(x[1])
        } else {
            paste(x, collapse = ",")
        }
    }), by = id] %>%
        dplyr::mutate(
            trait_N = trait %>% str_split(',') %>% sapply(length),
            SNPID_N = SNPID %>% str_split(',') %>% sapply(length)
        ) %>%
        dplyr::arrange(chr, start)

    # Save
    message(paste0('INFO: Saving to ', OUTPUT_GENES))
    res_dt %>% fwrite(OUTPUT_GENES, col.names = TRUE, sep = '\t')

    message(paste0('INFO: Saving collapsed to ', OUTPUT_COLLAPSED))
    res_dt_collapsed %>% fwrite(OUTPUT_COLLAPSED, col.names = TRUE, sep = '\t')
}

message('INFO: Complete')
