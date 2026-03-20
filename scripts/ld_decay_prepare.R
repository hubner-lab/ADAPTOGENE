#!/usr/bin/env Rscript
# Prepare sample lists and per-chromosome VCFs for PopLDdecay LD decay analysis
#
# Creates:
#   - One sample list file per group (group_name.txt, one sample per line)
#   - manifest.tsv with group names and sample counts
#   - Per-chromosome VCFs (if scope includes per_chromosome)
#   - chromosomes.txt listing all chromosome names

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
################
VCF_PATH        <- args[1]  # Filtered VCF
METADATA_PATH   <- args[2]  # metadata.tsv (site, sample, lat, lon, ...)
GROUP_BY        <- args[3]  # "site" or "cluster"
MIN_SAMPLES     <- as.integer(args[4])
SCOPE           <- args[5]  # "genome_wide", "per_chromosome", or "both"
SAMPLE_LISTS_DIR <- args[6]
CHR_VCFS_DIR    <- args[7]
CLUSTERS_PATH   <- args[8]  # clusters_K{k}.tsv or "NULL"
################

message("INFO: Preparing LD decay analysis")
message(paste0("INFO: Group by: ", GROUP_BY))
message(paste0("INFO: Min samples per group: ", MIN_SAMPLES))
message(paste0("INFO: Scope: ", SCOPE))

dir.create(SAMPLE_LISTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHR_VCFS_DIR, recursive = TRUE, showWarnings = FALSE)

#=============================================================================
# READ METADATA AND ASSIGN GROUPS
#=============================================================================
meta <- fread(METADATA_PATH, colClasses = c("site" = "character", "sample" = "character"))
# Standardize column names
colnames(meta)[1:2] <- c("site", "sample")

if (GROUP_BY == "site") {
    groups <- meta[, .(sample, group = site)]
} else if (GROUP_BY == "cluster") {
    if (CLUSTERS_PATH == "NULL") stop("Clusters file required when group_by='cluster'")
    clusters <- fread(CLUSTERS_PATH)
    # clusters_K{k}.tsv has: sample, cluster (Q-matrix assignment)
    # Find max-Q cluster per sample
    q_cols <- setdiff(colnames(clusters), "sample")
    clusters[, cluster := q_cols[max.col(.SD, ties.method = "first")], .SDcols = q_cols]
    groups <- clusters[, .(sample, group = cluster)]
}

message(paste0("INFO: Total samples: ", nrow(groups)))

#=============================================================================
# FILTER GROUPS BY MIN_SAMPLES
#=============================================================================
group_sizes <- groups[, .N, by = group]
valid_groups <- group_sizes[N >= MIN_SAMPLES, group]
skipped_groups <- group_sizes[N < MIN_SAMPLES]

if (nrow(skipped_groups) > 0) {
    for (i in seq_len(nrow(skipped_groups))) {
        message(paste0("WARNING: Skipping group '", skipped_groups$group[i],
                       "' (", skipped_groups$N[i], " samples < min_samples=", MIN_SAMPLES, ")"))
    }
}

if (length(valid_groups) == 0) {
    message("WARNING: No groups meet the min_samples threshold. Only 'All' group will be analyzed.")
}

#=============================================================================
# WRITE SAMPLE LIST FILES
#=============================================================================
# Write "All" sample list (all samples)
all_samples <- groups$sample
writeLines(all_samples, file.path(SAMPLE_LISTS_DIR, "All.txt"))
message(paste0("INFO: Written All.txt (", length(all_samples), " samples)"))

# Write per-group sample lists
for (grp in valid_groups) {
    grp_samples <- groups[group == grp, sample]
    # Clean group name for filename (replace spaces/special chars)
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", grp)
    writeLines(grp_samples, file.path(SAMPLE_LISTS_DIR, paste0(safe_name, ".txt")))
    message(paste0("INFO: Written ", safe_name, ".txt (", length(grp_samples), " samples)"))
}

# Write manifest
manifest <- data.table(
    group = c("All", valid_groups),
    n_samples = c(length(all_samples), group_sizes[group %in% valid_groups, N])
)
fwrite(manifest, file.path(SAMPLE_LISTS_DIR, "manifest.tsv"), sep = "\t")
message(paste0("INFO: Written manifest.tsv (", nrow(manifest), " groups)"))

#=============================================================================
# SPLIT VCF BY CHROMOSOME (if needed)
#=============================================================================
if (SCOPE %in% c("per_chromosome", "both")) {
    message("INFO: Splitting VCF by chromosome")

    # Get chromosome list from VCF (using grep/awk — works without index)
    chr_cmd <- paste0("grep -v '^#' ", VCF_PATH, " | awk -F'\\t' '{print $1}' | sort -u")
    chr_list <- system(chr_cmd, intern = TRUE)
    chr_list <- chr_list[chr_list != ""]

    if (length(chr_list) == 0) stop("No chromosomes found in VCF")
    message(paste0("INFO: Found ", length(chr_list), " chromosomes: ", paste(chr_list, collapse = ", ")))

    # Write chromosomes.txt
    writeLines(chr_list, file.path(CHR_VCFS_DIR, "chromosomes.txt"))

    # Extract VCF header
    header_cmd <- paste0("grep '^#' ", VCF_PATH)
    vcf_header <- system(header_cmd, intern = TRUE)

    # Split VCF per chromosome using awk (no index needed)
    for (chr in chr_list) {
        out_vcf <- file.path(CHR_VCFS_DIR, paste0(chr, ".vcf"))
        writeLines(vcf_header, out_vcf)
        # Append data lines matching this chromosome
        cmd <- paste0("awk -F'\\t' '$1 == \"", chr, "\"' ", VCF_PATH, " >> ", out_vcf)
        system(cmd)
        n_lines <- as.integer(system(paste0("grep -cv '^#' ", out_vcf), intern = TRUE))
        message(paste0("INFO: Written ", chr, ".vcf (", n_lines, " SNPs)"))
    }
} else {
    # Write empty chromosomes.txt so downstream doesn't error
    writeLines(character(0), file.path(CHR_VCFS_DIR, "chromosomes.txt"))
}

message("INFO: LD decay preparation complete")
