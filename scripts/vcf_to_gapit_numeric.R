#!/usr/bin/env Rscript
# Convert imputed VCF to GAPIT numeric format (GD + GM files)
# GD: Taxa + numeric genotype columns (0/1/2)
# GM: Name, Chromosome, Position

library(vcfR)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
########################
VCF_FILE   <- args[1]  # Imputed full VCF
OUTPUT_GD  <- args[2]  # Output GD file
OUTPUT_GM  <- args[3]  # Output GM file
########################

message("INFO: Reading VCF: ", VCF_FILE)
vcf <- read.vcfR(VCF_FILE, verbose = FALSE)

# Extract sample names
samples <- colnames(vcf@gt)[-1]  # Remove FORMAT column
message("INFO: ", length(samples), " samples")

# Extract SNP info for GM
chrom <- vcf@fix[, "CHROM"]
pos <- as.integer(vcf@fix[, "POS"])
snp_names <- paste0(chrom, ":", pos)
message("INFO: ", length(snp_names), " SNPs")

# Build GM: Name, Chromosome, Position
gm <- data.table(Name = snp_names, Chromosome = chrom, Position = pos)

# Extract GT matrix and convert to numeric 0/1/2
gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# Convert GT strings to numeric dosage
gt_to_numeric <- function(gt) {
    # Handle various GT formats: 0/0, 0/1, 1/1, 0|0, 0|1, 1|1, ./.
    ifelse(is.na(gt) | gt == "./." | gt == ".|.", NA_integer_,
    ifelse(gt == "0/0" | gt == "0|0", 0L,
    ifelse(gt == "0/1" | gt == "1/0" | gt == "0|1" | gt == "1|0", 1L,
    ifelse(gt == "1/1" | gt == "1|1", 2L, NA_integer_))))
}

message("INFO: Converting genotypes to numeric")
gd_mat <- apply(gt_matrix, 2, gt_to_numeric)

# Build GD: Taxa + SNP columns (transposed: samples as rows, SNPs as columns)
gd <- data.table(Taxa = samples)
gd <- cbind(gd, as.data.table(t(gd_mat)))
setnames(gd, c("Taxa", snp_names))

message("INFO: GD dimensions: ", nrow(gd), " samples x ", ncol(gd) - 1, " SNPs")
message("INFO: GM dimensions: ", nrow(gm), " SNPs")

# Write outputs
fwrite(gd, OUTPUT_GD, sep = "\t")
fwrite(gm, OUTPUT_GM, sep = "\t")

message("INFO: Wrote GD to ", OUTPUT_GD)
message("INFO: Wrote GM to ", OUTPUT_GM)
message("INFO: VCF to GAPIT numeric conversion complete")
