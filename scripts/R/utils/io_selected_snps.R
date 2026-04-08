# io_selected_snps.R — read/write helpers for selected_snps tables
# Usage: source("/pipeline/scripts/R/utils/io_selected_snps.R")
# Requires: data.table loaded by the sourcing script

# Write the significant SNPs table produced by find_sig_snps.R.
write_selected_snps <- function(df, path) {
    data.table::fwrite(df, path, sep = '\t', quote = FALSE)
}

# Read a selected_snps table with correct column types.
# chr forced to character; SNPID as character.
read_selected_snps <- function(path) {
    data.table::fread(path, sep = '\t', header = TRUE,
                      colClasses = c("chr" = "character", "SNPID" = "character"))
}
