#!/usr/bin/env Rscript
# Assemble per-trait p-value files into a single wide table.
#
# Each per-trait file has columns: SNPID, chr, pos, <trait_name>
# The assembled wide table has: SNPID, chr, pos, trait1, trait2, ...
#
# Usage:
#   Rscript assemble_pvalues.R file1.tsv file2.tsv ... output.tsv
#
# The last argument is the output file; all preceding arguments are inputs.

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
n <- length(args)
if (n < 2) stop("Usage: assemble_pvalues.R <input1> [input2 ...] <output>")

trait_files <- args[seq_len(n - 1L)]
out_file    <- args[n]

message("INFO: Assembling ", length(trait_files), " per-trait file(s) → ", out_file)

tables <- lapply(trait_files, function(f) {
    if (!file.exists(f)) stop("Input file not found: ", f)
    dt <- data.table::fread(f, sep = "\t", header = TRUE,
                            colClasses = c(chr = "character", SNPID = "character"))
    fixed <- c("SNPID", "chr", "pos")
    if (!all(fixed %in% names(dt)))
        stop("Missing required columns (SNPID/chr/pos) in: ", f)
    dt
})

# Merge all on SNPID/chr/pos (full outer join preserves all SNPs even if any trait
# had a minor filtering difference)
result <- Reduce(function(a, b) {
    merge(a, b, by = c("SNPID", "chr", "pos"), all = TRUE)
}, tables)

result <- result[order(chr, pos)]

message("INFO: Assembled ", ncol(result) - 3L, " trait(s) × ",
        nrow(result), " SNPs")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(result, out_file, sep = "\t")
message("INFO: Wrote ", out_file)
