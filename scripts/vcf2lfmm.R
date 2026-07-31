library(LEA)

args = commandArgs(trailingOnly=TRUE)
VCF=args[1] # VCF filename
FORCE = ifelse(length(args) == 2, args[2], F) #  T or TRUE

vcf2lfmm(VCF, force = as.logical(FORCE))

# Record whether this LFMM has any missing (9-coded) genotypes at all -- computed
# once here, reused by impute.R downstream instead of re-scanning the matrix.
# Simulated datasets (e.g. SLiM output) can be genuinely 100% complete, a case
# LEA::impute() does not support (it errors when there is nothing to impute).
LFMM_PATH <- sub("\\.vcf$", ".lfmm", VCF)
NMISSING_PATH <- sub("\\.vcf$", ".lfmm_nmissing", VCF)
has_missing <- FALSE
con <- file(LFMM_PATH, open = "rt")
repeat {
    lines <- readLines(con, n = 1000)
    if (length(lines) == 0) break
    if (any(grepl("(^|\\s)9(\\s|$)", lines))) { has_missing <- TRUE; break }
}
close(con)
writeLines(as.character(as.integer(has_missing)), NMISSING_PATH)
