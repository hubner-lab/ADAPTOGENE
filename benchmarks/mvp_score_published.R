#!/usr/bin/env Rscript
# mvp_score_published.R -- score the MVP authors' OWN per-SNP p-values through OUR harness.
#
# The deposit's muts_full table ships the p-values Lotterhos 2023 actually analysed
# (`LEA3.2_lfmm2_mlog10P_tempenv` / `_salenv` for LFMM, `RDA_mlog10P` for RDA). Re-scoring
# them with benchmarks/eval_detection.R removes every scoring-convention difference, so any
# remaining gap against our pipeline's AUC-PR is a difference in the p-values themselves --
# i.e. in how the association was run -- not in how detection was counted.
#
# Emits wide p-value tables in the harness's format (SNPID/chr/pos/<traits>), using the SAME
# coordinate derivation and duplicate collapse as convert_mvp.R so the keys join to the same
# truth tables.
#
# Usage:
#   Rscript benchmarks/mvp_score_published.R --seed=1232548 [--raw-dir=...] [--out-dir=...]

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || is.na(x) || identical(x, "")) y else x
A <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
  a <- sub("^--", "", a)
  A[[gsub("-", "_", sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
}
SEED    <- A$seed
RAW_DIR <- A$raw_dir %||% "/pipeline/data/mvp/raw"
OUT_DIR <- A$out_dir %||% file.path("/pipeline/benchmarks/mvp_eval/published", paste0("MVP", SEED))
if (is.null(SEED)) stop("Usage: mvp_score_published.R --seed=N")

LG_LEN <- 50000L
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

muts_path <- file.path(RAW_DIR, "mutations", paste0(SEED, "_Rout_muts_full.txt.gz"))
con <- gzfile(muts_path, open = "rt"); h2 <- readLines(con, n = 2L); close(con)
hdr <- scan(text = h2[1], what = "", quiet = TRUE)
r1  <- scan(text = h2[2], what = "", quiet = TRUE)
cn  <- if (length(r1) == length(hdr) + 1L) c(".rowname", hdr) else hdr
muts <- fread(cmd = paste("zcat", shQuote(muts_path)), header = FALSE, skip = 1L,
              sep = " ", col.names = cn, showProgress = FALSE)

# Same coordinate derivation as convert_mvp.R: pos_pyslim is 1-based genome-global, the
# boundary locus belongs to the LOWER linkage group.
muts[, LG := as.integer(LG)]
muts[, chr := LG]
muts[, pos := as.integer(pos_pyslim - (LG - 1L) * LG_LEN)]
stopifnot(all(muts$pos >= 1L & muts$pos <= LG_LEN))

# Same duplicate collapse, keeping the causal row.
muts[, .row := .I][, .is_causal := as.logical(causal)]
setorder(muts, chr, pos, -.is_causal)
keep <- sort(muts[, .SD[1L], by = .(chr, pos)]$.row)
muts <- muts[order(.row)][.row %in% keep]
setorder(muts, chr, pos)
message("INFO: ", nrow(muts), " loci after collapse")

mlog10_to_p <- function(x) {
  p <- 10^(-as.numeric(x))
  p[!is.finite(p)] <- NA_real_
  pmin(pmax(p, .Machine$double.xmin), 1)
}

base <- muts[, .(SNPID = paste0(chr, ":", pos), chr, pos)]

lfmm <- copy(base)
lfmm[, bio_1 := mlog10_to_p(muts$`LEA3.2_lfmm2_mlog10P_tempenv`)]
lfmm[, bio_2 := mlog10_to_p(muts$`LEA3.2_lfmm2_mlog10P_salenv`)]
fwrite(lfmm, file.path(OUT_DIR, "published_LFMM_pvalues.tsv"), sep = "\t")

rda <- copy(base)
rda[, climate_multivariate := mlog10_to_p(muts$RDA_mlog10P)]
fwrite(rda, file.path(OUT_DIR, "published_RDA_pvalues.tsv"), sep = "\t")

rda_c <- copy(base)
rda_c[, climate_multivariate := mlog10_to_p(muts$RDA_mlog10P_corr)]
fwrite(rda_c, file.path(OUT_DIR, "published_RDA_structcorr_pvalues.tsv"), sep = "\t")

message("INFO: wrote published p-value tables to ", OUT_DIR)
for (f in c("published_LFMM_pvalues.tsv", "published_RDA_pvalues.tsv",
            "published_RDA_structcorr_pvalues.tsv")) {
  d <- fread(file.path(OUT_DIR, f))
  message("  ", f, ": ", nrow(d), " rows, non-NA per trait: ",
          paste(sapply(setdiff(names(d), c("SNPID", "chr", "pos")),
                       function(x) paste0(x, "=", sum(!is.na(d[[x]])))), collapse = ", "))
}
