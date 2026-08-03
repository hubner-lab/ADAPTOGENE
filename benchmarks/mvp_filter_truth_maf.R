#!/usr/bin/env Rscript
# =============================================================================
# mvp_filter_truth_maf.R -- build MAF-matched truth tables for a parallel MAF arm.
#
# WHY THIS EXISTS, AND WHY IT IS NOT OPTIONAL
# -------------------------------------------
# A MAF arm changes which SNPs exist in the scored table. If the truth table still lists loci
# the VCF no longer contains, every recall number is divided by a denominator that includes
# loci nothing could possibly have found. That is precisely the defect that disqualified the
# Laruson benchmark (_archive/laruson/README.md): truth defined at one MAF, VCF at another, so
# recall was capped by construction rather than by method quality. Re-running at Filter.maf
# 0.05 without re-deriving the truth would rebuild that defect on purpose.
#
# So the rule is: filter the truth table by the SAME MAF rule the pipeline applies, dropping
# every row -- causal, linked_neutral and background_neutral alike -- whose SNP will not
# survive. Dropping only the causal rows would leave the false-positive denominator counting
# background SNPs that are not in the VCF either.
#
# MAF SOURCE, in order of authority:
#   --vcf   the ACTUAL filtered VCF the pipeline produced. Use this whenever it exists; it
#           settles every boundary question by construction.
#   loci_freq.tsv  the converter's own per-locus MAF over the same 1000 individuals. Correct
#           in the interior, and the only option before mode=processing has run.
#
# The boundary matters and was measured, not assumed: plink's `--maf 0.05` EXCLUDES variants
# *below* the threshold, i.e. it keeps MAF >= 0.05. Filtering the truth with a strict `>`
# kept 7 823 loci where plink kept 7 841 -- 18 SNPs sit at exactly 0.05. So the comparison
# here is `>=`, matching plink, and --vcf is preferred when available.
#
# The truth-join gate in mvp_score_sweep.sh remains the backstop either way: it checks the
# filtered counts against what actually turned out testable, and a disagreement there is a
# real finding, not noise to suppress.
#
# Usage:
#   Rscript mvp_filter_truth_maf.R --seed=1232548 --maf=0.05 --out-suffix=M05
#       [--vcf=/path/to/filtered.vcf] [--data-dir=/pipeline/data/mvp]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args <- parse_kv_args(commandArgs(trailingOnly = TRUE))
req  <- function(k) {
    if (is.null(args[[k]]) || !nzchar(args[[k]])) stop("Missing required --", gsub("_", "-", k))
    args[[k]]
}
opt  <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]

SEED     <- req("seed")
MAF      <- as.numeric(req("maf"))
OUT_SFX  <- req("out_suffix")
DATA_DIR <- opt("data_dir", file.path(PIPELINE_ROOT, "data/mvp"))

SRC <- file.path(DATA_DIR, paste0("MVP", SEED))
DST <- file.path(DATA_DIR, paste0("MVP", SEED, OUT_SFX))
dir.create(DST, recursive = TRUE, showWarnings = FALSE)

freq <- fread(file.path(SRC, "loci_freq.tsv"), colClasses = c(chr = "character"))
stopifnot(all(c("chr", "pos", "maf") %in% names(freq)))
freq[, key := paste(chr, pos, sep = ":")]

# `>=`, not `>`: plink --maf excludes variants BELOW the threshold. See the header note.
keep_freq <- freq[maf >= MAF, key]
VCF <- opt("vcf", "")
if (nzchar(VCF)) {
    if (!file.exists(VCF)) stop("--vcf given but not found: ", VCF)
    v <- fread(cmd = paste("grep -v '^#'", shQuote(VCF), "| cut -f1,2"),
               header = FALSE, colClasses = c(V1 = "character"))
    keep <- paste(v$V1, v$V2, sep = ":")
    message(sprintf("INFO: SNP set taken from the filtered VCF: %d SNPs", length(keep)))
    d <- length(setdiff(keep_freq, keep)) + length(setdiff(keep, keep_freq))
    if (d) message(sprintf(
        "INFO: loci_freq MAF >= %.3f would have kept %d (%d SNPs differ) -- the VCF wins",
        MAF, length(keep_freq), d))
} else {
    keep <- keep_freq
    message(sprintf("INFO: MAF >= %.3f keeps %d / %d SNPs (from loci_freq; pass --vcf to be exact)",
                    MAF, length(keep), nrow(freq)))
}

rows <- list()
for (f in c("truth_temp.tsv", "truth_sal.tsv", "truth_any.tsv")) {
    src <- file.path(SRC, f)
    if (!file.exists(src)) { message("SKIP (absent): ", src); next }
    tr <- fread(src, colClasses = c(chr = "character"))
    tr[, key := paste(chr, pos, sep = ":")]
    before <- tr[, .N, by = category]
    out <- tr[key %in% keep][, key := NULL]
    after <- out[, .N, by = category]
    fwrite(out, file.path(DST, f), sep = "\t")
    m <- merge(before, after, by = "category", all = TRUE, suffixes = c("_before", "_after"))
    m[is.na(N_after), N_after := 0L]
    m[, `:=`(seed = SEED, table = f, maf = MAF)]
    rows[[length(rows) + 1L]] <- m
    message(sprintf("INFO: %-16s %d -> %d rows", f, nrow(tr), nrow(out)))
}

# loci_freq is copied through, MAF-filtered, so every downstream script that keys on it
# (mvp_wza_gate.R, mvp_profile_detections.R) sees the same SNP universe as the truth tables
# rather than silently reaching back into the unfiltered arm.
fwrite(freq[key %in% keep][, key := NULL], file.path(DST, "loci_freq.tsv"), sep = "\t")

prov <- rbindlist(rows)
setcolorder(prov, c("seed", "table", "maf", "category", "N_before", "N_after"))
fwrite(prov, file.path(DST, "truth_filter_provenance.tsv"), sep = "\t")
message("\n=== truth rows kept, by category ===")
print(prov)
message("\nINFO: wrote MAF-matched truth to ", DST)
