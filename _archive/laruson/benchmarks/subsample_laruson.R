#!/usr/bin/env Rscript
# =============================================================================
# subsample_laruson.R — build a balanced N-per-site subsample of the Láruson
# converted dataset as a NEW project, leaving data/laruson_converted/ untouched.
#
# Why: the full replicate is 10,000 samples (100 sites x 100). RDA's anova.cca
# was OOM-killed at 20 g and 150 g on that cohort and had to drop permutations
# 999 -> 99. At 1,000 samples the genotype matrix is ~10x smaller, the full
# hyperparameter ladder is affordable, and all 100 sites / both environmental
# axes are preserved — only within-deme replication shrinks.
#
# THE SEED IS NOT ARBITRARY. Re-applying Filter.maf to a 10x smaller cohort can
# drop causal loci that sat near the threshold. This script therefore measures,
# for each candidate seed, how many of the causal loci that pass MAF at the full
# cohort still pass on the subsample, and selects the first seed that retains
# ALL of them. That keeps the 1k and 10k truth sets identical, which is what
# makes hyperparameter rungs comparable across the two stages. If no candidate
# seed qualifies, the script stops rather than silently shipping a smaller
# truth set.
#
# Usage:
#   Rscript subsample_laruson.R [--n-per-site=10] [--maf=0.05] \
#       [--seeds=1,2,3,4,5,6,7,8,9,10] [--in-dir=...] [--out-dir=...] \
#       [--truth=...] [--vcf-name=laruson.vcf]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")

parse_args <- function(argv) {
    out <- list()
    for (a in argv) {
        if (!grepl("^--[A-Za-z0-9._-]+=", a)) stop("Malformed argument: ", a)
        out[[gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", a))]] <- sub("^--[^=]+=", "", a)
    }
    out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
opt  <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]

IN_DIR      <- opt("in_dir",  file.path(PIPELINE_ROOT, "data/laruson_converted"))
OUT_DIR     <- opt("out_dir", file.path(PIPELINE_ROOT, "data/laruson_1k"))
TRUTH       <- opt("truth",   file.path(IN_DIR, "causal_loci.tsv"))
VCF_NAME    <- opt("vcf_name", "laruson.vcf")
N_PER_SITE  <- as.integer(opt("n_per_site", "10"))
MAF_CUTOFF  <- as.numeric(opt("maf", "0.05"))
SEEDS       <- as.integer(strsplit(opt("seeds", paste(1:10, collapse = ",")), ",")[[1]])

VCF_IN <- file.path(IN_DIR, VCF_NAME)
for (f in c(VCF_IN, TRUTH, file.path(IN_DIR, "metadata.tsv"))) {
    if (!file.exists(f)) stop("Required input missing: ", f)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------- metadata
# Sample and site IDs are numeric-looking but must stay character or the joins
# against the VCF header break (CLAUDE.md R Script Conventions).
meta <- fread(file.path(IN_DIR, "metadata.tsv"),
              colClasses = c(site = "character", sample = "character"))
message("INFO: full cohort ", nrow(meta), " samples across ", uniqueN(meta$site), " sites")

site_sizes <- meta[, .N, by = site]
if (any(site_sizes$N < N_PER_SITE)) {
    stop("Some sites have fewer than ", N_PER_SITE, " samples: ",
         paste(site_sizes[N < N_PER_SITE, site], collapse = ", "))
}

# ------------------------------------------- causal genotypes from the raw VCF
truth  <- fread(TRUTH, colClasses = c(chr = "character"))
causal <- truth[category == "causal"]
message("INFO: ", nrow(causal), " causal loci in truth table")

targets <- tempfile(fileext = ".txt")
fwrite(causal[, .(chr, pos)], targets, sep = "\t", col.names = FALSE)

message("INFO: extracting causal-locus genotypes via bcftools query (single pass)")
gt_raw <- system2("bcftools",
                  c("query", "-T", shQuote(targets),
                    "-f", shQuote("%CHROM\\t%POS[\\t%GT]\\n"), shQuote(VCF_IN)),
                  stdout = TRUE)
if (!length(gt_raw)) stop("bcftools query returned no rows — check contig naming in ", VCF_IN)

gt <- fread(text = paste(gt_raw, collapse = "\n"), header = FALSE,
            colClasses = list(character = 1))
setnames(gt, c("chr", "pos", meta$sample))
message("INFO: recovered ", nrow(gt), " / ", nrow(causal), " causal loci from the VCF")
if (nrow(gt) != nrow(causal)) {
    stop("Causal-locus recovery incomplete — refusing to choose a seed on partial data")
}

# GT strings are unphased diploid ("0/1"); dosage = count of ALT alleles.
dose <- as.matrix(gt[, meta$sample, with = FALSE])
dose <- matrix(lengths(regmatches(dose, gregexpr("1", dose, fixed = TRUE))),
               nrow = nrow(gt), dimnames = list(paste(gt$chr, gt$pos, sep = ":"), meta$sample))

maf_of <- function(cols) {
    f <- rowSums(dose[, cols, drop = FALSE]) / (2 * length(cols))
    pmin(f, 1 - f)
}

maf_full  <- maf_of(meta$sample)
keep_full <- names(maf_full)[maf_full >= MAF_CUTOFF]
message("INFO: causal loci passing MAF >= ", MAF_CUTOFF, " at full cohort: ",
        length(keep_full), " / ", nrow(causal))

# ------------------------------------------------------------ seed selection
by_site <- split(meta$sample, meta$site)

draw <- function(seed) {
    set.seed(seed)
    unlist(lapply(by_site, function(s) sample(s, N_PER_SITE)), use.names = FALSE)
}

chosen <- NA_integer_
report <- data.table(seed = integer(), retained = integer(), passing = integer())
for (s in SEEDS) {
    sel <- draw(s)
    m   <- maf_of(sel)
    retained <- sum(m[keep_full] >= MAF_CUTOFF)
    passing  <- sum(m >= MAF_CUTOFF)
    report <- rbind(report, data.table(seed = s, retained = retained, passing = passing))
    message(sprintf("  seed %-3d retains %d/%d of the full-cohort causal set (%d causal pass overall)",
                    s, retained, length(keep_full), passing))
    if (is.na(chosen) && retained == length(keep_full)) chosen <- s
}
if (is.na(chosen)) {
    stop("No candidate seed retained all ", length(keep_full), " full-cohort causal loci. ",
         "Widen --seeds, or accept a smaller truth set explicitly (not done silently).")
}
message("INFO: selected seed ", chosen, " — retains all ", length(keep_full),
        " full-cohort causal loci")

keep_samples <- draw(chosen)
stopifnot(length(keep_samples) == N_PER_SITE * uniqueN(meta$site))

# Preserve original VCF column order so the subset VCF and metadata stay aligned.
keep_ordered <- meta$sample[meta$sample %in% keep_samples]

# --------------------------------------------------------------- write outputs
writeLines(keep_ordered, file.path(OUT_DIR, "keep_samples.txt"))
fwrite(meta[sample %in% keep_samples][match(keep_ordered, sample)],
       file.path(OUT_DIR, "metadata.tsv"), sep = "\t")

for (f in c("environments_present.tsv", "environments_future.tsv",
            "laruson_minimal.gff3", "causal_loci.tsv")) {
    src <- file.path(IN_DIR, f)
    if (file.exists(src)) {
        file.copy(src, file.path(OUT_DIR, f), overwrite = TRUE)
        message("INFO: copied ", f, " verbatim (site-level or position-level — unchanged by subsampling)")
    }
}

fwrite(report, file.path(OUT_DIR, "seed_selection.tsv"), sep = "\t")
writeLines(c(
    paste0("seed\t", chosen),
    paste0("n_per_site\t", N_PER_SITE),
    paste0("n_samples\t", length(keep_ordered)),
    paste0("n_sites\t", uniqueN(meta$site)),
    paste0("maf_cutoff\t", MAF_CUTOFF),
    paste0("causal_pass_full_cohort\t", length(keep_full)),
    paste0("causal_retained_subsample\t", report[seed == chosen, retained])
), file.path(OUT_DIR, "subsample_provenance.tsv"))

message("INFO: subsetting VCF with bcftools view (", length(keep_ordered), " samples)")
vcf_out <- file.path(OUT_DIR, VCF_NAME)
# -I/--no-update is REQUIRED, not an optimisation. Three rows in this replicate
# (1:84913, 1:109742, 1:167837 — none causal, none surviving the MAF filter)
# carry a GT allele index of 2 while declaring a single ALT, so bcftools'
# AC/AN recalculation aborts with "Incorrect allele (2)". This is a third
# pyslim VCF-writer quirk alongside the two in docs/laruson_dataset_notes.md
# section 3 (blank REF on causal rows, duplicate POS). Skipping the INFO update
# passes every record through verbatim, which is also what keeps this subset
# directly comparable to the full-cohort run that plink already consumed.
rc <- system2("bcftools",
              c("view", "-I", "-S", shQuote(file.path(OUT_DIR, "keep_samples.txt")),
                "-o", shQuote(vcf_out), "-O", "v", shQuote(VCF_IN)))
if (rc != 0) stop("bcftools view failed with exit code ", rc)

message("INFO: wrote ", vcf_out)
message("INFO: done — subsample project inputs are in ", OUT_DIR)
