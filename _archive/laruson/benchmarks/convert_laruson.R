#!/usr/bin/env Rscript
# convert_laruson.R — Convert one Láruson et al. 2022 SLiM simulation replicate
# (Dryad DOI 10.5061/dryad.x95x69pkk) into ADAPTOGENE's standard input set for
# mode=processing onward, using Climate.source: custom.
#
# Column/format facts used below are all confirmed empirically against the real
# archive — see docs/laruson_dataset_notes.md (Gate G1). Do not re-derive these
# from docs/laruson_automation_roadmap.md's Phase 4 spec, which used placeholder
# guesses written before the archive was available.
#
# Usage:
#   Rscript convert_laruson.R RAW_DIR OUT_DIR CASE_N SEED TP_WINDOW_KB
#
# Args:
#   1  RAW_DIR       Directory holding the extracted archive files (read-only),
#                    e.g. data/laruson/. Expects, for the given Case/seed:
#                      Case{N}_{SEED}.vcf.gz
#                      Case{N}_{SEED}_ind.txt
#                      Case{N}_{SEED}_causal_mutations_pos_filtered.txt
#                      Case{N}_{SEED}_ML_WF_CG_sum_Gen.txt
#   2  OUT_DIR       Output directory (created if absent), e.g. data/laruson_converted/
#   3  CASE_N        Case number, 1-4 (see docs/laruson_dataset_notes.md §1)
#   4  SEED          Replicate seed, e.g. 2889863491989 (see seeds_source_R.txt)
#   5  TP_WINDOW_KB  Half-window (kb) for linked_neutral classification around each
#                    causal locus (placeholder default 5 — see dataset notes §4)
#
# Outputs (all written to OUT_DIR; RAW_DIR is never modified):
#   laruson.vcf                 decompressed VCF, unmodified otherwise (already
#                                chr-normalized: single contig, name "1")
#   metadata.tsv                site, sample, latitude, longitude, bio_1, bio_2
#   environments_present.tsv    site, bio_1, bio_2
#   environments_future.tsv     site, bio_1, bio_2 (present + fixed synthetic shift —
#                                placeholder for maladaptation-machinery smoke test,
#                                NOT a real common-garden target; see dataset notes §4)
#   causal_loci.tsv              chr, pos, category (causal | linked_neutral | background_neutral)
#   fitness_timeseries.tsv       raw pass-through of _ML_WF_CG_sum_Gen.txt (population-average
#                                generation time series — NOT a per-site fitness table; kept for
#                                a future eval harness, not consumed by this adapter pass)
#   laruson_minimal.gff3         one synthetic "gene" record per causal locus (+-500bp)
#
# Sample-to-deme join: ind.txt row i (0-based, after header) == VCF sample column i
# (0-based, after the 9 fixed VCF columns). This is an assumption validated only by
# row-count equality (both files list the same simulation's individuals in the same
# PySlim run) -- no independent source-code cross-check was available in the archive.
# Hard-stops below if the counts ever disagree.

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: convert_laruson.R RAW_DIR OUT_DIR CASE_N SEED TP_WINDOW_KB")
}

RAW_DIR      <- args[1]
OUT_DIR      <- args[2]
CASE_N       <- args[3]
SEED         <- args[4]
TP_WINDOW_KB <- as.numeric(args[5])
TP_WINDOW_BP <- TP_WINDOW_KB * 1000

if (!dir.exists(RAW_DIR)) stop("ERROR: RAW_DIR not found: ", RAW_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PREFIX <- file.path(RAW_DIR, paste0("Case", CASE_N, "_", SEED))
vcf_gz_path   <- paste0(PREFIX, ".vcf.gz")
ind_path      <- paste0(PREFIX, "_ind.txt")
causal_path   <- paste0(PREFIX, "_causal_mutations_pos_filtered.txt")
cg_sum_path   <- paste0(PREFIX, "_ML_WF_CG_sum_Gen.txt")

for (f in c(vcf_gz_path, ind_path, causal_path, cg_sum_path)) {
    if (!file.exists(f)) stop("ERROR: Expected input file not found: ", f)
}

message("INFO: Converting Case", CASE_N, " seed ", SEED, " from ", RAW_DIR, " -> ", OUT_DIR)

#=============================================================================
# 1. VCF -- decompress, repair a pyslim writer quirk, drop duplicate-POS rows.
#    Already single-contig ("1"), no "chr" prefix, integer POS -- confirmed
#    against this archive, no further coordinate normalization needed.
#=============================================================================
vcf_out_path <- file.path(OUT_DIR, "laruson.vcf")
message("INFO: Decompressing + normalizing VCF -> ", vcf_out_path)

# Two known pyslim/tskit VCF-writer quirks in this archive, both confirmed empirically:
#
# (a) Causal-mutation rows: REF written blank, ALT holds a raw SLiM mutation ID
#     (e.g. "\t\t675377\t" or "\t\t674470,674654\t" when >=2 mutations stack at one
#     site) instead of a normal biallelic "0"/"1" pair. Confirmed on 102/102 causal
#     loci in the Case1 pilot replicate, 0/33541 neutral rows -- "blank REF" is a
#     100% reliable causal-row signal in this dataset (no other field is ever
#     legitimately blank). GT calls themselves are untouched/valid (0|0/0|1/1|0/1|1).
#
# (b) Duplicate POS: ~3% of positions carry >1 independently-simulated mutation at
#     the exact same site (e.g. one causal + one neutral, confirmed via distinct GT
#     patterns at pos 71772 in the Case1 replicate) -- not a formatting artifact, a
#     real recurrent-mutation coalescent outcome. This pipeline is biallelic-only
#     (see CLAUDE.md), so duplicate-POS rows must collapse to one row per position.
#     Naively keeping "whichever row comes first" risks silently discarding the
#     causal row's genotype signal in favor of a coincidentally co-located neutral
#     one -- so quirk (a)'s signal doubles as the tie-breaker for quirk (b): always
#     keep the causal (blank-REF) row when a duplicate group contains one.
MALFORMED_REF_PATTERN <- "^([^\t]*\t[^\t]*\t[^\t]*\t)\t[0-9]+(,[0-9]+)*\t"

# Decompress once to a plain-text temp file (native R gzfile() connection -- no
# external zcat/bash dependency, and no R.utils dependency for the fread calls below,
# which only support direct .gz reads via that optional package).
tmp_raw_vcf <- tempfile(tmpdir = OUT_DIR, fileext = ".raw.vcf")
in_con  <- gzfile(vcf_gz_path, open = "rt")
out_con <- file(tmp_raw_vcf, open = "wt")
while (length(chunk <- readLines(in_con, n = 100000)) > 0) writeLines(chunk, out_con)
close(in_con)
close(out_con)

# Pass 1 (lightweight): fread only POS + REF (columns 2 and 4) -- data.table skips
# parsing/storing the other 10,005 genotype columns entirely, far faster than a
# manual per-line regex scan of the full (very wide) lines.
pos_ref <- fread(tmp_raw_vcf, skip = "#CHROM", header = TRUE, select = c(2, 4),
                  col.names = c("pos", "ref"), colClasses = list(character = 4))
n_total_rows <- nrow(pos_ref)
line_dt <- data.table(idx = seq_len(n_total_rows), pos = pos_ref$pos,
                       malformed = is.na(pos_ref$ref) | pos_ref$ref == "")
setorder(line_dt, pos, -malformed)  # within a POS group, malformed (causal) row sorts first
keep_dt <- line_dt[, .SD[1L], by = pos]
keep_idx <- sort(keep_dt$idx)
n_dropped <- n_total_rows - length(keep_idx)
if (n_dropped > 0) {
    message("INFO: Dropping ", n_dropped, " duplicate-POS rows (", n_total_rows,
            " -> ", length(keep_idx), " variants) -- keeping the causal/blank-REF row ",
            "wherever a duplicate group contains one, see comment above")
}
keep_idx_set <- keep_idx  # sorted, unique -- safe for match() below

# Pass 2: re-read the decompressed temp file, apply the REF/ALT repair only to lines
# being kept, write only kept data lines (plus all header/meta lines) to the final VCF.
n_fixed <- 0L
row_counter <- 0L
in_con  <- file(tmp_raw_vcf, open = "rt")
out_con <- file(vcf_out_path, open = "wt")
while (length(chunk <- readLines(in_con, n = 100000)) > 0) {
    is_data_line <- !startsWith(chunk, "#")
    n_data_this_chunk <- sum(is_data_line)
    if (n_data_this_chunk > 0) {
        data_idx <- row_counter + seq_len(n_data_this_chunk)
        is_kept <- data_idx %in% keep_idx_set
        keep_mask <- rep(TRUE, length(chunk))
        keep_mask[is_data_line] <- is_kept
        chunk <- chunk[keep_mask]

        # Recompute is_data_line/malformed on the now-filtered chunk for the repair step
        is_data_line2 <- !startsWith(chunk, "#")
        is_malformed <- is_data_line2 & grepl("\t\t", chunk, fixed = TRUE)
        if (any(is_malformed)) {
            n_fixed <- n_fixed + sum(is_malformed)
            chunk[is_malformed] <- sub(MALFORMED_REF_PATTERN, "\\1<<FIXREF>>\t", chunk[is_malformed])
            chunk[is_malformed] <- sub("<<FIXREF>>\t", "0\t1\t", chunk[is_malformed], fixed = TRUE)
        }
        row_counter <- row_counter + n_data_this_chunk
    }
    writeLines(chunk, out_con)
}
close(in_con)
close(out_con)
unlink(tmp_raw_vcf)
if (!file.exists(vcf_out_path) || file.info(vcf_out_path)$size == 0) {
    stop("ERROR: VCF decompression failed for ", vcf_gz_path)
}
if (n_fixed > 0) {
    message("INFO: Repaired ", n_fixed, " malformed REF/ALT rows (blank REF -> \"0\", ",
            "SLiM mutation-ID ALT -> \"1\") -- see comment above")
}

# Re-scan the written file for any remaining blank field or duplicate-POS data lines --
# catches an incomplete repair/dedup and any other malformation this converter doesn't
# yet know about.
written_check <- fread(vcf_out_path, skip = "#CHROM", header = TRUE, select = c(2, 4),
                        col.names = c("pos", "ref"), colClasses = list(character = 4))
still_malformed <- sum(is.na(written_check$ref) | written_check$ref == "")
still_dup <- sum(duplicated(written_check$pos))
if (still_malformed > 0 || still_dup > 0) {
    stop("ERROR: ", still_malformed, " blank-field and ", still_dup, " duplicate-POS data ",
         "lines remain after repair/dedup -- inspect ", vcf_out_path, " manually before ",
         "trusting downstream output.")
}

# Read header line (the "#CHROM..." row) to get sample names, and the contig ID.
vcf_lines_meta <- readLines(vcf_out_path, n = 5000)
header_line <- vcf_lines_meta[startsWith(vcf_lines_meta, "#CHROM")]
if (length(header_line) == 0) stop("ERROR: No #CHROM header line found in ", vcf_out_path)
vcf_colnames <- strsplit(header_line, "\t")[[1]]
vcf_colnames[1] <- sub("^#", "", vcf_colnames[1])
vcf_samples <- vcf_colnames[-(1:9)]
n_vcf_samples <- length(vcf_samples)
message("INFO: VCF has ", n_vcf_samples, " samples")

contig_line <- vcf_lines_meta[grepl("^##contig=", vcf_lines_meta)][1]
CHR <- if (!is.na(contig_line)) sub(".*ID=([^,>]+).*", "\\1", contig_line) else "1"
message("INFO: Contig ID from VCF header: ", CHR)

#=============================================================================
# 2. ind.txt -- deme assignment per sample, via row-order join (see header note).
#=============================================================================
ind <- fread(ind_path, header = TRUE)
message("INFO: ind.txt has ", nrow(ind), " rows, columns: ", paste(names(ind), collapse = ", "))

if (nrow(ind) != n_vcf_samples) {
    stop("ERROR: ind.txt row count (", nrow(ind), ") != VCF sample count (", n_vcf_samples,
         ") -- the row-order join assumption (see docs/laruson_dataset_notes.md) does not hold ",
         "for this replicate. Stopping rather than silently misassigning samples to demes.")
}

ind[, sample := vcf_samples]
ind[, site := as.character(subpop)]

# Per-site environment: unique (site, opt0, opt1) -- one row per deme, confirmed
# constant within a deme in the source data. Fail loudly if that's ever violated.
site_env <- unique(ind[, .(site, opt0, opt1)])
if (anyDuplicated(site_env$site)) {
    stop("ERROR: opt0/opt1 not constant within at least one deme -- cannot derive a single ",
         "per-site environment table. Inspect ", ind_path, " manually.")
}
setnames(site_env, c("opt0", "opt1"), c("bio_1", "bio_2"))
setorder(site_env, site)
n_sites <- nrow(site_env)
message("INFO: ", n_sites, " distinct demes/sites")

#=============================================================================
# 3. Coordinates -- derive a valid, well-spaced lat/lon grid from the actual
#    unique bio_1/bio_2 values per site (rank-based row/col), not a hardcoded
#    formula -- robust across Cases whose grid orientation may differ.
#=============================================================================
bio1_levels <- sort(unique(site_env$bio_1))
bio2_levels <- sort(unique(site_env$bio_2))
site_env[, row := match(bio_1, bio1_levels) - 1L]
site_env[, col := match(bio_2, bio2_levels) - 1L]
site_env[, latitude  := row * 2.0]
site_env[, longitude := col * 3.0]
message("INFO: coordinate grid: ", length(bio1_levels), " x ", length(bio2_levels),
        " (lat step 2.0 deg, lon step 3.0 deg)")

#=============================================================================
# 4. metadata.tsv -- site, sample, latitude, longitude, bio_1, bio_2
#    Row order follows the VCF sample order (positional alignment expected by
#    downstream pipeline steps).
#=============================================================================
metadata <- ind[, .(site, sample)] %>%
    dplyr::left_join(site_env[, .(site, latitude, longitude, bio_1, bio_2)], by = "site")

if (anyNA(metadata$latitude) || anyNA(metadata$longitude)) {
    stop("ERROR: NA coordinates after join -- check site_env completeness")
}

metadata_path <- file.path(OUT_DIR, "metadata.tsv")
fwrite(metadata, metadata_path, sep = "\t")
message("INFO: Wrote ", metadata_path, " (", nrow(metadata), " samples, ", n_sites, " sites)")

#=============================================================================
# 5. environments_present.tsv / environments_future.tsv
#    Future = present + fixed synthetic shift. No real "future"/garden concept
#    exists in this static reciprocal-transplant landscape -- this is a
#    placeholder to smoke-test the maladaptation machinery only (see dataset notes §4).
#=============================================================================
present_path <- file.path(OUT_DIR, "environments_present.tsv")
fwrite(site_env[, .(site, bio_1, bio_2)], present_path, sep = "\t")
message("INFO: Wrote ", present_path)

FUTURE_SHIFT <- 0.5
future_env <- copy(site_env[, .(site, bio_1, bio_2)])
future_env[, bio_1 := bio_1 + FUTURE_SHIFT]
future_env[, bio_2 := bio_2 + FUTURE_SHIFT]
future_path <- file.path(OUT_DIR, "environments_future.tsv")
fwrite(future_env, future_path, sep = "\t")
message("INFO: Wrote ", future_path, " (placeholder +", FUTURE_SHIFT,
        " synthetic shift, not a real common-garden target)")

#=============================================================================
# 6. causal_loci.tsv -- chr, pos, category (causal | linked_neutral | background_neutral)
#    No effect_size available in the source file (contra the original roadmap
#    assumption) -- see dataset notes §3.
#=============================================================================
causal_pos <- as.integer(readLines(causal_path))
message("INFO: ", length(causal_pos), " causal loci (VCF POS values, confirmed exact matches)")

vcf_pos_dt <- fread(vcf_out_path, skip = "#CHROM", header = TRUE,
                     select = 2, col.names = "pos")
all_pos <- vcf_pos_dt$pos

missing_causal <- setdiff(causal_pos, all_pos)
if (length(missing_causal) > 0) {
    stop("ERROR: ", length(missing_causal), " causal positions not found in VCF POS -- ",
         "e.g. ", missing_causal[1], ". Aborting rather than emitting a wrong truth table.")
}

is_causal <- all_pos %in% causal_pos
is_linked <- !is_causal & sapply(all_pos, function(p) any(abs(p - causal_pos) <= TP_WINDOW_BP))
category <- ifelse(is_causal, "causal", ifelse(is_linked, "linked_neutral", "background_neutral"))

causal_loci <- data.table(chr = CHR, pos = all_pos, category = category)
causal_loci_path <- file.path(OUT_DIR, "causal_loci.tsv")
fwrite(causal_loci, causal_loci_path, sep = "\t")
message("INFO: Wrote ", causal_loci_path, " -- ",
        sum(category == "causal"), " causal, ",
        sum(category == "linked_neutral"), " linked_neutral, ",
        sum(category == "background_neutral"), " background_neutral ",
        "(TP window +-", TP_WINDOW_KB, " kb)")

#=============================================================================
# 7. fitness_timeseries.tsv -- raw pass-through, archived for a future eval
#    harness. NOT consumed by this adapter-only pipeline run (see dataset notes §3).
#=============================================================================
cg_sum <- fread(cg_sum_path, header = TRUE)
fitness_ts_path <- file.path(OUT_DIR, "fitness_timeseries.tsv")
fwrite(cg_sum, fitness_ts_path, sep = "\t")
message("INFO: Wrote ", fitness_ts_path, " (", nrow(cg_sum),
        " generation snapshots, population-average sympatry/allopatry -- not a per-site table)")

#=============================================================================
# 8. laruson_minimal.gff3 -- one synthetic "gene" record per causal locus (+-500bp)
#    so gene-annotation/enrichment rules resolve without needing a real GFF.
#=============================================================================
gff_lines <- c(
    "##gff-version 3",
    sprintf("##sequence-region %s 1 %d", CHR, max(all_pos) + 1000)
)
gene_starts <- pmax(1L, causal_pos - 500L)
gene_ends   <- causal_pos + 500L
for (i in seq_along(causal_pos)) {
    attrs <- sprintf("ID=gene_%s_%d;Name=causal_%s_%d;description=causal_%s_%d;biotype=protein_coding",
                      CHR, causal_pos[i], CHR, causal_pos[i], CHR, causal_pos[i])
    gff_lines <- c(gff_lines, paste(
        CHR, ".", "gene", gene_starts[i], gene_ends[i], ".", "+", ".", attrs,
        sep = "\t"
    ))
}
gff_path <- file.path(OUT_DIR, "laruson_minimal.gff3")
writeLines(gff_lines, gff_path)
message("INFO: Wrote ", gff_path, " (", length(causal_pos), " synthetic gene records)")

message("INFO: convert_laruson.R complete -- outputs in ", OUT_DIR)
