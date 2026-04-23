#!/usr/bin/env python3
"""
make_ld_windows.py — Create fixed-size non-overlapping LD-decay windows for WZA.

Replaces plink --blocks with genome-wide tiling at the empirical LD half-decay scale,
following Wu et al. 2026 (Science adz0777): "the genome was partitioned into 16,917
non-overlapping LD blocks spanning ~7 kb (decay scale) each".

Usage:
    python3 make_ld_windows.py <vcf> <ld_decay_tsv> <window_source> <group> <scope>
                               <fallback_kb> <min_window_kb> <output>

Arguments:
    vcf            : filtered VCF (read header for ##contig lengths, scan for SNP positions)
    ld_decay_tsv   : Structure/tables/ld_decay_half_distances.tsv
    window_source  : 'ld_decay' or 'fixed'
    group          : group row to use from TSV (default 'All')
    scope          : 'genome_wide' (single window size) or 'per_chromosome'
    fallback_kb    : window size in kb used when TSV unavailable or window_source=fixed
    min_window_kb  : minimum allowed window size in kb (floor)
    output         : output path for plink-compatible blocks.det

Output TSV columns (plink --blocks schema):
    CHR  F  BP1  BP2  KB  NSNPS  SNPS
"""

import sys
import os
import gzip
import math
import pandas as pd
from collections import defaultdict


def parse_vcf_contigs(vcf_path):
    """Extract chromosome lengths from VCF ##contig header lines."""
    lengths = {}
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##contig="):
                # ##contig=<ID=1,length=30427671>
                parts = line.strip().lstrip("##contig=<").rstrip(">").split(",")
                chrom = None
                length = None
                for p in parts:
                    k, _, v = p.partition("=")
                    if k.strip() == "ID":
                        chrom = v.strip()
                    elif k.strip() == "length":
                        try:
                            length = int(v.strip())
                        except ValueError:
                            pass
                if chrom and length:
                    lengths[chrom] = length
    return lengths


def scan_vcf_positions(vcf_path):
    """Stream VCF to collect all (chr, pos) SNP positions and max pos per chr."""
    positions = defaultdict(list)
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t", 5)
            if len(cols) < 2:
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            positions[chrom].append(pos)
    return positions


def get_window_sizes(ld_decay_tsv, window_source, group, scope, fallback_bp, min_bp, chroms):
    """Return a dict {chrom: window_size_bp} based on the LD decay TSV."""
    if window_source == "fixed":
        sys.stderr.write(f"INFO: window_source=fixed — using fallback {fallback_bp} bp for all chromosomes\n")
        return {c: fallback_bp for c in chroms}

    if not os.path.exists(ld_decay_tsv):
        sys.stderr.write(f"WARN: ld_decay_tsv not found ({ld_decay_tsv}) — using fallback {fallback_bp} bp\n")
        return {c: fallback_bp for c in chroms}

    try:
        df = pd.read_csv(ld_decay_tsv, sep="\t", dtype={"scope": str, "group": str})
    except Exception as e:
        sys.stderr.write(f"WARN: Could not read ld_decay_tsv ({e}) — using fallback {fallback_bp} bp\n")
        return {c: fallback_bp for c in chroms}

    df = df[df["group"] == group]
    if df.empty:
        sys.stderr.write(f"WARN: No rows for group='{group}' in ld_decay_tsv — using fallback\n")
        return {c: fallback_bp for c in chroms}

    result = {}

    if scope == "genome_wide":
        row = df[df["scope"] == "genome_wide"]
        if row.empty or pd.isna(row.iloc[0]["half_decay_bp"]):
            sys.stderr.write(f"WARN: No genome_wide row for group='{group}' — using fallback\n")
            w = fallback_bp
        else:
            w = int(row.iloc[0]["half_decay_bp"])
            w = max(w, min_bp)
            sys.stderr.write(f"INFO: Using genome_wide half_decay = {w} bp (group='{group}')\n")
        for c in chroms:
            result[c] = w

    else:  # per_chromosome
        for c in chroms:
            row = df[df["scope"] == str(c)]
            if row.empty or pd.isna(row.iloc[0]["half_decay_bp"]):
                sys.stderr.write(f"WARN: No per_chromosome row for group='{group}', chr='{c}' — using fallback\n")
                w = fallback_bp
            else:
                w = int(row.iloc[0]["half_decay_bp"])
                w = max(w, min_bp)
            result[c] = w
        sys.stderr.write(f"INFO: Per-chromosome window sizes: {result}\n")

    return result


def tile_chromosome(chrom, chrom_len, window_bp, snp_positions):
    """Generate non-overlapping windows for one chromosome; skip empty windows."""
    windows = []
    pos = 1
    while pos <= chrom_len:
        bp1 = pos
        bp2 = min(pos + window_bp - 1, chrom_len)
        pos = bp2 + 1

        # Count SNPs in this window
        n_snps = sum(1 for p in snp_positions if bp1 <= p <= bp2)
        if n_snps == 0:
            continue

        kb = (bp2 - bp1 + 1) / 1000.0
        windows.append({
            "CHR": chrom,
            "F": 0,
            "BP1": bp1,
            "BP2": bp2,
            "KB": round(kb, 1),
            "NSNPS": n_snps,
            "SNPS": "NA",
        })
    return windows


def run(vcf, ld_decay_tsv, window_source, group, scope, fallback_kb, min_window_kb, output):
    fallback_bp = int(float(fallback_kb) * 1000)
    min_bp = int(float(min_window_kb) * 1000)

    sys.stderr.write("INFO: Parsing VCF for chromosome lengths and SNP positions\n")
    chrom_lengths = parse_vcf_contigs(vcf)
    snp_positions = scan_vcf_positions(vcf)

    # Merge: use contig-header lengths; fall back to max observed position per chr
    all_chroms = sorted(set(list(chrom_lengths.keys()) + list(snp_positions.keys())))
    if not all_chroms:
        sys.exit("ERROR: No chromosomes found in VCF")

    chrom_len_map = {}
    for c in all_chroms:
        if c in chrom_lengths:
            chrom_len_map[c] = chrom_lengths[c]
        elif c in snp_positions:
            chrom_len_map[c] = max(snp_positions[c])
            sys.stderr.write(f"WARN: Chr {c} has no ##contig length; using max position {chrom_len_map[c]}\n")

    window_sizes = get_window_sizes(ld_decay_tsv, window_source, group, scope, fallback_bp, min_bp, all_chroms)

    # Tile chromosomes
    all_windows = []
    total_snps_assigned = 0
    for chrom in all_chroms:
        chrom_len = chrom_len_map[chrom]
        window_bp = window_sizes.get(chrom, fallback_bp)
        positions = sorted(snp_positions.get(chrom, []))
        wins = tile_chromosome(chrom, chrom_len, window_bp, positions)
        all_windows.extend(wins)
        n_chr_snps = sum(w["NSNPS"] for w in wins)
        total_snps_assigned += n_chr_snps
        sys.stderr.write(
            f"INFO: Chr {chrom}: {len(wins)} windows of {window_bp} bp, {n_chr_snps} SNPs assigned\n"
        )

    total_vcf_snps = sum(len(v) for v in snp_positions.values())
    if total_snps_assigned != total_vcf_snps:
        sys.stderr.write(
            f"WARN: Coverage mismatch — {total_snps_assigned} SNPs assigned vs "
            f"{total_vcf_snps} in VCF. Likely some positions outside chromosome boundaries.\n"
        )
    else:
        sys.stderr.write(f"INFO: Coverage OK — all {total_vcf_snps} SNPs assigned to windows\n")

    if not all_windows:
        sys.exit("ERROR: No windows generated — check VCF and LD decay parameters")

    out_df = pd.DataFrame(all_windows, columns=["CHR", "F", "BP1", "BP2", "KB", "NSNPS", "SNPS"])
    out_df.to_csv(output, sep="\t", index=False)

    n_wins = len(out_df)
    mean_snps = out_df["NSNPS"].mean()
    sys.stderr.write(
        f"INFO: Written {n_wins} windows, mean {mean_snps:.1f} SNPs/window "
        f"to {output}\n"
    )


if __name__ == "__main__":
    if len(sys.argv) != 9:
        sys.exit(
            f"Usage: {sys.argv[0]} vcf ld_decay_tsv window_source group scope "
            f"fallback_kb min_window_kb output"
        )
    vcf, ld_decay_tsv, window_source, group, scope = sys.argv[1:6]
    fallback_kb, min_window_kb, output = sys.argv[6], sys.argv[7], sys.argv[8]
    run(vcf, ld_decay_tsv, window_source, group, scope, fallback_kb, min_window_kb, output)
