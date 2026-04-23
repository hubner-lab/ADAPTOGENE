#!/usr/bin/env python3
"""
wza_run.py — Weighted-Z Analysis: aggregate per-SNP p-values to LD window level.

Implements the canonical Booker et al. 2024 (Mol Ecol Resources) WZA formula:
    z_i    = Phi^{-1}(1 - p_i)        [one-sided normal quantile; p clipped to (1e-15, 1-1e-3)]
    w_i    = MAF_i * (1 - MAF_i)      [heterozygosity weight; canonical Booker formula]
    Z_blk  = sum(w_i * z_i) / sqrt(sum(w_i^2))
    p_blk  = P(Z > Z_blk) = norm.sf(Z_blk)

SNP-count null correction (Booker spline, applied when >= 100 blocks):
    A cubic spline of mean(Z) and sd(Z) vs. n_snps is fit across all blocks.
    Z_adj = (Z - mean_hat(n)) / sd_hat(n)
    p_adj = norm.sf(Z_adj)
    When < 100 blocks (spline unstable) or sd_hat < 0.1, Z_adj = Z (no correction).

Follows Wu et al. 2026 (Science adz0777): blocks are fixed-size LD-decay-scale windows
(from make_ld_windows.py), so every SNP is assigned to exactly one block — no SNPs are
silently dropped due to gaps between plink Gabriel blocks.

Usage:
    python3 wza_run.py <pvalues_tsv> <blocks_det> <maf_frq> <output_tsv> <maf_filter> <min_snps>

Inputs:
    pvalues_tsv : TSV with columns SNPID, chr, pos, <trait1>, <trait2>, ...
    blocks_det  : output of make_ld_windows.py (CHR F BP1 BP2 KB NSNPS SNPS)
    maf_frq     : qc_maf_pos.tsv (CHR, POS, MAF)
    output_tsv  : output path for block-level p-values
    maf_filter  : minimum MAF for SNP inclusion (float)
    min_snps    : minimum valid SNPs per block to compute WZA (int)

Output TSV columns:
    block_id, chr, start, end, n_snps, <trait>_Z, <trait>_Z_adj, <trait>_p, <trait>_nsnps, ...
    where <trait>_p is the spline-adjusted p-value used in all downstream steps.
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.interpolate import UnivariateSpline


P_CLIP_LOW  = 1e-15
P_CLIP_HIGH = 1 - 1e-3
MIN_BLOCKS_FOR_SPLINE = 100


def assign_snps_to_blocks(pvals_chr, blocks_chr):
    """Return block_id array aligned to pvals_chr rows via interval containment."""
    if blocks_chr.empty:
        return np.full(len(pvals_chr), None)
    blk = blocks_chr.sort_values("BP1")
    bp1 = blk["BP1"].values
    bp2 = blk["BP2"].values
    bids = blk["block_id"].values
    pos = pvals_chr["pos"].values
    idx = np.searchsorted(bp1, pos, side="right") - 1
    clipped = np.clip(idx, 0, len(bp2) - 1)
    in_block = (idx >= 0) & (pos <= bp2[clipped])
    result = np.where(in_block, bids[clipped], None)
    n_unassigned = np.sum(result == None)  # noqa: E711
    if n_unassigned > 0:
        sys.stderr.write(
            f"WARN: {n_unassigned} SNPs on chr {pvals_chr['chr'].iloc[0]} "
            f"not assigned to any window (outside chromosome boundaries?)\n"
        )
    return result


def wza_block(p_arr, maf_arr, min_snps):
    """Compute WZA Z-statistic (pre-correction) for one block×trait.

    Returns (Z, p_raw, n_valid). p here is the raw parametric p before spline correction.
    """
    valid = np.isfinite(p_arr) & np.isfinite(maf_arr)
    n = valid.sum()
    if n < min_snps:
        return np.nan, np.nan, int(n)
    p_v = np.clip(p_arr[valid], P_CLIP_LOW, P_CLIP_HIGH)
    w_v = maf_arr[valid] * (1.0 - maf_arr[valid])  # canonical Booker weight: p*(1-p)
    z_v = norm.ppf(1.0 - p_v)
    denom = np.sqrt(np.sum(w_v ** 2))
    if denom == 0:
        return np.nan, np.nan, int(n)
    Z = float(np.sum(w_v * z_v) / denom)
    p = float(norm.sf(Z))
    return Z, p, int(n)


def adjust_wza_with_spline(block_df, trait_cols):
    """Apply SNP-count null correction (Booker spline) to all trait Z columns.

    Fits cubic spline of mean(Z) and sd(Z) vs n_snps across blocks.
    Sets Z_adj = (Z - mean_hat(n)) / sd_hat(n), p_adj = norm.sf(Z_adj).
    When < MIN_BLOCKS_FOR_SPLINE or spline is unstable, returns Z_adj = Z unchanged.
    """
    for trait in trait_cols:
        z_col = f"{trait}_Z"
        zadj_col = f"{trait}_Z_adj"
        p_col = f"{trait}_p"
        nsnps_col = f"{trait}_nsnps"

        valid_mask = block_df[z_col].notna() & block_df[nsnps_col].notna()
        n_valid_blocks = valid_mask.sum()

        if n_valid_blocks < MIN_BLOCKS_FOR_SPLINE:
            sys.stderr.write(
                f"WARN: Only {n_valid_blocks} valid blocks for trait '{trait}' "
                f"(< {MIN_BLOCKS_FOR_SPLINE}) — skipping spline null correction\n"
            )
            block_df[zadj_col] = block_df[z_col]
            block_df[p_col] = block_df[p_col]  # raw p unchanged
            continue

        sub = block_df[valid_mask].copy()
        x = sub[nsnps_col].values.astype(float)
        z = sub[z_col].values.astype(float)

        try:
            # Fit spline of mean and sd vs SNP count
            # Use binned approach: group by unique n_snps for robust fitting
            unique_ns = np.sort(np.unique(x))
            bin_means = np.array([z[x == n].mean() for n in unique_ns])
            bin_sds   = np.array([z[x == n].std(ddof=1) if (x == n).sum() > 1 else np.nan
                                   for n in unique_ns])

            # Drop bins with NaN sd
            valid_bins = np.isfinite(bin_means) & np.isfinite(bin_sds)
            if valid_bins.sum() < 4:
                raise ValueError("Too few valid bins for spline fit")

            xs = unique_ns[valid_bins]
            mean_spline = UnivariateSpline(xs, bin_means[valid_bins], s=len(xs), k=min(3, len(xs)-1))
            sd_spline   = UnivariateSpline(xs, bin_sds[valid_bins],   s=len(xs), k=min(3, len(xs)-1))

            mean_hat = mean_spline(x)
            sd_hat   = sd_spline(x)
            sd_hat   = np.maximum(sd_hat, 0.1)  # floor to avoid division by near-zero

            z_adj = (z - mean_hat) / sd_hat
            p_adj = norm.sf(z_adj)

            block_df.loc[valid_mask, zadj_col] = z_adj
            block_df.loc[~valid_mask, zadj_col] = np.nan
            block_df.loc[valid_mask, p_col]    = p_adj
            # p for rows with NaN Z stays NaN (already set)

            sys.stderr.write(
                f"INFO: Spline null correction applied for '{trait}' "
                f"({n_valid_blocks} blocks)\n"
            )

        except Exception as e:
            sys.stderr.write(
                f"WARN: Spline correction failed for trait '{trait}' ({e}) — using raw Z\n"
            )
            block_df[zadj_col] = block_df[z_col]

    return block_df


def run_wza(pvalues_tsv, blocks_det, maf_frq, output_tsv, maf_filter, min_snps):
    sys.stderr.write("INFO: Reading input files\n")
    pvals = pd.read_csv(pvalues_tsv, sep="\t", dtype={"chr": str})
    blocks = pd.read_csv(blocks_det, sep="\t", dtype={"CHR": str})
    maf = pd.read_csv(maf_frq, sep="\t", dtype={"CHR": str})

    trait_cols = [c for c in pvals.columns if c not in {"SNPID", "chr", "pos"}]
    pvals["pos"] = pd.to_numeric(pvals["pos"], errors="coerce")
    maf["POS"] = pd.to_numeric(maf["POS"], errors="coerce")
    maf = maf[maf["MAF"] >= maf_filter].rename(columns={"MAF": "maf_val", "CHR": "m_chr", "POS": "m_pos"})

    pvals = pvals.merge(
        maf, left_on=["chr", "pos"], right_on=["m_chr", "m_pos"], how="inner"
    ).drop(columns=["m_chr", "m_pos"])

    if pvals.empty:
        sys.stderr.write("WARN: No SNPs remain after MAF filter\n")
        _write_empty(output_tsv, trait_cols)
        return

    blocks["block_id"] = (
        blocks["CHR"].astype(str) + ":" +
        blocks["BP1"].astype(str) + "-" +
        blocks["BP2"].astype(str)
    )

    assigned = []
    for chrom, grp in pvals.groupby("chr"):
        blk_chr = blocks[blocks["CHR"] == chrom]
        ids = assign_snps_to_blocks(grp, blk_chr)
        assigned.append(pd.Series(ids, index=grp.index))
    pvals["block_id"] = pd.concat(assigned)

    n_before = len(pvals)
    pvals = pvals.dropna(subset=["block_id"])
    n_dropped = n_before - len(pvals)
    if n_dropped > 0:
        sys.stderr.write(f"WARN: {n_dropped} SNPs not assigned to any window (outside boundaries)\n")

    if pvals.empty:
        sys.stderr.write("WARN: No SNPs assigned to any LD window\n")
        _write_empty(output_tsv, trait_cols)
        return

    block_meta = blocks.set_index("block_id")[["CHR", "BP1", "BP2"]]

    # Compute WZA per block (raw Z, before spline correction)
    rows = []
    for block_id, grp in pvals.groupby("block_id"):
        meta = block_meta.loc[block_id]
        row = {
            "block_id": block_id,
            "chr": meta["CHR"],
            "start": int(meta["BP1"]),
            "end": int(meta["BP2"]),
            "n_snps": len(grp),
        }
        grp_maf = grp["maf_val"].values
        for trait in trait_cols:
            p_arr = grp[trait].values.astype(float)
            Z, p, n = wza_block(p_arr, grp_maf, min_snps)
            row[f"{trait}_Z"] = Z
            row[f"{trait}_p"] = p       # will be overwritten by spline-adjusted p
            row[f"{trait}_nsnps"] = n
        rows.append(row)

    out_df = pd.DataFrame(rows)

    # Initialize Z_adj columns = Z (overwritten by spline when enough blocks)
    for trait in trait_cols:
        out_df[f"{trait}_Z_adj"] = out_df[f"{trait}_Z"]

    # Apply spline null correction per trait
    out_df = adjust_wza_with_spline(out_df, trait_cols)

    # Drop rows where ALL traits have NaN adjusted p (empty blocks)
    p_cols = [f"{t}_p" for t in trait_cols]
    out_df = out_df.dropna(subset=p_cols, how="all")

    # Reorder columns: block_id, chr, start, end, n_snps, then per-trait groups
    fixed = ["block_id", "chr", "start", "end", "n_snps"]
    trait_order = []
    for t in trait_cols:
        trait_order += [f"{t}_Z", f"{t}_Z_adj", f"{t}_p", f"{t}_nsnps"]
    out_df = out_df[fixed + trait_order]

    out_df.to_csv(output_tsv, sep="\t", index=False)
    sys.stderr.write(
        f"INFO: WZA complete — {len(out_df)} windows written to {output_tsv}\n"
    )


def _write_empty(path, trait_cols):
    cols = ["block_id", "chr", "start", "end", "n_snps"]
    for t in trait_cols:
        cols += [f"{t}_Z", f"{t}_Z_adj", f"{t}_p", f"{t}_nsnps"]
    pd.DataFrame(columns=cols).to_csv(path, sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit(
            f"Usage: {sys.argv[0]} pvalues_tsv blocks_det maf_frq output_tsv maf_filter min_snps"
        )
    pvalues_tsv, blocks_det, maf_frq, output_tsv = sys.argv[1:5]
    maf_filter = float(sys.argv[5])
    min_snps = int(sys.argv[6])
    run_wza(pvalues_tsv, blocks_det, maf_frq, output_tsv, maf_filter, min_snps)
