#!/usr/bin/env python3
"""
wza_run.py — Weighted-Z Analysis: aggregate per-SNP p-values to LD block level.

Usage:
    python3 wza_run.py <pvalues_tsv> <blocks_det> <maf_frq> <output_tsv> <maf_filter> <min_snps>

Inputs:
    pvalues_tsv : TSV with columns SNPID, chr, pos, <trait1>, <trait2>, ...
    blocks_det  : plink --blocks output (CHR F BP1 BP2 KB NSNPS SNPS)
    maf_frq     : plink --freq output (CHR SNP A1 A2 MAF NCHROBS), SNP = CHR:POS
    output_tsv  : output path for block-level p-values
    maf_filter  : minimum MAF for SNP inclusion (float)
    min_snps    : minimum SNPs per block to compute WZA (int)

Output TSV columns:
    block_id, chr, start, end, n_snps, <trait>_Z, <trait>_p, <trait>_nsnps, ...

WZA formula (Booker et al. 2024, Mol Ecol Resources):
    z_i    = Phi^{-1}(1 - p_i)          [one-sided normal quantile]
    w_i    = sqrt(2 * MAF_i * (1-MAF_i)) [heterozygosity weight]
    Z_blk  = sum(w_i * z_i) / sqrt(sum(w_i^2))
    p_blk  = P(Z > Z_blk) = norm.sf(Z_blk)
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats import norm


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
    return np.where(in_block, bids[clipped], None)


def wza_block(p_arr, maf_arr, min_snps):
    """Compute WZA Z-statistic and one-sided p-value for one block×trait."""
    valid = np.isfinite(p_arr) & (p_arr > 0) & (p_arr < 1) & np.isfinite(maf_arr)
    n = valid.sum()
    if n < min_snps:
        return np.nan, np.nan, int(n)
    p_v = p_arr[valid]
    w_v = np.sqrt(2 * maf_arr[valid] * (1 - maf_arr[valid]))
    z_v = norm.ppf(1 - p_v)
    denom = np.sqrt(np.sum(w_v ** 2))
    if denom == 0:
        return np.nan, np.nan, int(n)
    Z = np.sum(w_v * z_v) / denom
    p = float(norm.sf(Z))
    return float(Z), p, int(n)


def run_wza(pvalues_tsv, blocks_det, maf_frq, output_tsv, maf_filter, min_snps):
    sys.stderr.write("INFO: Reading input files\n")
    pvals = pd.read_csv(pvalues_tsv, sep="\t", dtype={"chr": str})
    blocks = pd.read_csv(blocks_det, sep=r"\s+", dtype={"CHR": str})
    # maf_frq = qc_maf_pos.tsv: CHR (str), POS (int), MAF (float)
    maf = pd.read_csv(maf_frq, sep="\t", dtype={"CHR": str})

    trait_cols = [c for c in pvals.columns if c not in {"SNPID", "chr", "pos"}]
    pvals["pos"] = pd.to_numeric(pvals["pos"], errors="coerce")
    maf["POS"] = pd.to_numeric(maf["POS"], errors="coerce")
    maf = maf[maf["MAF"] >= maf_filter].rename(columns={"MAF": "maf_val", "CHR": "m_chr", "POS": "m_pos"})

    # Merge MAF by chr+pos (position-based join, independent of SNP ID format)
    pvals = pvals.merge(
        maf, left_on=["chr", "pos"], right_on=["m_chr", "m_pos"], how="inner"
    ).drop(columns=["m_chr", "m_pos"])

    if pvals.empty:
        sys.stderr.write("WARN: No SNPs remain after MAF filter\n")
        _write_empty(output_tsv, trait_cols)
        return

    # Build block_id (chr:start-end) and assign SNPs by position
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
    pvals = pvals.dropna(subset=["block_id"])

    if pvals.empty:
        sys.stderr.write("WARN: No SNPs assigned to any LD block\n")
        _write_empty(output_tsv, trait_cols)
        return

    # Block metadata lookup
    block_meta = blocks.set_index("block_id")[["CHR", "BP1", "BP2"]]

    # Compute WZA per block
    rows = []
    maf_vals = pvals["maf_val"].values
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
            row[f"{trait}_p"] = p
            row[f"{trait}_nsnps"] = n
        rows.append(row)

    out_df = pd.DataFrame(rows)
    out_df = out_df.dropna(subset=[f"{t}_p" for t in trait_cols], how="all")
    out_df.to_csv(output_tsv, sep="\t", index=False)
    sys.stderr.write(
        f"INFO: WZA complete — {len(out_df)} blocks written to {output_tsv}\n"
    )


def _write_empty(path, trait_cols):
    cols = ["block_id", "chr", "start", "end", "n_snps"]
    for t in trait_cols:
        cols += [f"{t}_Z", f"{t}_p", f"{t}_nsnps"]
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
