#!/usr/bin/env python3
"""Screen one Capblancq simulation run: architecture, truth set, structure band.

Emits one TSV row per run. The decisive columns are pc1_var / r2_pc_env:
`docs/gea_simulation_protocol.md` Part 1 requires structure that *partially*
masks the environment (R^2(PC~env) ~ 0.3-0.6). Capblancq's two gradients are
orthogonal by construction, so where this dataset actually sits is a measurement,
not an assumption.

Usage: capblancq_screen.py <run_dir> [maf_thresholds]
"""
import gzip
import pathlib
import sys

import numpy as np

MAF_THRESHOLDS = (0.01, 0.05)


def read_vcf(path):
    """Return (geno 0/1/2 int8 [n_var, n_samp], pos, mut_type).

    Genotype parsing is done on raw bytes rather than per-field string ops:
    the GT block is fixed-width ("0|1" + tab), so the two allele characters sit
    at a predictable stride and can be sliced out with numpy in one pass.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    pos, mt, gt_blocks = [], [], []
    with opener(path, "rb") as fh:
        for raw in fh:
            if raw.startswith(b"#"):
                continue
            f = raw.rstrip(b"\n").split(b"\t")
            pos.append(int(f[1]))
            t = 0
            for kv in f[7].split(b";"):
                if kv.startswith(b"MT="):
                    t = int(kv[3:].split(b",")[0])
                    break
            mt.append(t)
            gt_blocks.append(b"\t".join(f[9:]))

    n_var = len(gt_blocks)
    arr = np.frombuffer(b"\t".join(gt_blocks), dtype=np.uint8)
    # each sample occupies 4 bytes: allele, '|', allele, '\t' (last has no tab)
    stride = 4
    n_samp = (len(arr) // n_var + 1) // stride
    arr = np.frombuffer(
        b"".join(b + b"\t" for b in gt_blocks), dtype=np.uint8
    ).reshape(n_var, n_samp * stride)
    geno = ((arr[:, 0::stride] == ord("1")).astype(np.int8)
            + (arr[:, 2::stride] == ord("1")).astype(np.int8))
    return geno, np.asarray(pos), np.asarray(mt)


def pca_stats(geno_subset, E, n_pc=3):
    """PC variance fractions and R^2 of each PC on the environmental matrix."""
    G = geno_subset.astype(np.float64).T
    G -= G.mean(axis=0)
    sd = G.std(axis=0)
    sd[sd == 0] = 1.0
    G /= sd
    U, S, _ = np.linalg.svd(G, full_matrices=False)
    var_expl = S ** 2 / (S ** 2).sum()
    pcs = U[:, :n_pc] * S[:n_pc]
    return var_expl, [r2(pcs[:, i], E) for i in range(n_pc)]


def read_col(path):
    return np.loadtxt(path)


def r2(y, X):
    """R^2 of y regressed on columns of X (with intercept)."""
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    ss_tot = ((y - y.mean()) ** 2).sum()
    return 1.0 - (resid ** 2).sum() / ss_tot if ss_tot > 0 else np.nan


def main():
    run = pathlib.Path(sys.argv[1])
    alpha = run.parent.name
    label = f"{alpha}/{run.name}"

    geno, pos, mt = read_vcf(run / "genome_step1.vcf")
    n_var, n_samp = geno.shape

    freq = geno.mean(axis=1) / 2.0
    maf = np.minimum(freq, 1.0 - freq)

    is_m2, is_m3 = mt == 2, mt == 3
    is_causal = is_m2 | is_m3

    out = {
        "run": label, "alpha": alpha.replace("alpha", ""), "rep": run.name,
        "n_var": n_var, "n_samples": n_samp,
        "n_causal_m2": int(is_m2.sum()), "n_causal_m3": int(is_m3.sum()),
        "n_causal_total": int(is_causal.sum()), "n_neutral": int((mt == 1).sum()),
        "pos_min": int(pos.min()), "pos_max": int(pos.max()),
    }
    for thr in MAF_THRESHOLDS:
        keep = maf >= thr
        tag = str(thr).replace(".", "")
        out[f"n_var_maf{tag}"] = int(keep.sum())
        out[f"n_causal_maf{tag}"] = int((is_causal & keep).sum())
        out[f"n_causal_m2_maf{tag}"] = int((is_m2 & keep).sum())
        out[f"n_causal_m3_maf{tag}"] = int((is_m3 & keep).sum())
        denom = is_causal.sum()
        out[f"recall_ceiling_maf{tag}"] = round((is_causal & keep).sum() / denom, 4) if denom else np.nan
    out["causal_maf_median"] = round(float(np.median(maf[is_causal])), 5) if is_causal.any() else np.nan
    out["causal_maf_min"] = round(float(maf[is_causal].min()), 5) if is_causal.any() else np.nan

    env1 = read_col(run / "var1_step1")
    env2 = read_col(run / "var2_step1")
    env1f = read_col(run / "var1_pred")
    env2f = read_col(run / "var2_pred")
    fit_now = read_col(run / "fitness_step1")
    fit_fut = read_col(run / "fitness_pred")
    E = np.column_stack([env1, env2])

    # --- structure -----------------------------------------------------------
    # Two PCAs, deliberately. "all" is the conventional diagnostic, but here
    # causal loci are a large fraction of the MAF-filtered set, so an all-SNP
    # PCA is dominated by the adaptive signal itself and overstates confounding.
    # The neutral-only PCA (MT==1) is the honest measure of how much NEUTRAL
    # structure aligns with the environment - i.e. what structure correction
    # would actually remove.
    keep = maf >= 0.05
    var_all, r2_all = pca_stats(geno[keep], E)
    keep_neut = keep & (mt == 1)
    var_neut, r2_neut = pca_stats(geno[keep_neut], E)

    out["n_var_pca_all"] = int(keep.sum())
    out["n_var_pca_neutral"] = int(keep_neut.sum())
    for i in range(3):
        out[f"pc{i+1}_var"] = round(float(var_all[i]), 5)
    out["pc1_3_var"] = round(float(var_all[:3].sum()), 5)
    for i in range(3):
        out[f"r2_pc{i+1}_env"] = round(float(r2_all[i]), 4)
    for i in range(3):
        out[f"pc{i+1}_var_neutral"] = round(float(var_neut[i]), 5)
    for i in range(3):
        out[f"r2_pc{i+1}_env_neutral"] = round(float(r2_neut[i]), 4)
    out["r2_env1_env2"] = round(float(np.corrcoef(env1, env2)[0, 1] ** 2), 4)

    # --- offset-arm sanity --------------------------------------------------
    out["fit_now_mean"] = round(float(fit_now.mean()), 5)
    out["fit_fut_mean"] = round(float(fit_fut.mean()), 5)
    out["logfit_decline_mean"] = round(float(np.log(fit_fut + 1e-12).mean()
                                             - np.log(fit_now + 1e-12).mean()), 5)
    # null baseline: Euclidean climate transfer distance
    dist = np.sqrt((env1f - env1) ** 2 + (env2f - env2) ** 2)
    out["env_shift_mean"] = round(float(dist.mean()), 5)
    lf = np.log(fit_fut + 1e-12) - np.log(fit_now + 1e-12)
    out["r2_null_dist_vs_logfit"] = round(float(r2(lf, dist.reshape(-1, 1))), 4)

    print("\t".join(str(v) for v in out.values()))
    if len(sys.argv) > 2 and sys.argv[2] == "--header":
        print("\t".join(out.keys()), file=sys.stderr)


if __name__ == "__main__":
    main()
