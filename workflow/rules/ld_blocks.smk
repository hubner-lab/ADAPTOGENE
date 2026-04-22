# ld_blocks.smk — LD block computation and WZA block-level association rules.
#
# Activated only when GEA.block_mode=block or GWAS.block_mode=block.
# All rules are guarded by _BLOCK_SOURCE_REGEX derived from active block-mode sources.
#
# Block workflow per source×method:
#   vcf_filt → plink --blocks → blocks.det
#   blocks.det + {method}_pvalues_K{k}.tsv + maf_filtered.frq → wza_run.py → {method}_block_pvalues_K{k}.tsv
#   → find_sig_blocks.R  → {method}_block_pvalues_K{k}_sig_blocks_{adjust}.tsv
#   → combine_selected_blocks.R → selected_blocks.tsv
#   → select_within_block_snps.R → selected_snps.tsv  (same path as snp-mode; replaces that rule)
#   → create_regions.R (blocks branch) → regions_per_trait.tsv, regions_combined.tsv

# _BLOCK_SOURCES, _BLOCK_SOURCE_REGEX, _ALL_BLOCK_METHODS_REGEX defined in common.smk

if _BLOCK_SOURCES:

    # ─── STEP 0: Compute LD blocks (once per project) ───────────────────────

    rule compute_ld_blocks:
        """Compute genome-wide LD blocks using plink Gabriel 2002 algorithm."""
        input:
            vcf = W['vcf_filt']
        output:
            det = W['blocks_det']
        params:
            prefix  = f"{INTER}ld_blocks/plink",
            max_kb  = BLOCK_BLOCKS_MAX_KB
        log: f"{LOGDIR}processing/compute_ld_blocks.log"
        shell:
            """
            plink --vcf {input.vcf} \
                  --blocks no-pheno-req \
                  --blocks-max-kb {params.max_kb} \
                  --set-missing-var-ids @:# \
                  --allow-extra-chr \
                  --out {params.prefix} > {log} 2>&1
            cp {params.prefix}.blocks.det {output.det}
            """

    # ─── STEP 1: WZA — aggregate per-SNP p-values to per-block ─────────────

    rule assoc_run_wza:
        """Run WZA for all traits in a (source, method) p-values table."""
        input:
            pvalues    = f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_pvalues_K{K_BEST}.tsv",
            blocks_det = W['blocks_det'],
            maf        = O['qc_maf_pos']
        output:
            f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_block_pvalues_K{K_BEST}.tsv"
        wildcard_constraints:
            source = _BLOCK_SOURCE_REGEX,
            method = _ALL_BLOCK_METHODS_REGEX
        params:
            maf_filter = BLOCK_WZA_MAF_FILTER,
            min_snps   = BLOCK_WZA_MIN_SNPS
        log: f"{LOGDIR}{{source}}/wza_{{method}}.log"
        shell:
            """
            python3 /pipeline/scripts/wza_run.py \
                {input.pvalues} {input.blocks_det} {input.maf} \
                {output} {params.maf_filter} {params.min_snps} > {log} 2>&1
            """

    # ─── STEP 2: Find significant blocks ────────────────────────────────────

    rule assoc_find_sig_blocks:
        """Apply multiple-testing adjustment at block level for a (source, method, adjust) combination."""
        input:
            block_pvalues = f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_block_pvalues_K{K_BEST}.tsv"
        output:
            f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_block_pvalues_K{K_BEST}_sig_blocks_{{adjust}}.tsv"
        wildcard_constraints:
            source = _BLOCK_SOURCE_REGEX,
            method = _ALL_BLOCK_METHODS_REGEX,
            adjust = r"\w+_[\d.]+"
        params:
            method_label = lambda wc: wc.method
        log: f"{LOGDIR}{{source}}/find_sig_blocks_{{method}}_{{adjust}}.log"
        shell:
            """
            Rscript /pipeline/scripts/find_sig_blocks.R \
                {input.block_pvalues} {wildcards.adjust} \
                {params.method_label} {output} > {log} 2>&1
            """

    # ─── STEP 3: Combine significant blocks across methods ──────────────────

    rule assoc_combine_selected_blocks:
        """Combine significant blocks from all methods for a given source."""
        input:
            sig_blocks = lambda wc: [
                f"{OUTDIR}{wc.source}/tables/methods/{m}/{m}_block_pvalues_K{K_BEST}_sig_blocks_{a}.tsv"
                for m, a in _src(wc.source, "configs").items()
            ]
        output:
            f"{OUTDIR}{{source}}/tables/selected_blocks.tsv"
        wildcard_constraints:
            source = _BLOCK_SOURCE_REGEX
        params:
            sigblocks_str = lambda wc, input: " ".join(input.sig_blocks),
            method        = lambda wc: _src(wc.source, "combine_method"),
            gap           = lambda wc: _src(wc.source, "combine_gap"),
            predictors    = lambda wc: _src(wc.source, "predictors")
        log: f"{LOGDIR}{{source}}/combine_selected_blocks.log"
        shell:
            """
            Rscript /pipeline/scripts/combine_selected_blocks.R \
                "{params.sigblocks_str}" {params.method} {params.gap} \
                {params.predictors} {output} > {log} 2>&1
            """

    # ─── STEP 4: Within-block SNP selection ─────────────────────────────────
    # Produces selected_snps.tsv at the same path as snp-mode, feeding Manhattan + GF.

    rule assoc_select_within_block_snps:
        """Select per-SNP significant entities within significant blocks."""
        input:
            selected_blocks = f"{OUTDIR}{{source}}/tables/selected_blocks.tsv",
            pvalues_files   = lambda wc: [
                f"{OUTDIR}{wc.source}/tables/methods/{m}/{m}_pvalues_K{K_BEST}.tsv"
                for m in _src(wc.source, "configs").keys()
            ]
        output:
            f"{OUTDIR}{{source}}/tables/selected_snps.tsv"
        wildcard_constraints:
            source = _BLOCK_SOURCE_REGEX
        params:
            pvalues_str  = lambda wc, input: " ".join(input.pvalues_files),
            within_method= BLOCK_WITHIN_METHOD,
            multiplier   = BLOCK_WITHIN_MULTIPLIER,
            threshold    = BLOCK_WITHIN_THRESHOLD,
            top_n        = BLOCK_WITHIN_TOP_N,
            predictors   = lambda wc: _src(wc.source, "predictors")
        log: f"{LOGDIR}{{source}}/select_within_block_snps.log"
        shell:
            """
            Rscript /pipeline/scripts/select_within_block_snps.R \
                {input.selected_blocks} \
                "{params.pvalues_str}" \
                {params.within_method} \
                {params.multiplier} \
                {params.threshold} \
                {params.top_n} \
                {params.predictors} \
                {output} > {log} 2>&1
            """

    # ─── STEP 5: Create regions from significant blocks ──────────────────────
    # Replaces assoc_create_regions for block-mode sources.
    # Reads selected_blocks.tsv and reformats into regions_per_trait + regions_combined
    # using the same schema as snp-mode (so all downstream scripts are unchanged).

    rule assoc_create_regions_from_blocks:
        """Format significant blocks as per-trait and combined region tables."""
        input:
            selected_blocks = f"{OUTDIR}{{source}}/tables/selected_blocks.tsv",
            selected_snps   = f"{OUTDIR}{{source}}/tables/selected_snps.tsv"
        output:
            per_trait = f"{OUTDIR}{{source}}/tables/regions_per_trait.tsv",
            combined  = f"{OUTDIR}{{source}}/tables/regions_combined.tsv"
        wildcard_constraints:
            source = _BLOCK_SOURCE_REGEX
        log: f"{LOGDIR}{{source}}/create_regions_from_blocks.log"
        shell:
            """
            Rscript /pipeline/scripts/create_regions_from_blocks.R \
                {input.selected_blocks} {input.selected_snps} \
                {output.per_trait} {output.combined} > {log} 2>&1
            """
