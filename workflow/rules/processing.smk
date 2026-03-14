#=============================================================================
# MODULE 1: VCF PROCESSING
#=============================================================================

rule extract_samples:
    """Extract sample IDs from metadata for VCF subsetting."""
    input:  samples = f"{INDIR}{SAMPLES}"
    output: W['samples_list']
    log:    f"{LOGDIR}processing/extract_samples.log"
    shell:  "tail -n +2 {input.samples} | awk '{{print 0, $2}}' > {output} 2> {log}"

rule calculate_sample_missing:
    """Calculate per-sample missing genotype rate and filter samples."""
    input:
        vcf = f"{INDIR}{VCF_RAW}",
        samples = W['samples_list']
    output:
        stats = W['samples_missing_stats'],
        filtered = W['samples_filtered'],
        removed = W['samples_removed']
    params:
        threshold = SAMPLE_MISS,
        prefix = f"{INTER}samples/sample_miss_tmp"
    log: f"{LOGDIR}processing/calculate_sample_missing.log"
    threads: CPU
    shell:
        """
        # Calculate per-sample missingness using plink
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# --keep {input.samples} \
            --missing --out {params.prefix} > {log} 2>&1

        # Create stats file with header
        echo -e "FID\\tIID\\tMISS_PHENO\\tN_MISS\\tN_GENO\\tF_MISS" > {output.stats}
        tail -n +2 {params.prefix}.imiss >> {output.stats}

        # Filter samples: keep those with F_MISS <= threshold
        awk -v thresh={params.threshold} 'NR>1 && $6 <= thresh {{print $1, $2}}' {params.prefix}.imiss > {output.filtered}

        # List removed samples: those with F_MISS > threshold
        awk -v thresh={params.threshold} 'NR>1 && $6 > thresh {{print $1, $2, $6}}' {params.prefix}.imiss > {output.removed}

        # Log summary
        n_total=$(wc -l < {input.samples})
        n_kept=$(wc -l < {output.filtered})
        n_removed=$(wc -l < {output.removed})
        echo "INFO: Sample missingness filtering (threshold: {params.threshold})" >> {log}
        echo "INFO: Total samples: $n_total" >> {log}
        echo "INFO: Samples passing: $n_kept" >> {log}
        echo "INFO: Samples removed: $n_removed" >> {log}

        # Cleanup temp files
        rm -f {params.prefix}.*
        """

rule filter_vcf:
    """Filter VCF by MAF, missingness, and sample list (after sample missingness filter).
    Also normalizes chromosome names by removing 'chr' prefix to match LEA behavior."""
    input:
        vcf = f"{INDIR}{VCF_RAW}",
        samples = W['samples_filtered']  # Use filtered samples list
    output: W['vcf_filt']
    params: prefix = W['vcf_filt'].replace('.vcf', ''), maf = MAF, miss = MISS
    log:    f"{LOGDIR}processing/filter_vcf.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# --keep {input.samples} \
            --maf {params.maf} --geno {params.miss} \
            --recode vcf --out {params.prefix} > {log} 2>&1
        sed -i '/^#CHROM/s/\\t0_/\\t/g' {output}

        # Normalize chromosome names: strip 'chr' prefix (e.g., chr1 -> 1, chr2H -> 2H)
        # This ensures consistency with LEA's vcf2lfmm behavior
        sed -i 's/^chr//g' {output}

        echo "INFO: Normalized chromosome names (stripped 'chr' prefix)" >> {log}
        """

rule ld_prune:
    """LD prune VCF (PCA is done separately via LEA)."""
    input:  vcf = W['vcf_filt']
    output: vcf = W['vcf_ld'], prune = W['prune_in']
    params: prefix = W['vcf_ld'].replace('.vcf', ''), win = LD_WIN, step = LD_STEP, r2 = LD_R2
    log:    f"{LOGDIR}processing/ld_prune.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# \
            --indep-pairwise {params.win} {params.step} {params.r2} \
            --out {params.prefix} > {log} 2>&1

        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# --extract {output.prune} \
            --make-bed --recode vcf \
            --out {params.prefix} >> {log} 2>&1

        sed -i '/^#CHROM/s/\\t0_0_/\\t/g' {output.vcf}
        """

rule extract_vcf_sample_order:
    """Get sample order from VCF header for metadata alignment."""
    input:  vcf = W['vcf_filt']
    output: W['samples_order']
    log:    f"{LOGDIR}processing/extract_vcf_sample_order.log"
    shell:  "grep -m1 CHROM {input.vcf} | cut -f10- | tr '\\t' '\\n' > {output} 2> {log}"

rule align_metadata:
    """Align metadata rows to match VCF sample order."""
    input:  meta = f"{INDIR}{SAMPLES}", order = W['samples_order']
    output: O['metadata']
    log:    f"{LOGDIR}processing/align_metadata.log"
    shell:  "Rscript /pipeline/scripts/filter_arrange_metadata.R {input.meta} {input.order} {output} > {log} 2>&1"

rule normalize_gff:
    """Normalize GFF chromosome names by removing 'chr' prefix to match VCF.
    Creates a normalized GFF in intermediate directory used by all downstream analysis."""
    input:  gff = f"{INDIR}{GFF}" if GFF else []
    output: W['gff_normalized']
    log:    f"{LOGDIR}processing/normalize_gff.log"
    shell:
        """
        if [ -f "{input.gff}" ]; then
            # Copy GFF and normalize chromosome names (strip 'chr' prefix)
            grep '^#' {input.gff} > {output} 2> {log}
            grep -v '^#' {input.gff} | sed 's/^chr//g' >> {output} 2>> {log}
            echo "INFO: Normalized GFF chromosome names (stripped 'chr' prefix)" >> {log}
        else
            # Create empty file if no GFF provided
            touch {output}
            echo "INFO: No GFF provided, created empty normalized GFF" >> {log}
        fi
        """

rule vcf_to_lfmm:
    """Convert VCF to LEA formats (geno, lfmm)."""
    input:  vcf = W['vcf_ld']
    output: geno = W['geno'], lfmm = W['lfmm'], vcfsnp = W['vcfsnp'], removed = W['removed']
    log:    f"{LOGDIR}processing/vcf_to_lfmm.log"
    shell:  "Rscript /pipeline/scripts/vcf2lfmm.R {input.vcf} > {log} 2>&1"
