#=============================================================================
# MODULE 1: VCF PROCESSING + QC
#=============================================================================

rule raw_vcf_stats:
    """Run bcftools stats on raw VCF; output raw stats file parsed later by R."""
    input:  vcf = f"{INDIR}{VCF_RAW}"
    output: O['qc_raw_summary']
    log:    f"{LOGDIR}processing/raw_vcf_stats.log"
    shell:  "bcftools stats {input.vcf} > {output} 2> {log}"

rule extract_samples:
    """Extract sample IDs from metadata for VCF subsetting."""
    input:  samples = f"{INDIR}{SAMPLES}"
    output: W['samples_list']
    log:    f"{LOGDIR}processing/extract_samples.log"
    shell:  "tail -n +2 {input.samples} | awk '{{print 0, $2}}' > {output} 2> {log}"

rule calculate_sample_missing:
    """Calculate per-sample and per-SNP missing genotype rates. Filter samples by sample_miss threshold."""
    input:
        vcf = f"{INDIR}{VCF_RAW}",
        samples = W['samples_list']
    output:
        stats    = W['samples_missing_stats'],
        filtered = W['samples_filtered'],
        removed  = W['samples_removed'],
        lmiss    = O['qc_snp_miss_raw']
    params:
        threshold = SAMPLE_MISS,
        prefix    = f"{INTER}samples/sample_miss_tmp"
    log: f"{LOGDIR}processing/calculate_sample_missing.log"
    threads: CPU
    shell:
        """
        # Calculate per-sample and per-SNP missingness using plink
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# --keep {input.samples} \
            --missing --out {params.prefix} > {log} 2>&1

        # Create stats file with header (convert plink space-separated to TSV)
        echo -e "FID\tIID\tMISS_PHENO\tN_MISS\tN_GENO\tF_MISS" > {output.stats}
        tail -n +2 {params.prefix}.imiss | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' >> {output.stats}

        # Preserve per-SNP missingness for QC (convert plink space-sep to TSV)
        awk 'NR==1 {{print "CHR\tSNP\tN_MISS\tN_GENO\tF_MISS"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
            {params.prefix}.lmiss > {output.lmiss}

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

        # Cleanup temp files (keep .imiss and .lmiss — moved above)
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh {params.prefix}.imiss {params.prefix}.lmiss
        """

rule compute_het_raw:
    """Compute per-sample heterozygosity on the raw VCF (filtered samples).
    Used for het_vs_missingness plot and optional het outlier removal."""
    input:
        vcf     = f"{INDIR}{VCF_RAW}",
        samples = W['samples_filtered']
    output: W['het_stats_raw']
    params: prefix = f"{INTER}qc/het_raw_tmp"
    log:    f"{LOGDIR}processing/compute_het_raw.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# --keep {input.samples} \
            --het --out {params.prefix} > {log} 2>&1
        awk 'NR==1 {{print "FID\tIID\tO_HOM\tE_HOM\tN_NM\tF"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' \
            {params.prefix}.het > {output}
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh {params.prefix}.het
        """

if HET_OUTLIER_SD is not None:
    rule het_outlier_filter:
        """Remove samples with heterozygosity beyond het_outlier_sd SDs from mean.
        Produces updated samples list used by filter_vcf."""
        input:
            het     = W['het_stats_raw'],
            samples = W['samples_filtered']
        output: W['samples_het_filtered']
        params: sd_thresh = HET_OUTLIER_SD
        log:    f"{LOGDIR}processing/het_outlier_filter.log"
        shell:
            """
            # Compute obs_het per sample (N_NM - O_HOM) / N_NM, then filter beyond N SDs
            awk -v sd_thresh={params.sd_thresh} '
            BEGIN {{ OFS="\t" }}
            NR == 1 {{ next }}
            {{
                fid=$1; iid=$2; o_hom=$3; e_hom=$4; n_nm=$5
                het = (n_nm - o_hom) / n_nm
                samples[NR] = fid " " iid
                het_vals[NR] = het
                sum += het; count++
            }}
            END {{
                mean = sum / count
                for (i in het_vals) sq_sum += (het_vals[i] - mean)^2
                sd = sqrt(sq_sum / count)
                lo = mean - sd_thresh * sd
                hi = mean + sd_thresh * sd
                n_removed = 0
                for (i in samples) {{
                    if (het_vals[i] >= lo && het_vals[i] <= hi)
                        print samples[i]
                    else
                        n_removed++
                }}
                print "INFO: Het outlier filter: mean=" mean " sd=" sd " lo=" lo " hi=" hi > "/dev/stderr"
                print "INFO: Removed " n_removed " samples as het outliers" > "/dev/stderr"
            }}' {input.het} > {output} 2>> {log}
            """

rule compute_snp_freq_raw:
    """Compute allele frequencies on the raw VCF (filtered samples, before MAF filter)."""
    input:
        vcf     = f"{INDIR}{VCF_RAW}",
        samples = W['samples_het_filtered'] if HET_OUTLIER_SD is not None else W['samples_filtered']
    output: O['qc_maf_raw']
    params: prefix = f"{INTER}qc/maf_raw_tmp"
    log:    f"{LOGDIR}processing/compute_snp_freq_raw.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# --keep {input.samples} \
            --freq --out {params.prefix} > {log} 2>&1
        awk 'NR==1 {{print "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' \
            {params.prefix}.frq > {output}
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh {params.prefix}.frq
        """

rule compute_snp_density_raw:
    """Compute SNP density in 1Mb windows on the raw VCF."""
    input: f"{INDIR}{VCF_RAW}"
    output: O['qc_snp_density_raw']
    params: prefix = f"{INTER}qc/snp_density_raw_tmp"
    log:    f"{LOGDIR}processing/compute_snp_density_raw.log"
    shell:
        """
        /app/vcftools-0.1.16/bin/vcftools --vcf {input} --SNPdensity 1000000 --out {params.prefix} > {log} 2>&1
        cp {params.prefix}.snpden {output}
        rm -f {params.prefix}.snpden {params.prefix}.log
        """

if (MIN_DEPTH is not None or MAX_DEPTH is not None) and HAS_FORMAT_DP:
    _depth_expr_parts = []
    if MIN_DEPTH is not None:
        _depth_expr_parts.append(f"FMT/DP < {MIN_DEPTH}")
    if MAX_DEPTH is not None:
        _depth_expr_parts.append(f"FMT/DP > {MAX_DEPTH}")
    _DEPTH_EXPR = " | ".join(_depth_expr_parts)

    rule depth_filter:
        """Mask genotypes with FORMAT/DP outside the configured min/max range.
        Sets out-of-range genotypes to missing (.) before MAF/missingness filtering."""
        input:  vcf = f"{INDIR}{VCF_RAW}"
        output: W['vcf_depth_filtered']
        params: expr = _DEPTH_EXPR
        log:    f"{LOGDIR}processing/depth_filter.log"
        shell:
            """
            bcftools filter -S . -e '{params.expr}' {input.vcf} -o {output} > {log} 2>&1
            echo "INFO: Depth filter applied: {params.expr}" >> {log}
            n_masked=$(bcftools stats {output} | awk '/^SN.*missing genotypes/ {{print $4}}')
            echo "INFO: Total missing genotypes after depth masking: $n_masked" >> {log}
            """

rule filter_vcf:
    """Filter VCF by MAF, missingness, and sample list (after sample missingness filter).
    Also normalizes chromosome names by removing 'chr' prefix to match LEA behavior."""
    input:
        vcf = W['vcf_depth_filtered'] if ((MIN_DEPTH is not None or MAX_DEPTH is not None) and HAS_FORMAT_DP) else f"{INDIR}{VCF_RAW}",
        samples = W['samples_het_filtered'] if HET_OUTLIER_SD is not None else W['samples_filtered']
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
        sed -i '/^#CHROM/s/\t0_/\t/g' {output}

        # Normalize chromosome names: strip 'chr' prefix (e.g., chr1 -> 1, chr2H -> 2H)
        # This ensures consistency with LEA's vcf2lfmm behavior
        sed -i 's/^chr//g' {output}

        echo "INFO: Normalized chromosome names (stripped 'chr' prefix)" >> {log}
        """

rule compute_sample_het:
    """Compute per-sample heterozygosity on the filtered VCF (for het_vs_missingness plot)."""
    input:  vcf = W['vcf_filt']
    output: W['het_stats_filtered']
    params: prefix = f"{INTER}qc/het_filtered_tmp"
    log:    f"{LOGDIR}processing/compute_sample_het.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# \
            --het --out {params.prefix} > {log} 2>&1
        awk 'NR==1 {{print "FID\tIID\tO_HOM\tE_HOM\tN_NM\tF"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' \
            {params.prefix}.het > {output}
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh {params.prefix}.het
        """

rule compute_snp_freq_filtered:
    """Compute allele frequencies on the filtered VCF (post MAF+missingness filter)."""
    input:  vcf = W['vcf_filt']
    output:
        frq = O['qc_maf_filtered'],
        pos = O['qc_maf_pos']
    params: prefix = f"{INTER}qc/maf_filtered_tmp"
    log:    f"{LOGDIR}processing/compute_snp_freq_filtered.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# \
            --freq --make-bed --out {params.prefix} > {log} 2>&1
        awk 'NR==1 {{print "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' \
            {params.prefix}.frq > {output.frq}
        paste <(awk '{{print $1"\t"$4}}' {params.prefix}.bim) \
              <(awk 'NR>1 {{print $5}}' {params.prefix}.frq) | \
            awk 'BEGIN{{print "CHR\tPOS\tMAF"}} {{print}}' > {output.pos}
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh \
              {params.prefix}.frq {params.prefix}.bim {params.prefix}.bed {params.prefix}.fam
        """

rule compute_snp_density_filtered:
    """Compute SNP density in 1Mb windows on the filtered VCF."""
    input:  vcf = W['vcf_filt']
    output: O['qc_snp_density_filt']
    params: prefix = f"{INTER}qc/snp_density_filtered_tmp"
    log:    f"{LOGDIR}processing/compute_snp_density_filtered.log"
    shell:
        """
        /app/vcftools-0.1.16/bin/vcftools --vcf {input.vcf} --SNPdensity 1000000 --out {params.prefix} > {log} 2>&1
        cp {params.prefix}.snpden {output}
        rm -f {params.prefix}.snpden {params.prefix}.log
        """

if HAS_FORMAT_DP:
    rule compute_depth_stats:
        """Compute per-sample and per-site mean depth (only when FORMAT/DP field is present)."""
        input:  vcf = W['vcf_filt']
        output:
            sample = O['qc_depth_sample'],
            site   = O['qc_depth_site']
        params: prefix = f"{INTER}qc/depth_tmp"
        log:    f"{LOGDIR}processing/compute_depth_stats.log"
        shell:
            """
            /app/vcftools-0.1.16/bin/vcftools --vcf {input.vcf} --depth --out {params.prefix} > {log} 2>&1
            /app/vcftools-0.1.16/bin/vcftools --vcf {input.vcf} --site-mean-depth --out {params.prefix} >> {log} 2>&1
            cp {params.prefix}.idepth {output.sample}
            cp {params.prefix}.ldepth.mean {output.site}
            rm -f {params.prefix}.idepth {params.prefix}.ldepth.mean {params.prefix}.log
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

        sed -i '/^#CHROM/s/\t0_0_/\t/g' {output.vcf}
        """

rule extract_vcf_sample_order:
    """Get sample order from VCF header for metadata alignment."""
    input:  vcf = W['vcf_filt']
    output: W['samples_order']
    log:    f"{LOGDIR}processing/extract_vcf_sample_order.log"
    shell:  "grep -m1 CHROM {input.vcf} | cut -f10- | tr '\t' '\n' > {output} 2> {log}"

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

rule plot_qc_processing:
    """Generate all processing QC plots and derived tables."""
    input:
        raw_summary      = O['qc_raw_summary'],
        sample_miss      = W['samples_missing_stats'],
        snp_miss_raw     = O['qc_snp_miss_raw'],
        maf_raw          = O['qc_maf_raw'],
        maf_filtered     = O['qc_maf_filtered'],
        het_raw          = W['het_stats_raw'],
        het_filtered     = W['het_stats_filtered'],
        snp_density_raw  = O['qc_snp_density_raw'],
        snp_density_filt = O['qc_snp_density_filt'],
        metadata         = O['metadata'],
        vcf_filt         = W['vcf_filt'],
        vcf_ld           = W['vcf_ld'],
        depth_sample     = O['qc_depth_sample'] if HAS_FORMAT_DP else [],
        depth_site       = O['qc_depth_site'] if HAS_FORMAT_DP else []
    output:
        plot_sample_miss = O['qc_plot_sample_miss'],
        plot_het_miss    = O['qc_plot_het_miss'],
        plot_maf         = O['qc_plot_maf'],
        plot_snp_miss    = O['qc_plot_snp_miss'],
        plot_attrition   = O['qc_plot_attrition'],
        plot_snp_density = O['qc_plot_snp_density'],
        filtering_summary= O['qc_filtering_summary'],
        sample_het       = O['qc_sample_het'],
        **({"plot_depth": O['qc_plot_depth']} if HAS_FORMAT_DP else {})
    params:
        sample_miss_thresh = SAMPLE_MISS,
        snp_miss_thresh    = MISS,
        maf_thresh         = MAF,
        het_outlier_sd     = HET_OUTLIER_SD if HET_OUTLIER_SD is not None else 'NULL',
        has_dp             = 'TRUE' if HAS_FORMAT_DP else 'FALSE',
        depth_sample       = O['qc_depth_sample'] if HAS_FORMAT_DP else 'NULL',
        depth_site         = O['qc_depth_site'] if HAS_FORMAT_DP else 'NULL',
        plots_dir          = f"{MOD_PROCESSING}plots/",
        tables_dir         = f"{MOD_PROCESSING}tables/"
    log: f"{LOGDIR}processing/plot_qc_processing.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_qc_processing.R \
            {input.raw_summary} \
            {input.sample_miss} \
            {input.snp_miss_raw} \
            {input.maf_raw} \
            {input.maf_filtered} \
            {input.het_raw} \
            {input.het_filtered} \
            {input.snp_density_raw} \
            {input.snp_density_filt} \
            {input.metadata} \
            {input.vcf_filt} \
            {input.vcf_ld} \
            {params.sample_miss_thresh} \
            {params.snp_miss_thresh} \
            {params.maf_thresh} \
            {params.het_outlier_sd} \
            {params.has_dp} \
            {params.depth_sample} \
            {params.depth_site} \
            {params.plots_dir} \
            {params.tables_dir} > {log} 2>&1
        """
