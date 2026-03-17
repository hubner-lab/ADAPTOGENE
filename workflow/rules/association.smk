#=============================================================================
# MODULE 4: ASSOCIATION (EMMAX/LFMM)
#=============================================================================

# Full dataset processing (non-LD pruned) - needed for LFMM
rule vcf_to_lfmm_full:
    """Convert full filtered VCF to LEA formats for association analysis."""
    input:  vcf = W['vcf_filt']
    output:
        geno = W['geno_full'],
        lfmm = W['lfmm_full'],
        vcfsnp = W['vcfsnp_full'],
        removed = W['removed_full']
    log: f"{LOGDIR}association/vcf_to_lfmm_full.log"
    shell: "Rscript /pipeline/scripts/vcf2lfmm.R {input.vcf} > {log} 2>&1"

rule snmf_full:
    """Run sNMF on full (non-LD pruned) dataset for imputation."""
    input:  geno = W['geno_full']
    output: W['snmf_full']
    params: ks = K_START, ke = K_END, ploidy = PLOIDY, rep = REPEAT, mode = SNMF_PROJECT_MODE
    log:    f"{LOGDIR}association/snmf_full.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/snmf.R {input.geno} \
            {params.ks} {params.ke} {params.ploidy} {params.rep} \
            {threads} {params.mode} > {log} 2>&1
        """

rule impute_full:
    """Impute missing genotypes in full dataset using SNMF."""
    input:  snmf = W['snmf_full'], lfmm = W['lfmm_full']
    output: W['lfmm_imp_full']
    params: k = K_BEST
    log:    f"{LOGDIR}association/impute_full.log"
    shell:
        """
        Rscript /pipeline/scripts/impute.R {input.snmf} {input.lfmm} {params.k} {output} > {log} 2>&1
        """

rule lfmm2vcf_full:
    """Convert imputed full LFMM to VCF format."""
    input:  vcf = W['vcf_filt'], lfmm_imp = W['lfmm_imp_full']
    output: W['vcf_imp_full']
    params: blocksize = 20000
    log:    f"{LOGDIR}association/lfmm2vcf_full.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm2vcf.R {input.lfmm_imp} {input.vcf} {params.blocksize} > {log} 2>&1
        """

# TPED/kinship for association (shared by EMMAX and GAPIT)
rule tped_assoc:
    """Convert filtered VCF to TPED/TFAM for association EMMAX."""
    input: vcf = W['vcf_filt']
    output: tped = W['assoc_tped'], tfam = W['assoc_tfam']
    params: prefix = f"{WORK_FILT}emmax/{VCF_BASE}"
    log: f"{LOGDIR}association/tped_assoc.log"
    shell:
        """
        plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
            --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
        awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
            {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
        """

rule kinship_assoc:
    """Compute BN kinship matrix for association analysis."""
    input: tped = W['assoc_tped']
    output: W['assoc_kinship']
    params: prefix = f"{WORK_FILT}emmax/{VCF_BASE}"
    log: f"{LOGDIR}association/kinship_assoc.log"
    shell:
        "/pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"

# EMMAX analysis
rule emmax_analysis:
    """Run EMMAX association analysis for all traits."""
    input:
        vcf = W['vcf_filt'],
        tped = W['assoc_tped'],
        kinship = W['assoc_kinship'],
        traits = O['climate_site_scaled'],
        covariates = W['pca_projections'],  # LEA PCA projections
        metadata = O['metadata']
    output: assoc_pvalues("EMMAX")
    params:
        k = K_BEST,
        predictors = PREDICTORS_SELECTED,
        inter_dir = INTER,
        tables_dir = f"{MOD_ASSOC}tables/EMMAX/",
        tped_prefix = f"{WORK_FILT}emmax/{VCF_BASE}"
    log: f"{LOGDIR}association/emmax_analysis.log"
    shell:
        """
        Rscript /pipeline/scripts/emmax.R \
            {input.vcf} {params.k} {input.traits} {input.covariates} \
            {params.predictors} {params.inter_dir} {input.metadata} \
            {params.tables_dir} {params.tped_prefix} {input.kinship} > {log} 2>&1
        """

# LFMM analysis
rule lfmm_analysis:
    """Run LFMM association analysis for all traits."""
    input:
        lfmm_ld = W['lfmm_imp'],       # LD-pruned imputed (for model training)
        lfmm_full = W['lfmm_imp_full'], # Full imputed (for testing)
        climate = O['climate_site_scaled'],
        vcfsnp = W['vcfsnp_full']
    output: assoc_pvalues("LFMM")
    params:
        k = K_BEST,
        predictors = PREDICTORS_SELECTED,
        tables_dir = f"{MOD_ASSOC}tables/LFMM/"
    log: f"{LOGDIR}association/lfmm_analysis.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm.R \
            {input.lfmm_ld} {input.lfmm_full} {input.climate} \
            {params.k} {params.predictors} {input.vcfsnp} \
            {params.tables_dir} > {log} 2>&1
        """

# VCF to GAPIT numeric (shared by GEA and GWAS GAPIT models)
if ASSOC_GAPIT_CONFIGS or PHENO_GAPIT_CONFIGS:

    rule vcf_to_numeric:
        """Convert imputed VCF to GAPIT numeric format (GD + GM)."""
        input: vcf = W['vcf_imp_full']
        output:
            gd = W['gapit_gd'],
            gm = W['gapit_gm']
        log: f"{LOGDIR}association/vcf_to_numeric.log"
        shell:
            "Rscript /pipeline/scripts/vcf_to_gapit_numeric.R {input.vcf} {output.gd} {output.gm} > {log} 2>&1"

# GAPIT GEA analysis (GLM, MLM, CMLM, ECMLM, SUPER, MLMM, FarmCPU, BLINK)
if ASSOC_GAPIT_CONFIGS:

    rule gapit_analysis:
        """Run GAPIT association analysis for all GAPIT models and traits."""
        input:
            gd = W['gapit_gd'],
            gm = W['gapit_gm'],
            traits = O['climate_site_scaled'],
            pca = W['pca_projections'],
            kinship = W['assoc_kinship'],
            metadata = O['metadata']
        output: [assoc_pvalues(model) for model in ASSOC_GAPIT_CONFIGS]
        params:
            k = K_BEST,
            models = ','.join(ASSOC_GAPIT_CONFIGS.keys()),
            workdir = W['gapit_work'],
            tables_dir = f"{MOD_ASSOC}tables/",
            predictors = PREDICTORS_SELECTED,
            native_outdir = f"{MOD_ASSOC}GAPIT/"
        log: f"{LOGDIR}association/gapit_analysis.log"
        shell:
            """
            Rscript /pipeline/scripts/gapit.R \
                {input.gd} {input.gm} {input.traits} {input.pca} \
                {input.kinship} {params.k} {params.models} \
                {params.workdir} {params.tables_dir} {params.predictors} \
                {input.metadata} {params.native_outdir} NULL > {log} 2>&1
            """

# Manhattan plots - simple version (runs early, no region dependency)
# Produces both PNG and SVG in a single run
rule manhattan_plot:
    """Generate simple Manhattan plot for a specific trait and method."""
    input: assoc = lambda wc: assoc_pvalues(wc.method)
    output:
        png = f"{MOD_ASSOC}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}.png",
        svg = f"{MOD_ASSOC}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}.svg"
    wildcard_constraints:
        method = ASSOC_METHOD_REGEX,
        trait = r"bio_\d+",
        adjust = r"\w+_[\d.]+"
    params:
        k = K_BEST,
        plot_dir = lambda wc: f"{MOD_ASSOC}plots/manhattan/{wc.method}/",
        regions = "NULL",  # No regions for simple plot
        selected_snps = "NULL",  # No selected SNPs for simple plot
        scattermore_threshold = SCATTERMORE_THRESHOLD
    log: f"{LOGDIR}association/manhattan_{{method}}_{{trait}}_{{adjust}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_manhattan.R \
            {input.assoc} {wildcards.adjust} {params.k} {wildcards.method} \
            {wildcards.trait} {params.plot_dir} {params.regions} {params.selected_snps} \
            {params.scattermore_threshold} > {log} 2>&1
        """

# Manhattan plots with regions highlighted (runs after regions are created)
# Produces both PNG and SVG in a single run
# For combined methods (Sum/Overlap/PairOverlap), shows all selected SNPs
# with different shapes: circle=current method, triangle=other method, diamond=both
rule manhattan_plot_regions:
    """Generate Manhattan plot with per-trait regions highlighted."""
    input:
        assoc = lambda wc: assoc_pvalues(wc.method),
        regions = O['regions_per_trait'],
        # All sigSNPs files from all methods for per-trait method attribution
        sigsnps = lambda wc: [assoc_sigsnps(method, adjust) for method, adjust in ASSOC_CONFIGS.items()]
    output:
        png = f"{MOD_ASSOC}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}_regions.png",
        svg = f"{MOD_ASSOC}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}_regions.svg"
    wildcard_constraints:
        method = ASSOC_METHOD_REGEX,
        trait = r"bio_\d+",
        adjust = r"\w+_[\d.]+"
    params:
        k = K_BEST,
        plot_dir = lambda wc: f"{MOD_ASSOC}plots/manhattan/{wc.method}/",
        sigsnps_str = lambda wc, input: ','.join(input.sigsnps),
        scattermore_threshold = SCATTERMORE_THRESHOLD
    log: f"{LOGDIR}association/manhattan_regions_{{method}}_{{trait}}_{{adjust}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_manhattan.R \
            {input.assoc} {wildcards.adjust} {params.k} {wildcards.method} \
            {wildcards.trait} {params.plot_dir} {input.regions} "{params.sigsnps_str}" \
            {params.scattermore_threshold} > {log} 2>&1
        """

# Find significant SNPs - one rule per method
rule find_sig_snps:
    """Find significant SNPs for a specific method."""
    input: assoc = lambda wc: assoc_pvalues(wc.method)
    output: f"{MOD_ASSOC}tables/{{method}}/{{method}}_pvalues_K{K_BEST}_sig_snps_{{adjust}}.tsv"
    wildcard_constraints:
        method = ASSOC_METHOD_REGEX,
        adjust = r"\w+_[\d.]+"
    params: snp_dist = SNP_DISTANCE
    log: f"{LOGDIR}association/find_sig_snps_{{method}}_{{adjust}}.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/find_sig_snps.R \
            {input.assoc} {wildcards.adjust} {params.snp_dist} \
            {wildcards.method} {threads} {output} > {log} 2>&1
        """

# Combine significant SNPs from different methods into Selected_SNPs table
rule combine_selected_snps:
    """Combine significant SNPs from all association methods."""
    input:
        sigsnps = lambda wc: [assoc_sigsnps(method, adjust) for method, adjust in ASSOC_CONFIGS.items()]
    output: O['selected_snps']
    params:
        sigsnps_str = lambda wc, input: ' '.join(input.sigsnps),
        method = SIGSNPS_METHOD,
        gap = SIGSNPS_GAP,
        predictors = PREDICTORS_SELECTED
    log: f"{LOGDIR}association/combine_selected_snps.log"
    shell:
        """
        Rscript /pipeline/scripts/combine_selected_snps.R \
            "{params.sigsnps_str}" {params.method} {params.gap} \
            {params.predictors} {output} > {log} 2>&1
        """

# Merge SNPs into regions based on distance
rule create_regions:
    """Merge nearby significant SNPs into regions.
    Produces two outputs: per-trait regions and combined climate regions."""
    input: selected_snps = O['selected_snps']
    output:
        per_trait = O['regions_per_trait'],
        combined = O['regions_combined']
    params: region_dist = REGION_DISTANCE
    log: f"{LOGDIR}association/create_regions.log"
    shell:
        """
        Rscript /pipeline/scripts/create_regions.R \
            {input.selected_snps} {params.region_dist} \
            {output.per_trait} {output.combined} > {log} 2>&1
        """

# Find genes around regions
rule find_genes_around_regions:
    """Find genes overlapping per-trait regions."""
    input:
        regions = O['regions_per_trait'],
        gff = W['gff_normalized'],
        vcfsnp = W['vcfsnp_full']
    output:
        genes = O['genes_per_region'],
        collapsed = O['genes_per_region_collapsed']
    params:
        feature = GFF_FEATURE,
        promoter_len = PROMOTER_LENGTH,
        top_regions = TOP_REGIONS
    log: f"{LOGDIR}association/find_genes_around_regions.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/find_genes_around_regions.R \
            {input.gff} {input.regions} {params.feature} \
            {params.promoter_len} {input.vcfsnp} {threads} {params.top_regions} \
            {output.genes} {output.collapsed} > {log} 2>&1
        """

# Find genes for combined climate regions (reference only, no enrichment)
rule find_genes_combined_regions:
    """Find genes overlapping combined climate regions (all traits) for reference."""
    input:
        regions = O['regions_combined'],
        gff = W['gff_normalized'],
        vcfsnp = W['vcfsnp_full']
    output:
        genes = O['genes_combined_regions']
    params:
        feature = GFF_FEATURE,
        promoter_len = PROMOTER_LENGTH,
        top_regions = TOP_REGIONS
    log: f"{LOGDIR}association/find_genes_combined_regions.log"
    threads: CPU
    shell:
        """
        # Use same script but only keep the redundant output (genes file)
        # Combined regions don't need collapsed output since we don't run enrichment on them
        TEMP_COLLAPSED=$(mktemp)
        Rscript /pipeline/scripts/find_genes_around_regions.R \
            {input.gff} {input.regions} {params.feature} \
            {params.promoter_len} {input.vcfsnp} {threads} {params.top_regions} \
            {output.genes} $TEMP_COLLAPSED > {log} 2>&1
        rm -f $TEMP_COLLAPSED
        """

# GO enrichment analysis (only if GO_FIELD is specified)
rule run_enrichment:
    """Run per-region GO enrichment analysis for genes around significant regions.

    Creates per-region enrichment files in trait-specific subdirectories:
    - intermediate/enrichment/{trait}/Region_{region_id}_enrichment.qs (enrichResult object)
    - tables/association/enrichment/{trait}/Region_{region_id}_enrichment.tsv (summary table)
    - tables/association/enrichment/{trait}/Enrichment_{trait}_summary.tsv

    Output files are determined dynamically based on regions in Genes_per_region.tsv.
    """
    input:
        genes = O['genes_per_region'],  # Use redundant table (one row per region-gene pair)
        gff = W['gff_normalized']
    output: W['enrichment_done']
    params:
        go_field = GO_FIELD,
        feature = GFF_FEATURE,
        tables_dir = f"{MOD_ASSOC}tables/enrichment/",
        intermediate_dir = f"{INTER}enrichment/association/"
    log: f"{LOGDIR}association/run_enrichment.log"
    shell:
        """
        # Run enrichment (creates per-region and per-trait summary files)
        Rscript /pipeline/scripts/run_enrichment.R \
            {input.genes} {input.gff} {params.go_field} {params.feature} \
            {params.tables_dir} {params.intermediate_dir} > {log} 2>&1

        # Touch done file to indicate completion
        touch {output}
        """

# Enrichment visualization plots
rule plot_enrichment:
    """Generate enrichment visualization plots for each region.

    Creates three plot types for each region with significant enrichment:
    - emapplot: GO term similarity network
    - cnetplot: Gene-concept network
    - dotplot: Standard enrichment dotplot

    Plots are organized by trait in plots/association/enrichment/{trait}/
    Reads .qs enrichment objects from intermediate/enrichment/
    """
    input:
        enrichment_done = W['enrichment_done'],
        genes = O['genes_per_region'],
        regions = O['regions_per_trait']
    output:
        W['enrichment_plots_done']
    params:
        intermediate_dir = f"{INTER}enrichment/association/",
        plots_dir = f"{MOD_ASSOC}plots/enrichment/",
        top_terms = ENRICHMENT_TOP_TERMS,
        width = ENRICHMENT_PLOT_WIDTH,
        height = ENRICHMENT_PLOT_HEIGHT,
        cnet_label = ENRICHMENT_CNET_LABEL,
        top_plot_regions = ENRICHMENT_TOP_PLOT_REGIONS
    log: f"{LOGDIR}association/plot_enrichment.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_enrichment.R \
            {params.intermediate_dir} {params.plots_dir} {params.top_terms} \
            {params.width} {params.height} {params.cnet_label} \
            {input.genes} {params.top_plot_regions} {input.regions} > {log} 2>&1

        # Touch done file
        touch {output}
        """

# Combined Manhattan plots (all traits and methods)
rule manhattan_combined:
    """Generate combined Manhattan plots showing all traits and methods together.
    Produces two versions: simple and with regions highlighted."""
    input:
        assoc_tables = [assoc_pvalues(method) for method in ASSOC_CONFIGS],
        regions = O['regions_combined']
    output:
        simple_png = O['manhattan_combined'],
        simple_svg = f"{MOD_ASSOC}plots/manhattan_combined_K{K_BEST}.svg",
        regions_png = O['manhattan_combined_regions'],
        regions_svg = f"{MOD_ASSOC}plots/manhattan_combined_K{K_BEST}_regions.svg"
    params:
        assoc_str = ','.join([
            f"{method}:{adjust}:{assoc_pvalues(method)}"
            for method, adjust in ASSOC_CONFIGS.items()
        ]),
        predictors = PREDICTORS_SELECTED,
        k = K_BEST,
        plot_dir = f"{MOD_ASSOC}plots/"
    log: f"{LOGDIR}association/manhattan_combined.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_manhattan_combined.R \
            "{params.assoc_str}" {params.predictors} {params.k} \
            {input.regions} {params.plot_dir} > {log} 2>&1
        """
