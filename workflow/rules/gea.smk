#=============================================================================
# MODULE 4: ASSOCIATION (GEA)
# Non-GAPIT methods are declared dynamically from GEA_METHODS registry.
# GAPIT models are declared as a single multi-output rule (gapit_analysis).
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
    log: f"{LOGDIR}gea/vcf_to_lfmm_full.log"
    shell: "Rscript /pipeline/scripts/vcf2lfmm.R {input.vcf} > {log} 2>&1"

rule snmf_full:
    """Run sNMF on full (non-LD pruned) dataset for imputation."""
    input:  geno = W['geno_full']
    output: W['snmf_full']
    params: ks = K_BEST, ke = K_BEST, ploidy = PLOIDY, rep = REPEAT, mode = SNMF_PROJECT_MODE
    log:    f"{LOGDIR}gea/snmf_full.log"
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
    log:    f"{LOGDIR}gea/impute_full.log"
    shell:
        """
        Rscript /pipeline/scripts/impute.R {input.snmf} {input.lfmm} {params.k} {output} > {log} 2>&1
        """

rule lfmm2vcf_full:
    """Convert imputed full LFMM to VCF format."""
    input:  vcf = W['vcf_filt'], lfmm_imp = W['lfmm_imp_full']
    output: W['vcf_imp_full']
    params: blocksize = 20000
    log:    f"{LOGDIR}gea/lfmm2vcf_full.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm2vcf.R {input.lfmm_imp} {input.vcf} {params.blocksize} > {log} 2>&1
        """

# TPED/kinship for association (shared by EMMAX and GAPIT)
rule tped_gea:
    """Convert filtered VCF to TPED/TFAM for association EMMAX."""
    input: vcf = W['vcf_filt']
    output: tped = W['assoc_tped'], tfam = W['assoc_tfam']
    params: prefix = f"{WORK_FILT}emmax/{VCF_BASE}"
    log: f"{LOGDIR}gea/tped_assoc.log"
    shell:
        """
        plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
            --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
        awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
            {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
        """

rule kinship_gea:
    """Compute BN kinship matrix for association analysis."""
    input: tped = W['assoc_tped']
    output: W['assoc_kinship']
    params: prefix = f"{WORK_FILT}emmax/{VCF_BASE}"
    log: f"{LOGDIR}gea/kinship_assoc.log"
    shell:
        "/pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"

# --- EMMAX climate-coord subset (climate coordinate NA handling) ---
# EMMAX binds climate values to genotypes positionally (cbind on TFAM row order),
# so when samples with missing lat/lon are dropped from the climate table
# (filter_coord_samples), the genotype-side TPED/kinship must be rebuilt on the
# same coord-valid sample set. Mirrors subset_vcf_gwas/tped_gwas_trait exactly.
if "EMMAX" in GEA_OTHER_CONFIGS:

    rule subset_vcf_gea_climate:
        """Subset filtered VCF to coord-valid samples for GEA EMMAX."""
        input:
            vcf = W['vcf_filt'],
            samples = W['coord_valid_samples']
        output: W['vcf_filt_climate']
        params: prefix = f"{WORK_FILT}climate/{VCF_BASE}"
        log: f"{LOGDIR}gea/subset_vcf_gea_climate.log"
        shell:
            """
            plink --vcf {input.vcf} --keep {input.samples} --const-fid 0 --allow-extra-chr \
                --recode vcf --out {params.prefix} > {log} 2>&1
            sed -i '/^#CHROM/s/\\t0_/\\t/g' {output}
            """

    rule tped_gea_climate:
        """Convert coord-valid-subset VCF to TPED/TFAM for association EMMAX."""
        input: vcf = W['vcf_filt_climate']
        output: tped = W['assoc_tped_climate'], tfam = W['assoc_tfam_climate']
        params: prefix = f"{WORK_FILT}climate/emmax/{VCF_BASE}"
        log: f"{LOGDIR}gea/tped_gea_climate.log"
        shell:
            """
            plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
                --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
            awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
                {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
            """

    rule kinship_gea_climate:
        """Compute BN kinship matrix for the coord-valid subset."""
        input: tped = W['assoc_tped_climate']
        output: W['assoc_kinship_climate']
        params: prefix = f"{WORK_FILT}climate/emmax/{VCF_BASE}"
        log: f"{LOGDIR}gea/kinship_gea_climate.log"
        shell:
            "/pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"

# =============================================================================
# GEA Association rules — per-factor caching (Phase 2)
#
# Each bioclimatic factor is computed independently and persisted as:
#   _intermediate/gea_per_trait/{method}/{trait}_pvalues_K{k}.tsv
#
# The wide table consumed by downstream rules is assembled from the selected
# factors. Changing the selection only reruns the assembly + any NEW factors;
# factors already computed are reused. Changing K or filter params changes the
# imputation input path, cascading invalidation across all factors.
# =============================================================================

# --- EMMAX per-trait ---
if "EMMAX" in GEA_OTHER_CONFIGS:

    rule assoc_emmax_gea_trait:
        """Run EMMAX for a single GEA bioclimatic trait (per-factor caching).
        Uses the coord-valid-subset VCF/TPED/kinship (see subset_vcf_gea_climate) since
        climate values (traits) exclude samples with missing lat/lon; PCA covariates
        stay on the full cohort and are row-subset internally by load_pca_covariates()
        via samples_order (see emmax.R)."""
        input:
            vcf        = W['vcf_filt_climate'],
            tped       = W['assoc_tped_climate'],
            kinship    = W['assoc_kinship_climate'],
            traits     = O['climate_site_scaled'],
            covariates = W['pca_projections'],
            metadata   = W['metadata_climate'],
            samples_order = W['samples_order'],
        output: f"{INTER}gea_per_trait/EMMAX/{{trait}}_pvalues_K{K_BEST}.tsv"
        wildcard_constraints:
            trait = r"bio_\d+"
        params:
            k           = K_BEST,
            inter_dir   = INTER,
            tped_prefix = f"{WORK_FILT}climate/emmax/{VCF_BASE}",
        log: f"{LOGDIR}gea/assoc_emmax_{{trait}}.log"
        shell:
            """
            Rscript /pipeline/scripts/emmax.R \
                {input.vcf} {params.k} {input.traits} {input.covariates} \
                {wildcards.trait} {params.inter_dir} {input.metadata} \
                {output} {params.tped_prefix} {input.kinship} {input.samples_order} > {log} 2>&1
            """

    rule assoc_emmax_gea_assemble:
        """Assemble per-trait EMMAX files into the wide table consumed by downstream."""
        input: [assoc_pvalues_trait("EMMAX", t) for t in get_predictors_list()]
        output: assoc_pvalues("EMMAX")
        params:
            traits_str = lambda wc, input: " ".join(input)
        log: f"{LOGDIR}gea/assoc_emmax_assemble.log"
        shell:
            "Rscript /pipeline/scripts/assemble_pvalues.R {params.traits_str} {output} > {log} 2>&1"

# --- LFMM per-trait ---
if "LFMM" in GEA_OTHER_CONFIGS:

    rule assoc_lfmm_gea_trait:
        """Run LFMM for a single GEA bioclimatic trait (per-factor caching).
        Uses coord-valid-subset genotype matrices (see subset_lfmm_climate in
        processing.smk) since climate (traits) excludes samples with missing lat/lon
        and LFMM binds climate rows to genotype-matrix rows positionally."""
        input:
            lfmm_ld  = W['lfmm_imp_climate'],
            lfmm_full = W['lfmm_imp_full_climate'],
            climate  = O['climate_site_scaled'],
            vcfsnp   = W['vcfsnp_full'],
        output: f"{INTER}gea_per_trait/LFMM/{{trait}}_pvalues_K{K_BEST}.tsv"
        wildcard_constraints:
            trait = r"bio_\d+"
        params:
            k = K_BEST,
        log: f"{LOGDIR}gea/assoc_lfmm_{{trait}}.log"
        shell:
            """
            Rscript /pipeline/scripts/lfmm.R \
                {input.lfmm_ld} {input.lfmm_full} {input.climate} \
                {params.k} {wildcards.trait} {input.vcfsnp} \
                {output} > {log} 2>&1
            """

    rule assoc_lfmm_gea_assemble:
        """Assemble per-trait LFMM files into the wide table consumed by downstream."""
        input: [assoc_pvalues_trait("LFMM", t) for t in get_predictors_list()]
        output: assoc_pvalues("LFMM")
        params:
            traits_str = lambda wc, input: " ".join(input)
        log: f"{LOGDIR}gea/assoc_lfmm_assemble.log"
        shell:
            "Rscript /pipeline/scripts/assemble_pvalues.R {params.traits_str} {output} > {log} 2>&1"

# VCF to GAPIT numeric (shared by GEA and GWAS GAPIT models)
if GEA_GAPIT_CONFIGS or GWAS_GAPIT_CONFIGS:

    rule vcf_to_numeric:
        """Convert imputed VCF to GAPIT numeric format (GD + GM)."""
        input: vcf = W['vcf_imp_full']
        output:
            gd = W['gapit_gd'],
            gm = W['gapit_gm']
        log: f"{LOGDIR}gea/vcf_to_numeric.log"
        shell:
            "Rscript /pipeline/scripts/vcf_to_gapit_numeric.R {input.vcf} {output.gd} {output.gm} > {log} 2>&1"

# --- GAPIT GEA per-trait ---
if GEA_GAPIT_CONFIGS:

    rule gapit_gea_trait:
        """Run GAPIT for a single GEA bioclimatic trait, all configured models in one call.
        Output files land under _intermediate/gea_per_trait/{model}/{trait}_pvalues_K{k}.tsv.
        GAPIT writes {TABLES_DIR}/{model}/{OUTPUT_PREFIX}_pvalues_K{k}.tsv (existing behaviour
        when OUTPUT_PREFIX is set), so we reuse gapit.R unchanged.
        Uses full-dataset GD/kinship (climate coordinate NA handling); gapit.R subsets
        internally via coord_samples, same mechanism as GWAS DROP mode."""
        input:
            gd       = W['gapit_gd'],
            gm       = W['gapit_gm'],
            traits   = O['climate_site_scaled'],
            pca      = W['pca_projections'],
            kinship  = W['assoc_kinship'],
            metadata = O['metadata'],
            coord_samples = W['coord_valid_samples'],
        output:
            [f"{INTER}gea_per_trait/{model}/{{trait}}_pvalues_K{K_BEST}.tsv"
             for model in GEA_GAPIT_CONFIGS]
        wildcard_constraints:
            trait = r"bio_\d+"
        params:
            k             = K_BEST,
            models        = ','.join(GEA_GAPIT_CONFIGS.keys()),
            workdir       = lambda wc: f"{INTER}gapit/gea/{wc.trait}/",
            tables_dir    = f"{INTER}gea_per_trait/",
            native_outdir = f"{MOD_GEA}GAPIT_native_output/",
        log: f"{LOGDIR}gea/gapit_gea_{{trait}}.log"
        shell:
            """
            Rscript /pipeline/scripts/gapit.R \
                {input.gd} {input.gm} {input.traits} {input.pca} \
                {input.kinship} {params.k} {params.models} \
                {params.workdir} {params.tables_dir} {wildcards.trait} \
                {input.metadata} {params.native_outdir} {input.coord_samples} \
                {wildcards.trait} > {log} 2>&1
            """

    for _gapit_model in GEA_GAPIT_CONFIGS:
        _gapit_trait_files = [assoc_pvalues_trait(_gapit_model, t) for t in get_predictors_list()]
        rule:
            name:   f"assoc_{_gapit_model.lower()}_gea_assemble"
            input:  _gapit_trait_files
            output: assoc_pvalues(_gapit_model)
            params: traits_str = lambda wc, input: " ".join(input)
            log:    f"{LOGDIR}gea/assoc_{_gapit_model.lower()}_assemble.log"
            shell:
                "Rscript /pipeline/scripts/assemble_pvalues.R {params.traits_str} {output} > {log} 2>&1"
