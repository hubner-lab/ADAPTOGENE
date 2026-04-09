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

# Non-GAPIT GEA rules — one rule per method, driven by registry.
# Shell strings are inlined per engine (not from a variable) to avoid Snakemake's
# deferred shell evaluation capturing the loop variable by reference.
for _method in GEA_OTHER_CONFIGS:
    _engine  = GEA_METHODS[_method]["engine"]
    _inputs  = _gea_inputs(_engine)
    _params  = _gea_params(_engine, _method)
    _output  = assoc_pvalues(_method)
    _logpath = f"{LOGDIR}gea/assoc_{_method.lower()}.log"

    if _engine == "emmax":
        rule:
            name:   f"assoc_{_method.lower()}"
            input:  **_inputs
            output: _output
            params: **_params
            log:    _logpath
            shell:
                "Rscript /pipeline/scripts/emmax.R "
                "{input.vcf} {params.k} {input.traits} {input.covariates} "
                "{params.predictors} {params.inter_dir} {input.metadata} "
                "{params.tables_dir} {params.tped_prefix} {input.kinship} > {log} 2>&1"
    elif _engine == "lfmm":
        rule:
            name:   f"assoc_{_method.lower()}"
            input:  **_inputs
            output: _output
            params: **_params
            log:    _logpath
            shell:
                "Rscript /pipeline/scripts/lfmm.R "
                "{input.lfmm_ld} {input.lfmm_full} {input.climate} "
                "{params.k} {params.predictors} {input.vcfsnp} "
                "{params.tables_dir} > {log} 2>&1"

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

# GAPIT GEA analysis — one rule emits all configured GAPIT model outputs
if GEA_GAPIT_CONFIGS:

    rule gapit_analysis:
        """Run GAPIT association analysis for all GAPIT models and traits."""
        input:
            gd = W['gapit_gd'],
            gm = W['gapit_gm'],
            traits = O['climate_site_scaled'],
            pca = W['pca_projections'],
            kinship = W['assoc_kinship'],
            metadata = O['metadata']
        output: [assoc_pvalues(model) for model in GEA_GAPIT_CONFIGS]
        params:
            k = K_BEST,
            models = ','.join(GEA_GAPIT_CONFIGS.keys()),
            workdir = W['gapit_work'],
            tables_dir = f"{MOD_GEA}tables/methods/",
            predictors = PREDICTORS_SELECTED,
            native_outdir = f"{MOD_GEA}GAPIT_native_output/"
        log: f"{LOGDIR}gea/gapit_analysis.log"
        shell:
            """
            Rscript /pipeline/scripts/gapit.R \
                {input.gd} {input.gm} {input.traits} {input.pca} \
                {input.kinship} {params.k} {params.models} \
                {params.workdir} {params.tables_dir} {params.predictors} \
                {input.metadata} {params.native_outdir} NULL > {log} 2>&1
            """
