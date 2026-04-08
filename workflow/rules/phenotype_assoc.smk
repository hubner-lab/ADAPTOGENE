#=============================================================================
# MODULE 4B: PHENOTYPE ASSOCIATION (GWAS)
#=============================================================================

# --- PATH A: MEAN/MEDIAN mode (single sample set, all traits together) ---
if PHENO_ASSOC_CONFIGS and PHENO_MISSING != 'DROP':

    rule prepare_phenotypes:
        """Extract phenotype traits, impute missing values."""
        input: metadata = O['metadata']
        output:
            summary = O['pheno_missing_summary'],
            all_pheno = W['pheno_all_phenotypes'],
            site_means = [f"{WORK_FILT}phenotypes/{trait}_site_means.tsv" for trait in PHENO_TRAITS]
        params:
            strategy = PHENO_MISSING,
            work_dir = f"{WORK_FILT}phenotypes/"
        log: f"{LOGDIR}phenotype_association/prepare_phenotypes.log"
        shell:
            """
            Rscript /pipeline/scripts/prepare_phenotypes.R \
                {input.metadata} {params.strategy} {params.work_dir} {output.summary} > {log} 2>&1
            """

    rule tped_pheno:
        """Convert filtered VCF to TPED/TFAM for phenotype EMMAX."""
        input: vcf = W['vcf_filt']
        output: tped = W['pheno_tped'], tfam = W['pheno_tfam']
        params: prefix = f"{WORK_FILT}phenotypes/emmax/{VCF_BASE}"
        log: f"{LOGDIR}phenotype_association/tped_pheno.log"
        shell:
            """
            plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
                --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
            awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
                {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
            """

    rule kinship_pheno:
        """Compute BN kinship matrix for phenotype EMMAX."""
        input: tped = W['pheno_tped']
        output: W['pheno_kinship']
        params: prefix = f"{WORK_FILT}phenotypes/emmax/{VCF_BASE}"
        log: f"{LOGDIR}phenotype_association/kinship_pheno.log"
        shell:
            "/pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"

    # Non-GAPIT phenotype Path A rules — one rule per method, shell inlined per engine
    for _method in PHENO_OTHER_CONFIGS:
        _engine  = GWAS_METHODS[_method]["engine"]
        _inputs  = _pheno_a_inputs(_engine)
        _params  = _pheno_a_params(_engine, _method)
        _output  = pheno_pvalues(_method)
        _logpath = f"{LOGDIR}phenotype_association/pheno_a_{_method.lower()}.log"

        if _engine == "emmax":
            rule:
                name:   f"pheno_a_{_method.lower()}"
                input:  **_inputs
                output: _output
                params: **_params
                log:    _logpath
                shell:
                    "Rscript /pipeline/scripts/emmax_phenotypes.R "
                    "{input.vcf} {params.tped_prefix} {input.kinship} {input.pca} "
                    "{params.k} {input.phenotypes} {params.tables_dir} NULL {output} > {log} 2>&1"

    if PHENO_GAPIT_CONFIGS:
        rule gapit_pheno:
            """Run GAPIT for all phenotype traits (MEAN/MEDIAN mode)."""
            input:
                gd = W['gapit_gd'],
                gm = W['gapit_gm'],
                pca = W['pca_projections'],
                kinship = W['pheno_kinship'],
                phenotypes = W['pheno_all_phenotypes'],
                metadata = O['metadata']
            output: [pheno_pvalues(model) for model in PHENO_GAPIT_CONFIGS]
            params:
                k = K_BEST,
                models = ','.join(PHENO_GAPIT_CONFIGS.keys()),
                workdir = W['pheno_gapit_work'],
                tables_dir = f"{MOD_PHENO}tables/methods/",
                predictors = PHENO_PREDICTORS,
                native_outdir = f"{MOD_PHENO}GAPIT_native_output/"
            log: f"{LOGDIR}phenotype_association/gapit_pheno.log"
            shell:
                """
                Rscript /pipeline/scripts/gapit.R \
                    {input.gd} {input.gm} {input.phenotypes} {input.pca} \
                    {input.kinship} {params.k} {params.models} \
                    {params.workdir} {params.tables_dir} {params.predictors} \
                    {input.metadata} {params.native_outdir} NULL > {log} 2>&1
                """

# --- PATH B: DROP mode (per-trait sample sets) ---
if PHENO_ASSOC_CONFIGS and PHENO_MISSING == 'DROP':

    rule prepare_phenotypes:
        """Extract phenotype traits, create per-trait sample lists for DROP mode."""
        input: metadata = O['metadata']
        output:
            summary = O['pheno_missing_summary'],
            samples = [f"{WORK_FILT}phenotypes/{trait}_samples.list" for trait in PHENO_TRAITS],
            phenotypes = [f"{WORK_FILT}phenotypes/{trait}_phenotype.tsv" for trait in PHENO_TRAITS],
            site_means = [f"{WORK_FILT}phenotypes/{trait}_site_means.tsv" for trait in PHENO_TRAITS]
        params:
            strategy = PHENO_MISSING,
            work_dir = f"{WORK_FILT}phenotypes/"
        log: f"{LOGDIR}phenotype_association/prepare_phenotypes.log"
        shell:
            """
            Rscript /pipeline/scripts/prepare_phenotypes.R \
                {input.metadata} {params.strategy} {params.work_dir} {output.summary} > {log} 2>&1
            """

    rule subset_vcf_pheno:
        """Subset filtered VCF to per-trait samples."""
        input:
            vcf = W['vcf_filt'],
            samples = f"{WORK_FILT}phenotypes/{{pheno_trait}}_samples.list"
        output: f"{WORK_FILT}phenotypes/{{pheno_trait}}/{VCF_BASE}.vcf"
        wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
        params: prefix = f"{WORK_FILT}phenotypes/{{pheno_trait}}/{VCF_BASE}"
        log: f"{LOGDIR}phenotype_association/subset_vcf_pheno_{{pheno_trait}}.log"
        shell:
            """
            plink --vcf {input.vcf} --keep {input.samples} --const-fid 0 --allow-extra-chr \
                --recode vcf --out {params.prefix} > {log} 2>&1
            sed -i '/^#CHROM/s/\\t0_/\\t/g' {output}
            """

    rule tped_pheno_trait:
        """Convert per-trait VCF to TPED/TFAM."""
        input: vcf = f"{WORK_FILT}phenotypes/{{pheno_trait}}/{VCF_BASE}.vcf"
        output:
            tped = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.tped",
            tfam = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.tfam"
        wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
        params: prefix = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}"
        log: f"{LOGDIR}phenotype_association/tped_pheno_trait_{{pheno_trait}}.log"
        shell:
            """
            plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
                --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
            awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
                {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
            """

    rule kinship_pheno_trait:
        """Compute BN kinship matrix for per-trait sample set."""
        input: tped = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.tped"
        output: f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.aBN.kinf"
        wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
        params: prefix = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}"
        log: f"{LOGDIR}phenotype_association/kinship_pheno_trait_{{pheno_trait}}.log"
        shell:
            "/pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"

    # Non-GAPIT phenotype Path B: per-trait rule + combine rule, one pair per method.
    # Shell strings are inlined per engine to avoid deferred-evaluation closure capture.
    for _method in PHENO_OTHER_CONFIGS:
        _engine           = GWAS_METHODS[_method]["engine"]
        _out_path         = f"{MOD_PHENO}tables/methods/{_method}/{{pheno_trait}}_pvalues_K{K_BEST}.tsv"
        _logpath          = f"{LOGDIR}phenotype_association/pheno_b_{_method.lower()}_{{pheno_trait}}.log"
        _tped_pfx         = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}"
        _per_trait_inputs = expand(
            f"{MOD_PHENO}tables/methods/{_method}/{{trait}}_pvalues_K{K_BEST}.tsv",
            trait=PHENO_TRAITS,
        )

        if _engine == "emmax":
            rule:
                name:   f"pheno_b_{_method.lower()}_trait"
                input:
                    vcf       = f"{WORK_FILT}phenotypes/{{pheno_trait}}/{VCF_BASE}.vcf",
                    kinship   = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.aBN.kinf",
                    pca       = W['pca_projections'],
                    phenotype = f"{WORK_FILT}phenotypes/{{pheno_trait}}_phenotype.tsv",
                output: _out_path
                wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
                params:
                    tped_prefix   = _tped_pfx,
                    k             = K_BEST,
                    tables_dir    = f"{MOD_PHENO}tables/methods/{_method}/",
                    samples_order = W['samples_order'],
                log: _logpath
                shell:
                    "Rscript /pipeline/scripts/emmax_phenotypes.R "
                    "{input.vcf} {params.tped_prefix} {input.kinship} {input.pca} "
                    "{params.k} {input.phenotype} {params.tables_dir} "
                    "{params.samples_order} {output} > {log} 2>&1"

        rule:
            name:  f"combine_pheno_pvalues_{_method.lower()}"
            input: _per_trait_inputs
            output:
                pvals = pheno_pvalues(_method),
                qvals = pheno_qvalues(_method),
            params: files_str = lambda wc, input: ' '.join(input)
            log: f"{LOGDIR}phenotype_association/combine_pheno_pvalues_{_method.lower()}.log"
            shell:
                "Rscript /pipeline/scripts/combine_pheno_pvalues.R "
                '"{params.files_str}" {output.pvals} {output.qvals} > {log} 2>&1'

    if PHENO_GAPIT_CONFIGS:
        rule gapit_pheno_trait:
            """Run GAPIT for a single phenotype trait (DROP mode).
            Uses full-dataset GD and kinship; gapit.R subsets internally."""
            input:
                gd = W['gapit_gd'],
                gm = W['gapit_gm'],
                pca = W['pca_projections'],
                kinship = W['assoc_kinship'],
                phenotype = f"{WORK_FILT}phenotypes/{{pheno_trait}}_phenotype.tsv",
                samples = f"{WORK_FILT}phenotypes/{{pheno_trait}}_samples.list",
                metadata = O['metadata']
            output: [f"{MOD_PHENO}tables/methods/{model}/{{pheno_trait}}_pvalues_K{K_BEST}.tsv" for model in PHENO_GAPIT_CONFIGS]
            wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
            params:
                k = K_BEST,
                models = ','.join(PHENO_GAPIT_CONFIGS.keys()),
                workdir = lambda wc: f"{INTER}gapit/phenotype_association/{wc.pheno_trait}/",
                tables_dir = f"{MOD_PHENO}tables/methods/",
                native_outdir = f"{MOD_PHENO}GAPIT_native_output/"
            log: f"{LOGDIR}phenotype_association/gapit_pheno_trait_{{pheno_trait}}.log"
            shell:
                """
                Rscript /pipeline/scripts/gapit.R \
                    {input.gd} {input.gm} {input.phenotype} {input.pca} \
                    {input.kinship} {params.k} {params.models} \
                    {params.workdir} {params.tables_dir} {wildcards.pheno_trait} \
                    {input.metadata} {params.native_outdir} {input.samples} \
                    {wildcards.pheno_trait} > {log} 2>&1
                """

        rule combine_gapit_pheno_pvalues:
            """Merge per-trait GAPIT p-value files into combined table (DROP mode)."""
            input: lambda wc: expand(f"{MOD_PHENO}tables/methods/{wc.method}/{{trait}}_pvalues_K{K_BEST}.tsv", trait=PHENO_TRAITS)
            output:
                pvals = f"{MOD_PHENO}tables/methods/{{method}}/{{method}}_pvalues_K{K_BEST}.tsv",
                qvals = f"{MOD_PHENO}tables/methods/{{method}}/{{method}}_qvalues_K{K_BEST}.tsv"
            wildcard_constraints:
                method = '|'.join(PHENO_GAPIT_CONFIGS.keys())
            params: files_str = lambda wc, input: ' '.join(input)
            log: f"{LOGDIR}phenotype_association/combine_pheno_pvalues_{{method}}.log"
            shell:
                """
                Rscript /pipeline/scripts/combine_pheno_pvalues.R \
                    "{params.files_str}" {output.pvals} {output.qvals} > {log} 2>&1
                """

# --- Downstream rules (shared by both paths, reuse existing scripts) ---
if PHENO_ASSOC_CONFIGS:

  rule find_sig_snps_pheno:
    """Find significant SNPs for phenotype association."""
    input: assoc = lambda wc: pheno_pvalues(wc.method)
    output: f"{MOD_PHENO}tables/methods/{{method}}/{{method}}_pvalues_K{K_BEST}_sig_snps_{{adjust}}.tsv"
    wildcard_constraints:
        method = PHENO_METHOD_REGEX,
        adjust = r"\w+_[\d.]+"
    params: snp_dist = PHENO_SNP_DISTANCE
    log: f"{LOGDIR}phenotype_association/find_sig_snps_pheno_{{method}}_{{adjust}}.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/find_sig_snps.R \
            {input.assoc} {wildcards.adjust} {params.snp_dist} \
            {wildcards.method} {threads} {output} > {log} 2>&1
        """

  rule combine_selected_snps_pheno:
      """Combine significant SNPs from phenotype association methods."""
      input:
          sigsnps = lambda wc: [pheno_sigsnps(method, adjust) for method, adjust in PHENO_ASSOC_CONFIGS.items()]
      output: O['pheno_selected_snps']
      params:
          sigsnps_str = lambda wc, input: ' '.join(input.sigsnps),
          method = PHENO_COMBINE_METHOD,
          gap = PHENO_COMBINE_GAP,
          predictors = PHENO_PREDICTORS
      log: f"{LOGDIR}phenotype_association/combine_selected_snps_pheno.log"
      shell:
          """
          Rscript /pipeline/scripts/combine_selected_snps.R \
              "{params.sigsnps_str}" {params.method} {params.gap} \
              {params.predictors} {output} > {log} 2>&1
          """
  
  rule create_regions_pheno:
      """Merge nearby significant SNPs into phenotype regions."""
      input:
          selected_snps = O['pheno_selected_snps'],
          ld_decay = ld_decay_input(PHENO_REGION_DISTANCE_AUTO)
      output:
          per_trait = O['pheno_regions_per_trait'],
          combined = O['pheno_regions_combined']
      params:
          region_dist = PHENO_REGION_DISTANCE,
          ld_decay_path = O.get('ld_decay_table', 'NULL') if PHENO_REGION_DISTANCE_AUTO else 'NULL'
      log: f"{LOGDIR}phenotype_association/create_regions_pheno.log"
      shell:
          """
          Rscript /pipeline/scripts/create_regions.R \
              {input.selected_snps} {params.region_dist} \
              {output.per_trait} {output.combined} \
              {params.ld_decay_path} > {log} 2>&1
          """
  
  rule find_genes_pheno:
      """Find genes overlapping phenotype per-trait regions."""
      input:
          regions = O['pheno_regions_per_trait'],
          gff = W['gff_normalized'],
          vcfsnp = W['vcfsnp_full']
      output:
          genes = O['pheno_genes_per_region'],
          collapsed = O['pheno_genes_collapsed']
      params:
          feature = GFF_FEATURE,
          promoter_len = PHENO_PROMOTER_LENGTH
      log: f"{LOGDIR}phenotype_association/find_genes_pheno.log"
      threads: CPU
      shell:
          """
          Rscript /pipeline/scripts/find_genes_around_regions.R \
              {input.gff} {input.regions} {params.feature} \
              {params.promoter_len} {input.vcfsnp} {threads} \
              {output.genes} {output.collapsed} > {log} 2>&1
          """

  rule find_genes_combined_pheno:
      """Find genes overlapping combined phenotype regions."""
      input:
          regions = O['pheno_regions_combined'],
          gff = W['gff_normalized'],
          vcfsnp = W['vcfsnp_full']
      output: genes = O['pheno_genes_combined']
      params:
          feature = GFF_FEATURE,
          promoter_len = PHENO_PROMOTER_LENGTH
      log: f"{LOGDIR}phenotype_association/find_genes_combined_pheno.log"
      threads: CPU
      shell:
          """
          TEMP_COLLAPSED=$(mktemp)
          Rscript /pipeline/scripts/find_genes_around_regions.R \
              {input.gff} {input.regions} {params.feature} \
              {params.promoter_len} {input.vcfsnp} {threads} \
              {output.genes} $TEMP_COLLAPSED > {log} 2>&1
          rm -f $TEMP_COLLAPSED
          """
  
  rule run_enrichment_pheno:
      """Run GO enrichment for phenotype association regions."""
      input:
          genes = O['pheno_genes_per_region'],
          gff = W['gff_normalized']
      output: W['pheno_enrichment_done']
      params:
          go_field = GO_FIELD,
          feature = GFF_FEATURE,
          tables_dir = f"{MOD_PHENO}tables/enrichment/",
          intermediate_dir = f"{INTER}enrichment_phenotypes/"
      log: f"{LOGDIR}phenotype_association/run_enrichment_pheno.log"
      shell:
          """
          Rscript /pipeline/scripts/run_enrichment.R \
              {input.genes} {input.gff} {params.go_field} {params.feature} \
              {params.tables_dir} {params.intermediate_dir} > {log} 2>&1
          touch {output}
          """
  
  # Manhattan plots for phenotype traits
  rule manhattan_pheno:
      """Generate Manhattan plot and QQ plot for a phenotype trait."""
      input: assoc = lambda wc: pheno_pvalues(wc.method)
      output:
          png = f"{MOD_PHENO}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}.png",
          svg = f"{MOD_PHENO}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}.svg",
          qq_png = f"{MOD_PHENO}plots/manhattan/{{method}}/qq_{{trait}}_K{K_BEST}_{{adjust}}.png",
          qq_svg = f"{MOD_PHENO}plots/manhattan/{{method}}/qq_{{trait}}_K{K_BEST}_{{adjust}}.svg",
          background = f"{MOD_PHENO}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}_background.png",
          coords_json = f"{MOD_PHENO}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}_coords.json"
      wildcard_constraints:
          method = PHENO_METHOD_REGEX,
          trait = r"[a-zA-Z]\w*",
          adjust = r"\w+_[\d.]+"
      params:
          k = K_BEST,
          plot_dir = lambda wc: f"{MOD_PHENO}plots/manhattan/{wc.method}/",
          predictors = PHENO_PREDICTORS
      log: f"{LOGDIR}phenotype_association/manhattan_pheno_{{method}}_{{trait}}_{{adjust}}.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_manhattan.R \
              {input.assoc} {wildcards.adjust} {params.k} {wildcards.method} \
              {wildcards.trait} {params.plot_dir} {params.predictors} > {log} 2>&1
          """

  rule manhattan_combined_pheno:
      """Generate combined Manhattan and QQ plots for all phenotype traits."""
      input:
          assoc_tables = [pheno_pvalues(method) for method in PHENO_ASSOC_CONFIGS]
      output:
          simple_png = O['pheno_manhattan_combined'],
          simple_svg = f"{MOD_PHENO}plots/manhattan/combined/manhattan_combined_K{K_BEST}.svg",
          qq_png = O['pheno_qq_combined'],
          qq_svg = f"{MOD_PHENO}plots/manhattan/combined/qq_combined_K{K_BEST}.svg",
          background = f"{MOD_PHENO}plots/manhattan/combined/manhattan_combined_K{K_BEST}_background.png",
          coords_json = f"{MOD_PHENO}plots/manhattan/combined/manhattan_combined_K{K_BEST}_coords.json"
      params:
          assoc_str = ','.join([
              f"{method}:{adjust}:{pheno_pvalues(method)}"
              for method, adjust in PHENO_ASSOC_CONFIGS.items()
          ]),
          predictors = PHENO_PREDICTORS,
          k = K_BEST,
          plot_dir = f"{MOD_PHENO}plots/manhattan/combined/"
      log: f"{LOGDIR}phenotype_association/manhattan_combined_pheno.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_manhattan_combined.R \
              "{params.assoc_str}" {params.predictors} {params.k} \
              {params.plot_dir} > {log} 2>&1
          touch {output.simple_png} {output.simple_svg} {output.qq_png} {output.qq_svg} {output.background} {output.coords_json}
          """
  
  rule piemap_pheno:
      """Generate pie map for phenotype trait with trait-controlled pie sizes."""
      input:
          raster = W['climate_raster'],
          metadata = O['metadata'],
          clusters = clusters_table(K_BEST),
          trait = f"{WORK_FILT}phenotypes/{{pheno_trait}}_site_means.tsv"
      output:
          png = f"{MOD_PHENO}plots/piemap/phenomap_{{pheno_trait}}.png",
          svg = f"{MOD_PHENO}plots/piemap/phenomap_{{pheno_trait}}.svg"
      wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
      params:
          raster_layer = get_predictors_list()[0] if get_predictors_list() else "1",
          pie_alpha = PIEMAP_ALPHA,
          pop_label = PIEMAP_SHOW_LABELS,
          pop_label_size = PIEMAP_LABEL_SIZE,
          plot_dir = f"{MOD_PHENO}plots/piemap/",
          inter_dir = INTER,
          regionmap_extent = REGIONMAP_EXTENT,
          pie_scale = PIEMAP_PIE_SCALE,
          use_points = PIEMAP_USE_POINTS
      log: f"{LOGDIR}phenotype_association/piemap_pheno_{{pheno_trait}}.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_piemap.R \
              {input.raster} {params.raster_layer} {params.raster_layer} \
              {input.metadata} {input.clusters} \
              {input.trait} "{wildcards.pheno_trait}" \
              {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
              {params.plot_dir} {params.inter_dir} \
              phenomap_{wildcards.pheno_trait} {params.regionmap_extent} {params.pie_scale} Clusters {params.use_points} > {log} 2>&1
          """
