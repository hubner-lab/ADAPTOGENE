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
        """Compute kinship matrix for phenotype EMMAX."""
        input: tped = W['pheno_tped']
        output: W['pheno_kinship']
        params: prefix = f"{WORK_FILT}phenotypes/emmax/{VCF_BASE}"
        log: f"{LOGDIR}phenotype_association/kinship_pheno.log"
        shell:
            "/pipeline/scripts/emmax-kin-intel64 -v -s -d 10 -x {params.prefix} > {log} 2>&1"

    rule emmax_pheno:
        """Run EMMAX for all phenotype traits (MEAN/MEDIAN mode)."""
        input:
            vcf = W['vcf_filt'],
            tped = W['pheno_tped'],
            kinship = W['pheno_kinship'],
            pca = W['pca_projections'],
            phenotypes = W['pheno_all_phenotypes']
        output: pheno_pvalues("EMMAX")
        params:
            tped_prefix = f"{WORK_FILT}phenotypes/emmax/{VCF_BASE}",
            k = K_BEST,
            tables_dir = f"{MOD_PHENO}tables/"
        log: f"{LOGDIR}phenotype_association/emmax_pheno.log"
        shell:
            """
            Rscript /pipeline/scripts/emmax_phenotypes.R \
                {input.vcf} {params.tped_prefix} {input.kinship} {input.pca} \
                {params.k} {input.phenotypes} {params.tables_dir} NULL {output} > {log} 2>&1
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
        """Compute kinship matrix for per-trait sample set."""
        input: tped = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.tped"
        output: f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.aIBS.kinf"
        wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
        params: prefix = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}"
        log: f"{LOGDIR}phenotype_association/kinship_pheno_trait_{{pheno_trait}}.log"
        shell:
            "/pipeline/scripts/emmax-kin-intel64 -v -s -d 10 -x {params.prefix} > {log} 2>&1"

    rule emmax_pheno_trait:
        """Run EMMAX for a single phenotype trait (DROP mode)."""
        input:
            vcf = f"{WORK_FILT}phenotypes/{{pheno_trait}}/{VCF_BASE}.vcf",
            kinship = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.aIBS.kinf",
            pca = W['pca_projections'],
            phenotype = f"{WORK_FILT}phenotypes/{{pheno_trait}}_phenotype.tsv"
        output: f"{MOD_PHENO}tables/EMMAX/{{pheno_trait}}_pvalues_K{K_BEST}.tsv"
        wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
        params:
            tped_prefix = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}",
            k = K_BEST,
            tables_dir = f"{MOD_PHENO}tables/EMMAX/",
            samples_order = W['samples_order']
        log: f"{LOGDIR}phenotype_association/emmax_pheno_trait_{{pheno_trait}}.log"
        shell:
            """
            Rscript /pipeline/scripts/emmax_phenotypes.R \
                {input.vcf} {params.tped_prefix} {input.kinship} {input.pca} \
                {params.k} {input.phenotype} {params.tables_dir} {params.samples_order} \
                {output} > {log} 2>&1
            """

    rule combine_pheno_pvalues:
        """Merge per-trait p-value files into combined table (DROP mode)."""
        input: expand(f"{MOD_PHENO}tables/EMMAX/{{trait}}_pvalues_K{K_BEST}.tsv", trait=PHENO_TRAITS)
        output: pheno_pvalues("EMMAX")
        params: files_str = lambda wc, input: ' '.join(input)
        log: f"{LOGDIR}phenotype_association/combine_pheno_pvalues.log"
        shell:
            """
            Rscript /pipeline/scripts/combine_pheno_pvalues.R \
                "{params.files_str}" {output} > {log} 2>&1
            """

# --- Downstream rules (shared by both paths, reuse existing scripts) ---
if PHENO_ASSOC_CONFIGS:

  rule find_sig_snps_pheno:
    """Find significant SNPs for phenotype association."""
    input: assoc = lambda wc: pheno_pvalues(wc.method)
    output: f"{MOD_PHENO}tables/{{method}}/{{method}}_phenotypes_pvalues_K{K_BEST}_sig_snps_{{adjust}}.tsv"
    wildcard_constraints:
        method = r"EMMAX",
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
      input: selected_snps = O['pheno_selected_snps']
      output:
          per_trait = O['pheno_regions_per_trait'],
          combined = O['pheno_regions_combined']
      params: region_dist = PHENO_REGION_DISTANCE
      log: f"{LOGDIR}phenotype_association/create_regions_pheno.log"
      shell:
          """
          Rscript /pipeline/scripts/create_regions.R \
              {input.selected_snps} {params.region_dist} \
              {output.per_trait} {output.combined} > {log} 2>&1
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
          promoter_len = PHENO_PROMOTER_LENGTH,
          top_regions = PHENO_TOP_REGIONS
      log: f"{LOGDIR}phenotype_association/find_genes_pheno.log"
      threads: CPU
      shell:
          """
          Rscript /pipeline/scripts/find_genes_around_regions.R \
              {input.gff} {input.regions} {params.feature} \
              {params.promoter_len} {input.vcfsnp} {threads} {params.top_regions} \
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
          promoter_len = PHENO_PROMOTER_LENGTH,
          top_regions = PHENO_TOP_REGIONS
      log: f"{LOGDIR}phenotype_association/find_genes_combined_pheno.log"
      threads: CPU
      shell:
          """
          TEMP_COLLAPSED=$(mktemp)
          Rscript /pipeline/scripts/find_genes_around_regions.R \
              {input.gff} {input.regions} {params.feature} \
              {params.promoter_len} {input.vcfsnp} {threads} {params.top_regions} \
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
  
  rule plot_enrichment_pheno:
      """Generate enrichment plots for phenotype regions."""
      input:
          enrichment_done = W['pheno_enrichment_done'],
          genes = O['pheno_genes_per_region'],
          regions = O['pheno_regions_per_trait']
      output: W['pheno_enrichment_plots_done']
      params:
          intermediate_dir = f"{INTER}enrichment_phenotypes/",
          plots_dir = f"{MOD_PHENO}plots/enrichment/",
          top_terms = ENRICHMENT_TOP_TERMS,
          width = ENRICHMENT_PLOT_WIDTH,
          height = ENRICHMENT_PLOT_HEIGHT,
          cnet_label = ENRICHMENT_CNET_LABEL,
          top_plot_regions = ENRICHMENT_TOP_PLOT_REGIONS
      log: f"{LOGDIR}phenotype_association/plot_enrichment_pheno.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_enrichment.R \
              {params.intermediate_dir} {params.plots_dir} {params.top_terms} \
              {params.width} {params.height} {params.cnet_label} \
              {input.genes} {params.top_plot_regions} {input.regions} > {log} 2>&1
          touch {output}
          """
  
  # Manhattan plots for phenotype traits
  rule manhattan_pheno:
      """Generate Manhattan plot for a phenotype trait."""
      input: assoc = lambda wc: pheno_pvalues(wc.method)
      output:
          png = f"{MOD_PHENO}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}.png",
          svg = f"{MOD_PHENO}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}.svg"
      wildcard_constraints:
          method = r"EMMAX",
          trait = r"[a-zA-Z]\w*",
          adjust = r"\w+_[\d.]+"
      params:
          k = K_BEST,
          plot_dir = lambda wc: f"{MOD_PHENO}plots/manhattan/{wc.method}/",
          regions = "NULL",
          selected_snps = "NULL",
          scattermore_threshold = PHENO_SCATTERMORE_THRESHOLD
      log: f"{LOGDIR}phenotype_association/manhattan_pheno_{{method}}_{{trait}}_{{adjust}}.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_manhattan.R \
              {input.assoc} {wildcards.adjust} {params.k} {wildcards.method} \
              {wildcards.trait} {params.plot_dir} {params.regions} {params.selected_snps} \
              {params.scattermore_threshold} > {log} 2>&1
          """
  
  rule manhattan_pheno_regions:
      """Generate Manhattan plot with phenotype regions highlighted."""
      input:
          assoc = lambda wc: pheno_pvalues(wc.method),
          regions = O['pheno_regions_per_trait'],
          sigsnps = lambda wc: [pheno_sigsnps(method, adjust) for method, adjust in PHENO_ASSOC_CONFIGS.items()]
      output:
          png = f"{MOD_PHENO}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}_regions.png",
          svg = f"{MOD_PHENO}plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}_regions.svg"
      wildcard_constraints:
          method = r"EMMAX",
          trait = r"[a-zA-Z]\w*",
          adjust = r"\w+_[\d.]+"
      params:
          k = K_BEST,
          plot_dir = lambda wc: f"{MOD_PHENO}plots/manhattan/{wc.method}/",
          sigsnps_str = lambda wc, input: ','.join(input.sigsnps),
          scattermore_threshold = PHENO_SCATTERMORE_THRESHOLD
      log: f"{LOGDIR}phenotype_association/manhattan_pheno_regions_{{method}}_{{trait}}_{{adjust}}.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_manhattan.R \
              {input.assoc} {wildcards.adjust} {params.k} {wildcards.method} \
              {wildcards.trait} {params.plot_dir} {input.regions} "{params.sigsnps_str}" \
              {params.scattermore_threshold} > {log} 2>&1
          """
  
  rule manhattan_combined_pheno:
      """Generate combined Manhattan plots for all phenotype traits."""
      input:
          assoc_tables = [pheno_pvalues(method) for method in PHENO_ASSOC_CONFIGS],
          regions = O['pheno_regions_combined']
      output:
          simple_png = O['pheno_manhattan_combined'],
          simple_svg = f"{MOD_PHENO}plots/manhattan_combined_K{K_BEST}.svg",
          regions_png = O['pheno_manhattan_combined_regions'],
          regions_svg = f"{MOD_PHENO}plots/manhattan_combined_K{K_BEST}_regions.svg"
      params:
          assoc_str = ','.join([
              f"{method}:{adjust}:{pheno_pvalues(method)}"
              for method, adjust in PHENO_ASSOC_CONFIGS.items()
          ]),
          predictors = PHENO_PREDICTORS,
          k = K_BEST,
          plot_dir = f"{MOD_PHENO}plots/"
      log: f"{LOGDIR}phenotype_association/manhattan_combined_pheno.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_manhattan_combined.R \
              "{params.assoc_str}" {params.predictors} {params.k} \
              {input.regions} {params.plot_dir} > {log} 2>&1
          touch {output.simple_png} {output.simple_svg} {output.regions_png} {output.regions_svg}
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
