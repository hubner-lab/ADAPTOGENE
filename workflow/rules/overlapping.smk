#=============================================================================
# MODULE 4C: OVERLAPPING REGIONS (GEA + GWAS COMBINED)
#=============================================================================
if PHENO_ASSOC_CONFIGS and ASSOC_CONFIGS:

    rule combine_overlap_snps:
      """Combine GEA and GWAS significant SNPs, compute overlaps, create new combined regions."""
      input:
          gea_selected = O['selected_snps'],
          gwas_selected = O['pheno_selected_snps'],
          gea_regions = O['regions_combined'],
          gwas_regions = O['pheno_regions_combined'],
          ld_decay = ld_decay_input(OVERLAP_REGION_DISTANCE_AUTO)
      output:
          selected = O['overlap_selected_snps'],
          per_trait = O['overlap_regions_per_trait'],
          combined = O['overlap_regions_combined'],
          overlap = O['overlap_summary']
      params:
          overlap_rdist = OVERLAP_REGION_DISTANCE,
          ld_decay_path = O.get('ld_decay_table', 'NULL') if OVERLAP_REGION_DISTANCE_AUTO else 'NULL'
      log: f"{LOGDIR}overlapping/combine_overlap_snps.log"
      shell:
          """
          Rscript /pipeline/scripts/combine_overlap_snps.R \
              {input.gea_selected} {input.gwas_selected} \
              {input.gea_regions} {input.gwas_regions} \
              {params.overlap_rdist} \
              {output.selected} {output.per_trait} {output.combined} \
              {output.overlap} {params.ld_decay_path} > {log} 2>&1
          """
  
    rule find_genes_overlap:
      """Find genes overlapping combined per-trait regions from GEA + GWAS."""
      input:
          regions = O['overlap_regions_per_trait'],
          gff = W['gff_normalized'],
          vcfsnp = W['vcfsnp_full']
      output:
          genes = O['overlap_genes_per_region'],
          collapsed = O['overlap_genes_collapsed']
      params:
          feature = GFF_FEATURE,
          promoter_len = OVERLAP_PROMOTER_LENGTH,
          top_regions = OVERLAP_TOP_REGIONS
      log: f"{LOGDIR}overlapping/find_genes_overlap.log"
      threads: CPU
      shell:
          """
          Rscript /pipeline/scripts/find_genes_around_regions.R \
              {input.gff} {input.regions} {params.feature} \
              {params.promoter_len} {input.vcfsnp} {threads} {params.top_regions} \
              {output.genes} {output.collapsed} > {log} 2>&1
          """
  
    rule find_genes_combined_overlap:
      """Find genes overlapping combined regions from GEA + GWAS (reference)."""
      input:
          regions = O['overlap_regions_combined'],
          gff = W['gff_normalized'],
          vcfsnp = W['vcfsnp_full']
      output: genes = O['overlap_genes_combined']
      params:
          feature = GFF_FEATURE,
          promoter_len = OVERLAP_PROMOTER_LENGTH,
          top_regions = OVERLAP_TOP_REGIONS
      log: f"{LOGDIR}overlapping/find_genes_combined_overlap.log"
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
  
    rule run_enrichment_overlap:
      """Run GO enrichment for overlapping analysis regions."""
      input:
          genes = O['overlap_genes_per_region'],
          gff = W['gff_normalized']
      output: W['overlap_enrichment_done']
      params:
          go_field = GO_FIELD,
          feature = GFF_FEATURE,
          tables_dir = f"{MOD_OVERLAP}tables/enrichment/",
          intermediate_dir = f"{INTER}enrichment_overlapping/"
      log: f"{LOGDIR}overlapping/run_enrichment_overlap.log"
      shell:
          """
          Rscript /pipeline/scripts/run_enrichment.R \
              {input.genes} {input.gff} {params.go_field} {params.feature} \
              {params.tables_dir} {params.intermediate_dir} > {log} 2>&1
          touch {output}
          """
  
    rule plot_enrichment_overlap:
      """Generate enrichment plots for overlapping analysis regions."""
      input:
          enrichment_done = W['overlap_enrichment_done'],
          genes = O['overlap_genes_per_region'],
          regions = O['overlap_regions_per_trait']
      output: W['overlap_enrichment_plots_done']
      params:
          intermediate_dir = f"{INTER}enrichment_overlapping/",
          plots_dir = f"{MOD_OVERLAP}plots/enrichment/",
          top_terms = ENRICHMENT_TOP_TERMS,
          width = ENRICHMENT_PLOT_WIDTH,
          height = ENRICHMENT_PLOT_HEIGHT,
          cnet_label = ENRICHMENT_CNET_LABEL,
          top_plot_regions = ENRICHMENT_TOP_PLOT_REGIONS
      log: f"{LOGDIR}overlapping/plot_enrichment_overlap.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_enrichment.R \
              {params.intermediate_dir} {params.plots_dir} {params.top_terms} \
              {params.width} {params.height} {params.cnet_label} \
              {input.genes} {params.top_plot_regions} {input.regions} > {log} 2>&1
          touch {output}
          """

    rule miami_plot:
      """Static Miami plot combining GEA (top) and GWAS (bottom)."""
      input:
          gea_tables = [assoc_pvalues(m) for m in ASSOC_CONFIGS],
          gwas_tables = [pheno_pvalues(m) for m in PHENO_ASSOC_CONFIGS]
      output:
          png = O['overlap_miami'],
          svg = O['overlap_miami_svg'],
          background = f"{MOD_OVERLAP}plots/miami_combined_K{K_BEST}_background.png",
          coords_json = f"{MOD_OVERLAP}plots/miami_combined_K{K_BEST}_coords.json"
      params:
          gea_str = ','.join([
              f"{method}:{adjust}:{assoc_pvalues(method)}"
              for method, adjust in ASSOC_CONFIGS.items()
          ]),
          gwas_str = ','.join([
              f"{method}:{adjust}:{pheno_pvalues(method)}"
              for method, adjust in PHENO_ASSOC_CONFIGS.items()
          ]),
          gea_preds = PREDICTORS_SELECTED,
          gwas_preds = PHENO_PREDICTORS,
          kbest = K_BEST,
          regions = "NULL",
          sigsnps = "NULL",
          plot_dir = f"{MOD_OVERLAP}plots/",
          scattermore = SCATTERMORE_THRESHOLD
      log: f"{LOGDIR}overlapping/miami_plot.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_miami.R \
              "{params.gea_str}" "{params.gwas_str}" \
              {params.gea_preds} {params.gwas_preds} \
              {params.kbest} {params.regions} {params.sigsnps} \
              {params.plot_dir} {params.scattermore} > {log} 2>&1
          """

    rule miami_plot_regions:
      """Static Miami plot with region rectangles."""
      input:
          gea_tables = [assoc_pvalues(m) for m in ASSOC_CONFIGS],
          gwas_tables = [pheno_pvalues(m) for m in PHENO_ASSOC_CONFIGS],
          regions = O['overlap_regions_combined'],
          sigsnps = O['overlap_selected_snps']
      output:
          png = O['overlap_miami_regions'],
          svg = O['overlap_miami_regions_svg'],
          regions_background = f"{MOD_OVERLAP}plots/miami_combined_K{K_BEST}_regions_background.png"
      params:
          gea_str = ','.join([
              f"{method}:{adjust}:{assoc_pvalues(method)}"
              for method, adjust in ASSOC_CONFIGS.items()
          ]),
          gwas_str = ','.join([
              f"{method}:{adjust}:{pheno_pvalues(method)}"
              for method, adjust in PHENO_ASSOC_CONFIGS.items()
          ]),
          gea_preds = PREDICTORS_SELECTED,
          gwas_preds = PHENO_PREDICTORS,
          kbest = K_BEST,
          plot_dir = f"{MOD_OVERLAP}plots/",
          scattermore = SCATTERMORE_THRESHOLD
      log: f"{LOGDIR}overlapping/miami_plot_regions.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_miami.R \
              "{params.gea_str}" "{params.gwas_str}" \
              {params.gea_preds} {params.gwas_preds} \
              {params.kbest} {input.regions} {input.sigsnps} \
              {params.plot_dir} {params.scattermore} > {log} 2>&1
          """

# Pairwise trait overlap — works with GEA-only, GWAS-only, or both
if ASSOC_CONFIGS or PHENO_ASSOC_CONFIGS:

    rule compute_pairwise_overlaps:
      """Compute all pairwise trait overlaps from method-collapsed sig SNPs."""
      input:
          selected = [f for f in [
              O.get('selected_snps')      if ASSOC_CONFIGS      else None,
              O.get('pheno_selected_snps') if PHENO_ASSOC_CONFIGS else None,
          ] if f is not None]
      output:
          collapsed = O['pairwise_collapsed_snps'],
          pairwise  = O['pairwise_overlap_table']
      params:
          gea_snps  = O.get('selected_snps',       'NULL') if ASSOC_CONFIGS      else 'NULL',
          gwas_snps = O.get('pheno_selected_snps',  'NULL') if PHENO_ASSOC_CONFIGS else 'NULL',
          window    = PAIRWISE_WINDOW_SIZE,
          min_snps  = PAIRWISE_MIN_SNPS
      log: f"{LOGDIR}overlapping/compute_pairwise_overlaps.log"
      shell:
          """
          Rscript /pipeline/scripts/compute_pairwise_overlaps.R \
              {params.gea_snps} {params.gwas_snps} \
              {params.window} {params.min_snps} \
              {output.collapsed} {output.pairwise} > {log} 2>&1
          """
