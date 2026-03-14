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
          gwas_regions = O['pheno_regions_combined']
      output:
          selected = O['overlap_selected_snps'],
          per_trait = O['overlap_regions_per_trait'],
          combined = O['overlap_regions_combined'],
          overlap = O['overlap_summary']
      params:
          overlap_rdist = OVERLAP_REGION_DISTANCE
      log: f"{LOGDIR}overlapping/combine_overlap_snps.log"
      shell:
          """
          Rscript /pipeline/scripts/combine_overlap_snps.R \
              {input.gea_selected} {input.gwas_selected} \
              {input.gea_regions} {input.gwas_regions} \
              {params.overlap_rdist} \
              {output.selected} {output.per_trait} {output.combined} \
              {output.overlap} > {log} 2>&1
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
