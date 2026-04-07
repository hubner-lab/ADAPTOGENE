#=============================================================================
# MODULE 4C: OVERLAPPING REGIONS (GEA + GWAS COMBINED)
#=============================================================================
if PHENO_ASSOC_CONFIGS and ASSOC_CONFIGS:

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
          plot_dir = f"{MOD_OVERLAP}plots/"
      log: f"{LOGDIR}overlapping/miami_plot.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_miami.R \
              "{params.gea_str}" "{params.gwas_str}" \
              {params.gea_preds} {params.gwas_preds} \
              {params.kbest} {params.plot_dir} > {log} 2>&1
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
