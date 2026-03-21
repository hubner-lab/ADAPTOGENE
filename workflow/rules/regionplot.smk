#=============================================================================
# MODULE 5: REGIONPLOT
#=============================================================================

rule gff2topr:
    """Convert GFF to topr-compatible gene annotation format."""
    input: gff = W['gff_normalized']
    output: O['gff_topr']
    params:
        feature = GFF_FEATURE,
        genename = GFF_GENE_NAME,
        biotype = GFF_BIOTYPE
    log: f"{LOGDIR}association/gff2topr.log"
    shell:
        """
        python3 /pipeline/scripts/gff2topr.py \
            {input.gff} {params.feature} {params.genename} {params.biotype} \
            {output} > {log} 2>&1
        """

rule regionplot:
    """Generate regional Manhattan plots for top regions with all methods overlaid."""
    input:
        regions = O['regions_per_trait'],
        gff_topr = O['gff_topr'],
        assoc_tables = [assoc_pvalues(method) for method in ASSOC_CONFIGS]
    output: touch(O['regionplot_done'])
    params:
        assoc_str = ','.join([
            f"{method}:{adjust}:{assoc_pvalues(method)}"
            for method, adjust in ASSOC_CONFIGS.items()
        ]),
        top_regions = TOP_REGIONS,
        genes = GENES_TO_HIGHLIGHT,
        plot_dir = f"{MOD_REGPLOT}",
        custom_region = REGIONPLOT_REGION,
        custom_traits = REGIONPLOT_TRAITS,
        custom_methods = REGIONPLOT_METHOD
    log: f"{LOGDIR}association/regionplot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_regionplot.R \
            {input.regions} {input.gff_topr} {params.assoc_str} \
            {params.top_regions} {params.genes} {params.plot_dir} \
            {params.custom_region} {params.custom_traits} \
            {params.custom_methods} > {log} 2>&1
        """
