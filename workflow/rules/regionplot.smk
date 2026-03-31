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

# rule regionplot: DEPRECATED — regionplots are now generated on-demand in the Shiny app
# (Association tab > Region Explorer > Generate Region Plot)
