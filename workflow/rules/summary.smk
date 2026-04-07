#=============================================================================
# MODULE 9: PIPELINE SUMMARY
#=============================================================================

# Summary rule uses ruleorder to resolve ambiguity - only the matching mode runs
if MODE == 'processing':
    rule write_summary:
        """Write processing mode summary to Pipeline_summary.tsv."""
        input:
            vcf_filt           = W['vcf_filt'],
            vcf_ld             = W['vcf_ld'],
            samples_list       = W['samples_list'],
            samples_filtered   = W['samples_filtered'],
            samples_removed    = W['samples_removed'],
            qc_raw_summary     = O['qc_raw_summary'],
            filtering_summary  = O['qc_filtering_summary'],
            plot_done          = O['qc_plot_attrition'],
            depth_summary      = O['qc_depth_site'] if HAS_FORMAT_DP else []
        output: W['summary_done']
        params:
            summary_tsv    = O['summary'],
            het_outlier_sd = HET_OUTLIER_SD if HET_OUTLIER_SD is not None else 'NULL',
            has_dp         = 'TRUE' if HAS_FORMAT_DP else 'FALSE',
            depth_summary  = f"{MOD_PROC}tables/depth_summary.tsv" if HAS_FORMAT_DP else 'NULL'
        log: f"{LOGDIR}processing/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                processing {params.summary_tsv} \
                {input.vcf_filt} {input.vcf_ld} \
                {input.samples_list} {input.samples_filtered} {input.samples_removed} \
                {input.qc_raw_summary} {input.filtering_summary} \
                {params.het_outlier_sd} {params.has_dp} {params.depth_summary} > {log} 2>&1
            touch {output}
            """

elif MODE == 'structure':
    rule write_summary:
        """Write structure mode summary to Pipeline_summary.tsv."""
        input:
            cross_entropy = O['cross_entropy']
        output: W['summary_done']
        params: ks = K_START, ke = K_END, summary_tsv = O['summary']
        log: f"{LOGDIR}structure/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                structure {params.summary_tsv} \
                {params.ks} {params.ke} > {log} 2>&1
            touch {output}
            """

elif MODE == 'structure_K':
    rule write_summary:
        """Write structure_K mode summary to Pipeline_summary.tsv."""
        input:
            climate_site = O['climate_site'] if CLIMATE_ENABLED else [],
            ld_decay = O['ld_decay_table']
        output: W['summary_done']
        params:
            k = K_BEST,
            predictors = ALL_BIO_STR if CLIMATE_ENABLED else 'NULL',
            climate_site_path = O['climate_site'] if CLIMATE_ENABLED else 'NULL',
            summary_tsv = O['summary'],
            ld_decay_path = O['ld_decay_table'],
            ld_decay_group_by = LD_DECAY_GROUP_BY,
            ld_decay_scope = LD_DECAY_SCOPE
        log: f"{LOGDIR}structure_k/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                structure_K {params.summary_tsv} \
                {params.k} {params.climate_site_path} {params.predictors} \
                {params.ld_decay_path} {params.ld_decay_group_by} {params.ld_decay_scope} > {log} 2>&1
            touch {output}
            """

elif MODE == 'association':
    rule write_summary:
        """Write association mode summary to Pipeline_summary.tsv."""
        input:
            selected_snps = O['selected_snps'],
            regions_per_trait = O['regions_per_trait'],
            regions_combined = O['regions_combined'],
            genes = O['genes_per_region'],
            enrichment_done = W['enrichment_done'] if (GO_FIELD and GO_FIELD != 'NULL') else []
        output: W['summary_done']
        params:
            enrichment_path = W['enrichment_done'] if (GO_FIELD and GO_FIELD != 'NULL') else 'NULL',
            summary_tsv = O['summary']
        log: f"{LOGDIR}association/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                association {params.summary_tsv} \
                {input.selected_snps} {input.regions_per_trait} {input.regions_combined} \
                {input.genes} {params.enrichment_path} > {log} 2>&1
            touch {output}
            """

elif MODE == 'maladaptation':
    rule write_summary:
        """Write maladaptation mode summary to Pipeline_summary.tsv."""
        input:
            gf_adaptive = W['gf_adaptive'],
            offset_site = O['gf_offset_site_values']
        output: W['summary_done']
        params: summary_tsv = O['summary']
        log: f"{LOGDIR}maladaptation/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                maladaptation {params.summary_tsv} \
                {input.gf_adaptive} {input.offset_site} > {log} 2>&1
            touch {output}
            """

elif MODE in ('association_phenotypes', 'phenotype_association'):
    rule write_summary:
        """Write association_phenotypes mode summary to Pipeline_summary.tsv."""
        input:
            missing_summary = O['pheno_missing_summary'],
            selected_snps = O['pheno_selected_snps'],
            regions_per_trait = O['pheno_regions_per_trait'],
            regions_combined = O['pheno_regions_combined'],
            genes = O['pheno_genes_per_region'],
            enrichment_done = W['pheno_enrichment_done'] if (GO_FIELD and GO_FIELD != 'NULL') else []
        output: W['summary_done']
        params:
            enrichment_path = W['pheno_enrichment_done'] if (GO_FIELD and GO_FIELD != 'NULL') else 'NULL',
            summary_tsv = O['summary']
        log: f"{LOGDIR}phenotype_association/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                association_phenotypes {params.summary_tsv} \
                {input.missing_summary} {input.selected_snps} {input.regions_per_trait} \
                {input.regions_combined} {input.genes} {params.enrichment_path} > {log} 2>&1
            touch {output}
            """

elif MODE == 'overlapping':
    rule write_summary:
        """Overlapping mode: region analysis is fully interactive (Shiny). Just mark done."""
        input:
            [O['overlap_miami']] if (ASSOC_CONFIGS and PHENO_ASSOC_CONFIGS) else [],
            [O['pairwise_overlap_table']] if (ASSOC_CONFIGS or PHENO_ASSOC_CONFIGS) else [],
        output: W['summary_done']
        shell: "touch {output}"

elif MODE == 'haplotype_scan':
    _hap_scan_inputs = {}
    for i, tag in enumerate(HAP_TAGS):
        _hap_scan_inputs[f'scan_done_{i}'] = W[f'hap_scan_done_{tag}']
        _hap_scan_inputs[f'regions_{i}'] = O[f'hap_selected_regions_{tag}']

    rule write_summary:
        """Write haplotype_scan mode summary to Pipeline_summary.tsv."""
        input: **_hap_scan_inputs
        output: W['summary_done']
        params:
            summary_tsv = O['summary'],
            tags = ','.join(HAP_TAGS)
        log: f"{LOGDIR}haplotype_scan/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                haplotype_scan {params.summary_tsv} \
                {params.tags} > {log} 2>&1
            touch {output}
            """

elif MODE == 'haplotype':
    _hap_viz_inputs = {}
    for i, tag in enumerate(HAP_TAGS):
        _hap_viz_inputs[f'viz_done_{i}'] = W[f'hap_viz_done_{tag}']

    rule write_summary:
        """Write haplotype mode summary to Pipeline_summary.tsv."""
        input: **_hap_viz_inputs
        output: W['summary_done']
        params:
            summary_tsv = O['summary'],
            tags = ','.join(HAP_TAGS),
            epsilon = HAP_EPSILON_SELECTED
        log: f"{LOGDIR}haplotype/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                haplotype {params.summary_tsv} \
                {params.tags} {params.epsilon} > {log} 2>&1
            touch {output}
            """
