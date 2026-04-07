#=============================================================================
# MODULE 8: HAPLOTYPE ANALYSIS (crosshap)
#=============================================================================

def get_hap_regions_input(source):
    """Get the regions file path for a given haplotype source."""
    if source == 'association':
        return O.get('regions_combined', '')
    elif source == 'association_phenotypes':
        return O.get('pheno_regions_combined', '')
    elif source == 'custom':
        # Custom path is relative to pipeline root, convert to absolute
        return os.path.join('/pipeline', HAP_SCAN_REGIONS_FILE) if HAP_SCAN_REGIONS_FILE != 'NULL' else ''
    return ''

def get_hap_meta_type_from_tag(tag):
    """Extract metadata type ('site' or 'cluster_K{N}') from a combined tag."""
    # tag = "site_association" or "cluster_K5_association_phenotypes"
    for mt in HAP_META_TYPES:
        if tag.startswith(mt + '_'):
            return mt
    return 'site'

def get_hap_source_from_tag(tag):
    """Extract source from a combined tag."""
    for mt in HAP_META_TYPES:
        prefix = mt + '_'
        if tag.startswith(prefix):
            return tag[len(prefix):]
    return tag

# Generate haplotype rules for each (meta_tag, source) combination
for _hap_tag in HAP_TAGS:
    _hap_meta = get_hap_meta_type_from_tag(_hap_tag)
    _hap_source = get_hap_source_from_tag(_hap_tag)
    _hap_regions_input = get_hap_regions_input(_hap_source)

    # Rule: select top regions for haplotype analysis
    rule:
        name: f"haplotype_select_regions_{_hap_tag}"
        input:
            regions = _hap_regions_input,
            metadata = O['metadata']
        output: O[f'hap_selected_regions_{_hap_tag}']
        params:
            top_regions = HAP_SCAN_TOP_REGIONS,
            source = _hap_source,
            tag = _hap_tag
        log: f"{LOGDIR}haplotype_scan/haplotype_select_regions_{_hap_tag}.log"
        shell:
            """
            Rscript /pipeline/scripts/select_haplotype_regions.R \
                {input.regions} {output} {params.top_regions} \
                {params.source} > {log} 2>&1
            """

    # Rule: scan with epsilon range (crosshap haplotyping + clustree)
    _scan_inputs = {
        'selected_regions': O[f'hap_selected_regions_{_hap_tag}'],
        'vcf_filt': W['vcf_filt'],
        'metadata': O['metadata'],
    }
    # Add imputed VCF if available (from association mode), else use filtered
    if 'vcf_imp_full' in W:
        _scan_inputs['vcf_ld'] = W['vcf_imp_full']
    else:
        _scan_inputs['vcf_ld'] = W['vcf_filt']
    # Add clusters file for cluster metadata type
    if _hap_meta.startswith('cluster'):
        _scan_inputs['clusters'] = clusters_table(K_BEST)

    rule:
        name: f"haplotype_scan_{_hap_tag}"
        input: **_scan_inputs
        output: touch(W[f'hap_scan_done_{_hap_tag}'])
        params:
            meta_type = _hap_meta,
            clusters_file = clusters_table(K_BEST) if _hap_meta.startswith('cluster') else 'NULL',
            mgmin = HAP_SCAN_MGMIN,
            minhap = HAP_SCAN_MINHAP,
            min_snps = HAP_SCAN_MIN_SNPS,
            epsilon_range = ','.join(str(e) for e in HAP_SCAN_EPSILON_RANGE),
            inter_dir = f"{INTER}haplotype/{_hap_tag}/",
            plots_dir = f"{OUTDIR}haplotype_scan/{_hap_tag}/plots/clustree/",
            tag = _hap_tag
        log: f"{LOGDIR}haplotype_scan/haplotype_scan_{_hap_tag}.log"
        shell:
            """
            Rscript /pipeline/scripts/run_haplotype_scan.R \
                {input.selected_regions} {input.vcf_filt} {input.vcf_ld} \
                {input.metadata} {params.meta_type} {params.clusters_file} \
                {params.mgmin} {params.minhap} {params.epsilon_range} \
                {params.min_snps} {params.inter_dir} {params.plots_dir} > {log} 2>&1
            """

    # Rule: visualization at selected epsilon
    _viz_inputs = {
        'scan_done': W[f'hap_scan_done_{_hap_tag}'],
        'selected_regions': O[f'hap_selected_regions_{_hap_tag}'],
        'metadata': O['metadata'],
        'raster': W.get('climate_raster', ''),
        'clusters': clusters_table(K_BEST) if K_BEST else '',
    }

    rule:
        name: f"haplotype_viz_{_hap_tag}"
        input: **_viz_inputs
        output: touch(W[f'hap_viz_done_{_hap_tag}'])
        params:
            epsilon = HAP_EPSILON_SELECTED,
            meta_type = _hap_meta,
            tag = _hap_tag,
            inter_dir = f"{INTER}haplotype/{_hap_tag}/",
            plots_dir = f"{OUTDIR}haplotype/{_hap_tag}/plots/",
            tables_dir = f"{OUTDIR}haplotype/{_hap_tag}/tables/",
            raster_layer = get_predictors_list()[0] if get_predictors_list() else "1",
            pie_alpha = PIEMAP_ALPHA,
            pop_label = PIEMAP_SHOW_LABELS,
            pop_label_size = PIEMAP_LABEL_SIZE,
            regionmap_extent = REGIONMAP_EXTENT,
            pie_scale = PIEMAP_PIE_SCALE,
            use_points = PIEMAP_USE_POINTS
        log: f"{LOGDIR}haplotype/haplotype_viz_{_hap_tag}.log"
        shell:
            """
            Rscript /pipeline/scripts/run_haplotype_viz.R \
                {input.selected_regions} {input.metadata} \
                {params.epsilon} {params.meta_type} {params.tag} \
                {params.inter_dir} {params.plots_dir} {params.tables_dir} \
                {input.raster} {params.raster_layer} \
                {params.pie_alpha},{params.pop_label},{params.pop_label_size},{params.pie_scale},{params.use_points} \
                {params.regionmap_extent} {input.clusters} > {log} 2>&1
            """
