# ADAPTOGENE Pipeline - Refactored
# vim: filetype=python
import os
from pathlib import Path

#=============================================================================
# CONFIGURATION
#=============================================================================
#configfile: "/pipeline/config.yaml"

#=============================================================================
# HELPER FUNCTIONS
#=============================================================================
def check_numeric(value, name, allow_null=False):
    if allow_null and value is None: return True
    try: int(value); return True
    except: raise ValueError(f"{name} must be numeric, got: {value}")

def check_float(value, name, allow_null=False):
    if allow_null and value is None: return True
    try: float(value); return True
    except: raise ValueError(f"{name} must be float, got: {value}")

def check_in_list(value, allowed, name):
    if value not in allowed:
        raise ValueError(f"{name} must be one of {allowed}, got: {value}")

def check_file_exists(directory, filename, name):
    if not os.path.exists(os.path.join(directory, filename)):
        raise ValueError(f"File not found for {name}: {os.path.join(directory, filename)}")

def get_vcf_basename(vcf_path):
    """Extract basename from VCF file (handles .vcf and .vcf.gz)."""
    name = os.path.basename(vcf_path)
    if name.endswith('.vcf.gz'): return name[:-7]
    if name.endswith('.vcf'): return name[:-4]
    raise ValueError(f"VCF must end with .vcf or .vcf.gz: {vcf_path}")

def k_range(start, end):
    return list(range(int(start), int(end) + 1))

#=============================================================================
# PARSE AND VALIDATE CONFIGURATION
#=============================================================================
# Paths
INDIR = '/pipeline/' + config["INDIR"]
PROJECT = config["PROJECTNAME"]
OUTDIR = f'/pipeline/{PROJECT}_results/'
LOGDIR = f'/pipeline/{PROJECT}_logs/'

# Basic parameters
CPU = config['CPU']; check_numeric(CPU, 'CPU')
VCF_RAW = config['VCF_RAW']; check_file_exists(INDIR, VCF_RAW, 'VCF_RAW')
SAMPLES = config['SAMPLES']; check_file_exists(INDIR, SAMPLES, 'SAMPLES')
VCF_BASE = get_vcf_basename(VCF_RAW)

# VCF filtering
MAF = config['MAF']; check_float(MAF, 'MAF')
MISS = config['MISS']; check_float(MISS, 'MISS')
SAMPLE_MISS = config.get('SAMPLE_MISS', 0.5); check_float(SAMPLE_MISS, 'SAMPLE_MISS')
LD_WIN = config['LDwin']; check_numeric(LD_WIN, 'LDwin')
LD_STEP = config['LDstep']; check_numeric(LD_STEP, 'LDstep')
LD_R2 = config['LDr2']; check_float(LD_R2, 'LDr2')

# SNMF parameters
K_START = config['K_START']; check_numeric(K_START, 'K_START')
K_END = config['K_END']; check_numeric(K_END, 'K_END')
PLOIDY = config['PLOIDY']; check_numeric(PLOIDY, 'PLOIDY')
REPEAT = config['REPEAT']; check_numeric(REPEAT, 'REPEAT')
SNMF_PROJECT_MODE = config.get('PROJECT', 'new')

# structure_K parameters (used when mode=structure_K)
K_BEST = config.get('K_BEST', None)
PREDICTORS_SELECTED = config.get('PREDICTORS_SELECTED', '')
CROP_REGION = config.get('CROP_REGION', 'auto')
GAP = config.get('GAP', 0.5)
RESOLUTION = config.get('RESOLUTION', 0.5)
METRICS_WINSIZE = config.get('METRICS_WINSIZE', 10000)
CALC_POP_STATS = config.get('CALC_POP_STATS', False)  # Tajima's D, Pi, AMOVA, IBD - requires >= 3 samples per population
CUSTOM_TRAIT = config.get('CUSTOM_TRAIT', 'NULL')

# PieMap plot parameters (PIE_SIZE and PIE_RESCALE removed - now autoscaled based on map extent)
# Palette is now fixed to viridis plasma (colorblind-friendly)
PIEMAP_PIE_ALPHA = config.get('PieMap_PIE_ALPHA', 0.6)
PIEMAP_POP_LABEL = config.get('PieMap_POP_LABEL', 'F')
PIEMAP_POP_LABEL_SIZE = config.get('PieMap_POP_LABEL_SIZE', 10)

# Association parameters
GFF = config.get('GFF', '')
GFF_FEATURE = config.get('GFF_FEATURE', 'mRNA')
GENE_DISTANCE = config.get('GENE_DISTANCE', 1000000)
SNP_DISTANCE = config.get('SNP_DISTANCE', 100000)
PROMOTER_LENGTH = config.get('PROMOTER_LENGTH', 10000)
SIGSNPS_METHOD = config.get('sigSNPs_METHOD', 'EMMAX')
SIGSNPS_GAP = config.get('sigSNPs_GAP', 100000)
REGION_DISTANCE = config.get('REGION_DISTANCE', 2000000)
GO_FIELD = config.get('GO_FIELD', 'NULL')
TOP_REGIONS = config.get('TOP_REGIONS', 10)

# Regionplot parameters
GFF_GENENAME = config.get('GFF_GENENAME', 'description')
GFF_BIOTYPE = config.get('GFF_BIOTYPE', 'biotype')
GENES_TO_HIGHLIGHT = config.get('GENES_TO_HIGHLIGHT', 'all')
REGIONPLOT_REGION = config.get('REGIONPLOT_REGION', 'NULL')
REGIONPLOT_TRAITS = config.get('REGIONPLOT_TRAITS', 'NULL')
REGIONPLOT_ASSOCMETHOD = config.get('REGIONPLOT_ASSOCMETHOD', 'NULL')

# Maladaptation parameters
SSP = config.get('SSP', '585')
YEAR = config.get('YEAR', '2061-2080')
MODELS_STR = config.get('MODELS', '')
MODELS_LIST = [m.strip() for m in MODELS_STR.split(',') if m.strip()] if MODELS_STR else []
NTREE = config.get('NTREE', '500')
COR_THRESHOLD = config.get('COR_THRESHOLD', '0.5')
PCNM = config.get('PCNM', 'with')
GF_SUFFIX = config.get('GF_SUFFIX', '')
GF_RANDOM_MODEL = config.get('GF_RANDOM_MODEL', True)

def parse_association_configs(config):
    """Parse ASSOCIATION_CONFIGS into method -> adjust_threshold dict."""
    assoc_configs = config.get('ASSOCIATION_CONFIGS', [])
    configs = {}
    for cfg in assoc_configs:
        method = cfg['METHOD']
        adjust = f"{cfg['ADJUST']}_{cfg['THRESHOLD']}"
        if method in configs:
            raise ValueError(f"Method '{method}' appears multiple times in ASSOCIATION_CONFIGS")
        configs[method] = adjust
    return configs

ASSOC_CONFIGS = parse_association_configs(config)

#=============================================================================
# PATH DEFINITIONS
#=============================================================================
# Directory tags (easy to modify if adding new parameters)
FILT_TAG = f"maf{MAF}_miss{MISS}_smiss{SAMPLE_MISS}"
LD_TAG = f"ld{LD_R2}_win{LD_WIN}_step{LD_STEP}"
CLIMATE_TAG = f"climate_{RESOLUTION}"

# Base directories
WORK = f"{OUTDIR}work/"
PLOTS = f"{OUTDIR}plots/"
TABLES = f"{OUTDIR}tables/"
INTER = f"{OUTDIR}intermediate/"

# Working subdirectories (parameter-based structure)
WORK_FILT = f"{WORK}{FILT_TAG}/"
WORK_LD = f"{WORK_FILT}{LD_TAG}/"

# Working paths
W = {
    # Samples (in intermediate - shared across all filtering)
    'samples_list': f"{INTER}samples.list",
    'samples_filtered': f"{INTER}samples_filtered.list",  # After missingness filtering
    'samples_removed': f"{INTER}samples_removed.list",    # Samples removed by missingness
    'samples_missing_stats': f"{INTER}samples_missing_stats.tsv",  # Per-sample missing stats
    'samples_order': f"{INTER}samples_order.list",
    # Filtered VCF (in FILT_TAG directory)
    'vcf_filt': f"{WORK_FILT}{VCF_BASE}.vcf",
    # LD-pruned files (in LD_TAG subdirectory)
    'vcf_ld': f"{WORK_LD}{VCF_BASE}.vcf",
    'prune_in': f"{WORK_LD}{VCF_BASE}.prune.in",
    'geno': f"{WORK_LD}{VCF_BASE}.geno",
    'lfmm': f"{WORK_LD}{VCF_BASE}.lfmm",
    # LEA PCA outputs (created by pca_plot rule)
    # LEA::pca() strips extension and creates {basename}.pca/ directory
    'pca_projections': f"{WORK_LD}{VCF_BASE}.pca/{VCF_BASE}.projections",
    'pca_eigenvalues': f"{WORK_LD}{VCF_BASE}.pca/{VCF_BASE}.eigenvalues",
    'vcfsnp': f"{WORK_LD}{VCF_BASE}.vcfsnp",
    'removed': f"{WORK_LD}{VCF_BASE}.removed",
    'snmf': f"{WORK_LD}{VCF_BASE}.snmfProject",
}

# Output paths (organized results)
O = {
    'metadata': f"{TABLES}structure/metadata.tsv",
    'pca': f"{PLOTS}pca/pca.png",
    'pca_svg': f"{PLOTS}pca/pca.svg",
    'tracy': f"{PLOTS}pca/tracy_widom.png",
    'cross_entropy': f"{PLOTS}structure/cross_entropy_K{K_START}-{K_END}.png",
}

# K_BEST dependent paths (added dynamically when K_BEST is set)
def add_kbest_paths():
    """Add K_BEST dependent paths to W and O dictionaries."""
    if K_BEST is None:
        return
    
    # Imputed files in LD directory
    W['lfmm_imp'] = f"{WORK_LD}{VCF_BASE}_K{K_BEST}imp.lfmm"
    W['vcf_imp'] = f"{WORK_LD}{VCF_BASE}_K{K_BEST}imp.vcf"
    
    # Climate data
    W['climate_raster'] = f"{INTER}Climate_present_RasterStack.grd"

    # Tables - climate (resolution-specific)
    O['climate_site'] = f"{TABLES}{CLIMATE_TAG}/Climate_present_site.tsv"
    O['climate_site_scaled'] = f"{TABLES}{CLIMATE_TAG}/Climate_present_site_scaled.tsv"
    O['climate_all'] = f"{TABLES}{CLIMATE_TAG}/Climate_present_all.tsv"
    # Tables - structure/population stats
    O['tajima'] = f"{TABLES}structure/TajimaD_byPop.tsv"
    O['pi_div'] = f"{TABLES}structure/Pi_diversity_byPop.tsv"
    O['ibd_raw'] = f"{TABLES}structure/IBD_raw.tsv"
    O['ibd_pairs'] = f"{TABLES}structure/IBD_notIsolated.tsv"
    O['amova'] = f"{TABLES}structure/AMOVA.tsv"
    
    # Plots
    O['corr_heatmap'] = f"{PLOTS}{CLIMATE_TAG}/CorrelationHeatmap_present.png"
    O['mantel'] = f"{PLOTS}structure/MantelTest.png"
    O['amova_plot'] = f"{PLOTS}structure/AMOVA.png"

add_kbest_paths()

# Association paths (added when K_BEST is set and association mode is used)
def add_association_paths():
    """Add association-specific paths to W and O dictionaries."""
    if K_BEST is None or not ASSOC_CONFIGS:
        return

    # Full (non-LD pruned) files for association analysis
    W['geno_full'] = f"{WORK_FILT}{VCF_BASE}.geno"
    W['lfmm_full'] = f"{WORK_FILT}{VCF_BASE}.lfmm"
    W['vcfsnp_full'] = f"{WORK_FILT}{VCF_BASE}.vcfsnp"
    W['removed_full'] = f"{WORK_FILT}{VCF_BASE}.removed"
    W['snmf_full'] = f"{WORK_FILT}{VCF_BASE}.snmfProject"
    W['lfmm_imp_full'] = f"{WORK_FILT}{VCF_BASE}_K{K_BEST}imp.lfmm"
    W['vcf_imp_full'] = f"{WORK_FILT}{VCF_BASE}_K{K_BEST}imp.vcf"

    # EMMAX work directory (kinship, tped, tfam files)
    W['emmax_work'] = f"{WORK_FILT}emmax/"

    # Combined outputs - association
    O['selected_snps'] = f"{TABLES}association/Selected_SNPs.tsv"
    O['regions'] = f"{TABLES}association/Regions.tsv"
    O['genes_per_region'] = f"{TABLES}association/Genes_per_region.tsv"
    O['genes_per_region_collapsed'] = f"{TABLES}association/Genes_per_region_collapsed.tsv"
    O['enrichment'] = f"{TABLES}association/enrichment/Enrichment_combined.tsv"

    # Regionplot outputs
    O['gff_topr'] = f"{INTER}topr_gene_annotation.tsv"
    O['regionplot_done'] = f"{PLOTS}regionplot/.done"

add_association_paths()

# Maladaptation paths
def add_maladaptation_paths():
    """Add maladaptation-specific paths to W and O dictionaries."""
    if K_BEST is None or not MODELS_LIST:
        return

    SUFFIX = f"{GF_SUFFIX}_{PCNM}PCNM"

    # Per-model future climate rasters
    for model in MODELS_LIST:
        W[f'climate_future_{model}'] = f"{INTER}Climate_future_year{YEAR}_ssp{SSP}_{model}.grd"

    # Merged future climate
    W['climate_future_raster'] = f"{INTER}Climate_future_year{YEAR}_ssp{SSP}_RasterStack.grd"
    O['climate_future_all'] = f"{TABLES}{CLIMATE_TAG}/Climate_future_year{YEAR}_ssp{SSP}_all.tsv"
    O['climate_future_site'] = f"{TABLES}{CLIMATE_TAG}/Climate_future_year{YEAR}_ssp{SSP}_site.tsv"

    # Gradient Forest models
    W['gf_adaptive'] = f"{INTER}gradientForest_adaptive_{SUFFIX}.qs"
    W['gf_random'] = f"{INTER}gradientForest_random_{SUFFIX}.qs"

    # Genetic offset outputs
    W['gf_offset_raster'] = f"{INTER}GeneticOffset_{SUFFIX}.grd"
    O['gf_offset_map_values'] = f"{TABLES}gradientForest/GeneticOffset_map_{SUFFIX}.tsv"
    O['gf_offset_site_values'] = f"{TABLES}gradientForest/GeneticOffset_site_{SUFFIX}.tsv"

    # Gradient Forest plots
    O['gf_offset_piemap'] = f"{PLOTS}gradientForest/GeneticOffsetPieMap_{SUFFIX}.png"
    O['gf_offset_piemap_tajima'] = f"{PLOTS}gradientForest/GeneticOffsetPieMap_{SUFFIX}_TajimaD.png"
    O['gf_offset_piemap_diversity'] = f"{PLOTS}gradientForest/GeneticOffsetPieMap_{SUFFIX}_PiDiversity.png"
    O['gf_cumimp'] = f"{PLOTS}gradientForest/CumulativeImportance_{SUFFIX}.png"
    O['gf_importance'] = f"{PLOTS}gradientForest/OverallImportance_{SUFFIX}.png"

    # Future climate density plot
    O['density_future'] = f"{PLOTS}{CLIMATE_TAG}/DensityPlot_future_ssp{SSP}_{YEAR}.png"

add_maladaptation_paths()

# Templates for K-dependent outputs
def clusters_table(k): return f"{TABLES}structure/clusters_K{k}.tsv"
def structure_plot(k): return f"{PLOTS}structure/structure_K{k}.png"
def pca_struct_plot(k): return f"{PLOTS}pca/pca_structure_K{k}.png"
def pop_diff_plot(k): return f"{PLOTS}structure/pop_diff_K{k}.png"

# Templates for climate/trait-dependent outputs
DENSITY_PLOT_COMBINED = f"{PLOTS}{CLIMATE_TAG}/DensityPlot_present.png"
def piemap_tajima(bio): return f"{PLOTS}piemap/PieMap_{bio}_TajimaD.png"
def piemap_diversity(bio): return f"{PLOTS}piemap/PieMap_{bio}_PiDiversity.png"
def piemap_notrait(bio): return f"{PLOTS}piemap/PieMap_{bio}.png"

# Templates for association outputs
def assoc_pvalues(method): return f"{TABLES}association/{method}/{method}_pvalues_K{K_BEST}.tsv"
def assoc_sigsnps(method, adjust): return f"{TABLES}association/{method}/{method}_pvalues_K{K_BEST}_sigSNPs_{adjust}.tsv"
def assoc_genes(method, adjust): return f"{TABLES}association/{method}/{method}_pvalues_K{K_BEST}_{adjust}_genesAround{GENE_DISTANCE}.tsv"
def assoc_genes_collapsed(method, adjust): return f"{TABLES}association/{method}/{method}_pvalues_K{K_BEST}_{adjust}_genesAround{GENE_DISTANCE}_collapsed.tsv"
def manhattan_plot(method, trait, adjust): return f"{PLOTS}{method}/Manhattan_{trait}_K{K_BEST}_{adjust}.png"
def manhattan_plot_regions(method, trait, adjust): return f"{PLOTS}{method}/Manhattan_{trait}_K{K_BEST}_{adjust}_regions.png"

#=============================================================================
# CREATE DIRECTORIES
#=============================================================================
dirs_to_create = [WORK, WORK_FILT, WORK_LD, PLOTS,
          f"{PLOTS}pca/", f"{PLOTS}structure/", f"{PLOTS}{CLIMATE_TAG}/", f"{PLOTS}piemap/",
          TABLES, f"{TABLES}structure/", f"{TABLES}{CLIMATE_TAG}/",
          INTER, LOGDIR]

# Add association directories for each method
for method in ASSOC_CONFIGS:
    dirs_to_create.append(f"{PLOTS}{method}/")
    dirs_to_create.append(f"{TABLES}association/{method}/")

# Add enrichment directory
dirs_to_create.append(f"{TABLES}association/")
dirs_to_create.append(f"{TABLES}association/enrichment/")

# Add regionplot directory
dirs_to_create.append(f"{PLOTS}regionplot/")

# Add gradientForest directories
dirs_to_create.append(f"{PLOTS}gradientForest/")
dirs_to_create.append(f"{TABLES}gradientForest/")

for d in dirs_to_create:
    os.makedirs(d, exist_ok=True)

workdir: OUTDIR

#=============================================================================
# MODE AND TARGETS
#=============================================================================
MODE = config.get('mode', None)

def get_predictors_list():
    """Parse PREDICTORS_SELECTED into list."""
    if not PREDICTORS_SELECTED:
        return []
    return [p.strip() for p in PREDICTORS_SELECTED.split(',')]

def get_targets(mode):
    if mode == 'processing':
        return [
            W['samples_missing_stats'], W['samples_removed'],  # Sample missingness outputs
            W['vcf_filt'], W['vcf_ld'], W['geno'], W['lfmm'], O['metadata']
        ]
    
    elif mode == 'structure':
        ks = k_range(K_START, K_END)
        return (
            [clusters_table(k) for k in ks] +
            [structure_plot(k) for k in ks] +
            [pca_struct_plot(k) for k in ks] +
            [pop_diff_plot(k) for k in ks] +
            [O['pca'], O['tracy'], O['cross_entropy']]
        )
    
    elif mode == 'structure_K':
        check_numeric(K_BEST, 'K_BEST')
        predictors = get_predictors_list()
        if not predictors:
            raise ValueError("PREDICTORS_SELECTED must be set for structure_K mode")

        targets = (
            # Imputed data
            [W['lfmm_imp'], W['vcf_imp']] +
            # Climate data
            [O['climate_site'], O['climate_site_scaled'], O['climate_all'], W['climate_raster']] +
            # Combined density plot for all predictors
            [DENSITY_PLOT_COMBINED] +
            # Correlation heatmap (doesn't require pop stats)
            [O['corr_heatmap']]
        )

        # Simple PieMaps (uniform pie size) - always generated
        targets += [piemap_notrait(bio) for bio in predictors]

        # Population statistics (requires >= 3 samples per population)
        if CALC_POP_STATS:
            targets += (
                # Population metrics
                [O['tajima'], O['pi_div'], O['ibd_raw'], O['ibd_pairs']] +
                # Summary plots requiring pop stats
                [O['mantel'], O['amova'], O['amova_plot']] +
                # PieMaps with trait overlays (in addition to simple piemaps)
                [piemap_tajima(bio) for bio in predictors] +
                [piemap_diversity(bio) for bio in predictors]
            )

        return targets

    elif mode == 'association':
        check_numeric(K_BEST, 'K_BEST')
        predictors = get_predictors_list()
        if not predictors:
            raise ValueError("PREDICTORS_SELECTED must be set for association mode")
        if not ASSOC_CONFIGS:
            raise ValueError("ASSOCIATION_CONFIGS must be set for association mode")
        if GFF:
            check_file_exists(INDIR, GFF, 'GFF')

        targets = []

        # Per-method targets
        for method, adjust in ASSOC_CONFIGS.items():
            # P-values table
            targets.append(assoc_pvalues(method))
            # Significant SNPs per method
            targets.append(assoc_sigsnps(method, adjust))
            # Manhattan plots per trait (simple + with regions)
            # Each rule produces both PNG and SVG
            for trait in predictors:
                targets.append(manhattan_plot(method, trait, adjust))
                targets.append(manhattan_plot_regions(method, trait, adjust))

        # Combined analysis targets
        targets.append(O['selected_snps'])
        targets.append(O['regions'])
        targets.append(O['genes_per_region'])
        targets.append(O['genes_per_region_collapsed'])

        # Enrichment (if GO_FIELD is specified)
        if GO_FIELD and GO_FIELD != 'NULL':
            targets.append(O['enrichment'])

        return targets

    elif mode == 'regionplot':
        check_numeric(K_BEST, 'K_BEST')
        if not ASSOC_CONFIGS:
            raise ValueError("ASSOCIATION_CONFIGS must be set for regionplot mode")
        if GFF:
            check_file_exists(INDIR, GFF, 'GFF')

        return [O['regionplot_done']]

    elif mode == 'maladaptation':
        check_numeric(K_BEST, 'K_BEST')
        check_numeric(SSP, 'SSP')
        check_numeric(NTREE, 'NTREE')
        check_float(COR_THRESHOLD, 'COR_THRESHOLD')
        if not MODELS_LIST:
            raise ValueError("MODELS must be set for maladaptation mode")

        targets = [
            # Future climate
            O['climate_future_site'],
            O['climate_future_all'],
            W['climate_future_raster'],
            # Gradient Forest models
            W['gf_adaptive'],
            # Genetic offset
            O['gf_offset_map_values'],
            O['gf_offset_site_values'],
            # Plots
            O['gf_cumimp'],
            O['gf_importance'],
            O['gf_offset_piemap'],
            O['density_future'],
        ]

        if GF_RANDOM_MODEL:
            targets.append(W['gf_random'])

        # Add TajimaD and PiDiversity GO piemap variants if population stats were calculated
        if CALC_POP_STATS:
            targets.extend([
                O['gf_offset_piemap_tajima'],
                O['gf_offset_piemap_diversity'],
            ])

        return targets

    elif mode is None:
        raise ValueError("Specify mode: --config mode=processing or mode=structure or mode=structure_K or mode=association")
    else:
        raise ValueError(f"Unknown mode: {mode}")

#=============================================================================
# MAIN RULE
#=============================================================================
rule all:
    input: get_targets(MODE)

#=============================================================================
# MODULE 1: VCF PROCESSING
#=============================================================================

rule extract_samples:
    """Extract sample IDs from metadata for VCF subsetting."""
    input:  samples = f"{INDIR}{SAMPLES}"
    output: W['samples_list']
    log:    f"{LOGDIR}extract_samples.log"
    shell:  "tail -n +2 {input.samples} | awk '{{print 0, $2}}' > {output} 2> {log}"

rule calculate_sample_missing:
    """Calculate per-sample missing genotype rate and filter samples."""
    input:
        vcf = f"{INDIR}{VCF_RAW}",
        samples = W['samples_list']
    output:
        stats = W['samples_missing_stats'],
        filtered = W['samples_filtered'],
        removed = W['samples_removed']
    params:
        threshold = SAMPLE_MISS,
        prefix = f"{INTER}sample_miss_tmp"
    log: f"{LOGDIR}calculate_sample_missing.log"
    threads: CPU
    shell:
        """
        # Calculate per-sample missingness using plink
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# --keep {input.samples} \
            --missing --out {params.prefix} > {log} 2>&1

        # Create stats file with header
        echo -e "FID\\tIID\\tMISS_PHENO\\tN_MISS\\tN_GENO\\tF_MISS" > {output.stats}
        tail -n +2 {params.prefix}.imiss >> {output.stats}

        # Filter samples: keep those with F_MISS <= threshold
        awk -v thresh={params.threshold} 'NR>1 && $6 <= thresh {{print $1, $2}}' {params.prefix}.imiss > {output.filtered}

        # List removed samples: those with F_MISS > threshold
        awk -v thresh={params.threshold} 'NR>1 && $6 > thresh {{print $1, $2, $6}}' {params.prefix}.imiss > {output.removed}

        # Log summary
        n_total=$(wc -l < {input.samples})
        n_kept=$(wc -l < {output.filtered})
        n_removed=$(wc -l < {output.removed})
        echo "INFO: Sample missingness filtering (threshold: {params.threshold})" >> {log}
        echo "INFO: Total samples: $n_total" >> {log}
        echo "INFO: Samples passing: $n_kept" >> {log}
        echo "INFO: Samples removed: $n_removed" >> {log}

        # Cleanup temp files
        rm -f {params.prefix}.*
        """

rule filter_vcf:
    """Filter VCF by MAF, missingness, and sample list (after sample missingness filter)."""
    input:
        vcf = f"{INDIR}{VCF_RAW}",
        samples = W['samples_filtered']  # Use filtered samples list
    output: W['vcf_filt']
    params: prefix = W['vcf_filt'].replace('.vcf', ''), maf = MAF, miss = MISS
    log:    f"{LOGDIR}filter_vcf.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# --keep {input.samples} \
            --maf {params.maf} --geno {params.miss} \
            --recode vcf --out {params.prefix} > {log} 2>&1
        sed -i '/^#CHROM/s/\\t0_/\\t/g' {output}
        """

rule ld_prune:
    """LD prune VCF (PCA is done separately via LEA)."""
    input:  vcf = W['vcf_filt']
    output: vcf = W['vcf_ld'], prune = W['prune_in']
    params: prefix = W['vcf_ld'].replace('.vcf', ''), win = LD_WIN, step = LD_STEP, r2 = LD_R2
    log:    f"{LOGDIR}ld_prune.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# \
            --indep-pairwise {params.win} {params.step} {params.r2} \
            --out {params.prefix} > {log} 2>&1

        plink --vcf {input.vcf} --const-fid --allow-extra-chr \
            --set-missing-var-ids @:# --extract {output.prune} \
            --make-bed --recode vcf \
            --out {params.prefix} >> {log} 2>&1

        sed -i '/^#CHROM/s/\\t0_0_/\\t/g' {output.vcf}
        """

rule extract_vcf_sample_order:
    """Get sample order from VCF header for metadata alignment."""
    input:  vcf = W['vcf_filt']
    output: W['samples_order']
    log:    f"{LOGDIR}extract_vcf_sample_order.log"
    shell:  "grep -m1 CHROM {input.vcf} | cut -f10- | tr '\\t' '\\n' > {output} 2> {log}"

rule align_metadata:
    """Align metadata rows to match VCF sample order."""
    input:  meta = f"{INDIR}{SAMPLES}", order = W['samples_order']
    output: O['metadata']
    log:    f"{LOGDIR}align_metadata.log"
    shell:  "Rscript /pipeline/scripts/filter_arrange_metadata.R {input.meta} {input.order} {output} > {log} 2>&1"

rule vcf_to_lfmm:
    """Convert VCF to LEA formats (geno, lfmm)."""
    input:  vcf = W['vcf_ld']
    output: geno = W['geno'], lfmm = W['lfmm'], vcfsnp = W['vcfsnp'], removed = W['removed']
    log:    f"{LOGDIR}vcf_to_lfmm.log"
    shell:  "Rscript /pipeline/scripts/vcf2lfmm.R {input.vcf} > {log} 2>&1"

#=============================================================================
# MODULE 2: POPULATION STRUCTURE
#=============================================================================

rule snmf:
    """Run sNMF across K range."""
    input:  geno = W['geno']
    output: W['snmf']
    params: ks = K_START, ke = K_END, ploidy = PLOIDY, rep = REPEAT, mode = SNMF_PROJECT_MODE
    log:    f"{LOGDIR}snmf.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/SNMF.R {input.geno} \
            {params.ks} {params.ke} {params.ploidy} {params.rep} \
            {threads} {params.mode} > {log} 2>&1
        """

rule pca_plot:
    """Generate PCA and Tracy-Widom plots (also creates LEA PCA projections/eigenvalues)."""
    input:  lfmm = W['lfmm'], meta = O['metadata']
    output:
        pca = O['pca'],
        pca_svg = O['pca_svg'],
        tracy = O['tracy'],
        projections = W['pca_projections'],
        eigenvalues = W['pca_eigenvalues']
    params: plot_dir = f"{PLOTS}pca/", inter_dir = INTER
    log:    f"{LOGDIR}pca_plot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pca.R \
            {input.lfmm} {input.meta} {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule cross_entropy_plot:
    """Plot cross-entropy for K selection."""
    input:  snmf = W['snmf']
    output: O['cross_entropy']
    params: ks = K_START, ke = K_END, plot_dir = f"{PLOTS}structure/", inter_dir = INTER
    log:    f"{LOGDIR}cross_entropy.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_cross_entropy.R \
            {input.snmf} {params.ks} {params.ke} \
            {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule extract_clusters:
    """Extract Q-matrix (cluster assignments) for specific K."""
    input:  snmf = W['snmf'], meta = O['metadata']
    output: clusters_table("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k
    log:    f"{LOGDIR}extract_clusters_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/extract_clusters.R \
            {input.snmf} {input.meta} {params.k} {output} > {log} 2>&1
        """

rule structure_barplot:
    """Generate structure barplot for specific K."""
    input:  clusters = clusters_table("{k}")
    output: structure_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k, plot_dir = f"{PLOTS}structure/", inter_dir = INTER
    log:    f"{LOGDIR}structure_plot_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_structure.R \
            {input.clusters} {params.k} {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule pca_structure_plot:
    """Generate PCA with structure pie charts for specific K (uses LEA PCA)."""
    input:
        clusters = clusters_table("{k}"),
        projections = W['pca_projections'],
        eigenvalues = W['pca_eigenvalues']
    output: pca_struct_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k, plot_dir = f"{PLOTS}pca/", inter_dir = INTER
    log:    f"{LOGDIR}pca_structure_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pca_structure.R \
            {input.clusters} {input.projections} {input.eigenvalues} \
            {params.k} {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule pop_diff_test:
    """Population differentiation test for specific K."""
    input:  snmf = W['snmf']
    output: pop_diff_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k, ploidy = PLOIDY, plot_dir = f"{PLOTS}structure/", inter_dir = INTER
    log:    f"{LOGDIR}pop_diff_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pop_diff.R \
            {input.snmf} {params.k} {params.ploidy} \
            {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

#=============================================================================
# MODULE 3: STRUCTURE_K (requires K_BEST selection)
#=============================================================================

rule impute_ld:
    """Impute missing genotypes using SNMF results for K_BEST."""
    input:  snmf = W['snmf'], lfmm = W['lfmm']
    output: W['lfmm_imp']
    params: k = K_BEST
    log:    f"{LOGDIR}impute_ld.log"
    shell:
        """
        Rscript /pipeline/scripts/Impute.R {input.snmf} {input.lfmm} {params.k} {output} > {log} 2>&1
        """

rule lfmm2vcf_ld:
    """Convert imputed LFMM back to VCF format."""
    input:  vcf = W['vcf_ld'], lfmm_imp = W['lfmm_imp']
    output: W['vcf_imp']
    params: blocksize = 20000
    log:    f"{LOGDIR}lfmm2vcf_ld.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm2vcf.R {input.lfmm_imp} {input.vcf} {params.blocksize} > {log} 2>&1
        """
    # NOTE: lfmm2vcf.R already exists, outputs .vcf with same base as .lfmm input

rule download_climate_present:
    """Download and process present climate data for sampling locations."""
    input:  meta = O['metadata']
    output:
        site = O['climate_site'],
        site_scaled = O['climate_site_scaled'],
        all_vals = O['climate_all'],
        raster = W['climate_raster']
    params:
        crop = CROP_REGION,
        gap = GAP,
        resolution = RESOLUTION,
        data_dir = f"{INDIR}",
        inter_dir = INTER,
        tables_dir = f"{TABLES}{CLIMATE_TAG}/"
    log: f"{LOGDIR}download_climate_present.log"
    shell:
        """
        Rscript /pipeline/scripts/download_climate_present.R \
            {input.meta} {params.crop} {params.gap} {params.data_dir} {params.resolution} \
            {params.inter_dir} {params.tables_dir} > {log} 2>&1
        """

rule density_plot:
    """Generate combined density plot for all climate predictors."""
    input:  climate = O['climate_site']
    output: DENSITY_PLOT_COMBINED
    params:
        predictors = PREDICTORS_SELECTED,
        plot_dir = f"{PLOTS}{CLIMATE_TAG}/",
        inter_dir = INTER,
        prefix = "DensityPlot_present"
    log:    f"{LOGDIR}density_plot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_density.R \
            {input.climate} {params.predictors} {params.plot_dir} {params.inter_dir} \
            {params.prefix} > {log} 2>&1
        """

rule tajima_d:
    """Calculate Tajima's D per population."""
    input:  vcf = W['vcf_filt'], meta = O['metadata']
    output: O['tajima']
    params: winsize = METRICS_WINSIZE, inter_dir = INTER
    log:    f"{LOGDIR}tajima_d.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/TajimaD.R \
            {input.vcf} {input.meta} {params.winsize} {threads} \
            {output} {params.inter_dir} > {log} 2>&1
        """

rule pi_diversity:
    """Calculate nucleotide diversity (Pi) per population."""
    input:  vcf = W['vcf_filt'], meta = O['metadata']
    output: O['pi_div']
    params: winsize = METRICS_WINSIZE, inter_dir = INTER
    log:    f"{LOGDIR}pi_diversity.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/Pi_diversity.R \
            {input.vcf} {input.meta} {params.winsize} {threads} \
            {output} {params.inter_dir} > {log} 2>&1
        """

rule ibd:
    """Calculate isolation by distance between populations."""
    input:  clusters = clusters_table(K_BEST), meta = O['metadata']
    output: raw = O['ibd_raw'], pairs = O['ibd_pairs']
    params: tables_dir = f"{TABLES}structure/"
    log:    f"{LOGDIR}ibd.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/IBD.R \
            {input.clusters} {input.meta} {threads} {params.tables_dir} > {log} 2>&1
        """

rule correlation_heatmap:
    """Generate correlation heatmap of climate variables and traits."""
    input:  climate = O['climate_site'], meta = O['metadata']
    output: O['corr_heatmap']
    params: plot_dir = f"{PLOTS}{CLIMATE_TAG}/", inter_dir = INTER
    log:    f"{LOGDIR}correlation_heatmap.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_correlation_heatmap.R \
            {input.climate} {input.meta} {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule mantel_test:
    """Perform Mantel test for IBD/IBE."""
    input:  meta = O['metadata'], climate = O['climate_site'], clusters = clusters_table(K_BEST)
    output: O['mantel']
    params:
        predictors = PREDICTORS_SELECTED,
        plot_dir = f"{PLOTS}structure/",
        inter_dir = INTER
    log: f"{LOGDIR}mantel_test.log"
    shell:
        """
        Rscript /pipeline/scripts/mantel_test.R \
            {input.meta} {input.clusters} {input.climate} {params.predictors} \
            {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule amova:
    """Perform AMOVA analysis."""
    input:  vcf = W['vcf_imp'], meta = O['metadata']
    output: table = O['amova'], plot = O['amova_plot']
    params: plot_dir = f"{PLOTS}structure/", tables_dir = f"{TABLES}structure/", inter_dir = INTER
    log:    f"{LOGDIR}amova.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/amova.R \
            {input.vcf} {input.meta} {threads} \
            {params.plot_dir} {params.tables_dir} {params.inter_dir} > {log} 2>&1
        """

rule piemap_plot:
    """Generate PieMap plots for a climate predictor with trait-scaled pie sizes."""
    input:
        meta = O['metadata'],
        clusters = clusters_table(K_BEST),
        raster = W['climate_raster'],
        tajima = O['tajima'],
        diversity = O['pi_div']
    output:
        tajima_plot = piemap_tajima("{bio}"),
        diversity_plot = piemap_diversity("{bio}")
    wildcard_constraints: bio = r"bio_\d+"
    params:
        bio = lambda wc: wc.bio,
        pie_alpha = PIEMAP_PIE_ALPHA,
        pop_label = PIEMAP_POP_LABEL,
        pop_label_size = PIEMAP_POP_LABEL_SIZE,
        plot_dir = f"{PLOTS}piemap/",
        inter_dir = INTER
    log: f"{LOGDIR}piemap_{{bio}}.log"
    shell:
        """
        # TajimaD piemap
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            {input.tajima} "Tajima's D" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            PieMap_{params.bio}_TajimaD > {log} 2>&1

        # PiDiversity piemap
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            {input.diversity} "Pi Diversity" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            PieMap_{params.bio}_PiDiversity >> {log} 2>&1
        """

rule piemap_simple:
    """Generate basic PieMap plots with autoscaled uniform pie size (no trait overlay)."""
    input:
        meta = O['metadata'],
        clusters = clusters_table(K_BEST),
        raster = W['climate_raster']
    output: piemap_notrait("{bio}")
    wildcard_constraints: bio = r"bio_\d+"
    params:
        bio = lambda wc: wc.bio,
        pie_alpha = PIEMAP_PIE_ALPHA,
        pop_label = PIEMAP_POP_LABEL,
        pop_label_size = PIEMAP_POP_LABEL_SIZE,
        plot_dir = f"{PLOTS}piemap/",
        inter_dir = INTER
    log: f"{LOGDIR}piemap_simple_{{bio}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            NULL NULL \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            PieMap_{params.bio} > {log} 2>&1
        """

#=============================================================================
# MODULE 4: ASSOCIATION (EMMAX/LFMM)
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
    log: f"{LOGDIR}vcf_to_lfmm_full.log"
    shell: "Rscript /pipeline/scripts/vcf2lfmm.R {input.vcf} > {log} 2>&1"

rule snmf_full:
    """Run sNMF on full (non-LD pruned) dataset for imputation."""
    input:  geno = W['geno_full']
    output: W['snmf_full']
    params: ks = K_START, ke = K_END, ploidy = PLOIDY, rep = REPEAT, mode = SNMF_PROJECT_MODE
    log:    f"{LOGDIR}snmf_full.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/SNMF.R {input.geno} \
            {params.ks} {params.ke} {params.ploidy} {params.rep} \
            {threads} {params.mode} > {log} 2>&1
        """

rule impute_full:
    """Impute missing genotypes in full dataset using SNMF."""
    input:  snmf = W['snmf_full'], lfmm = W['lfmm_full']
    output: W['lfmm_imp_full']
    params: k = K_BEST
    log:    f"{LOGDIR}impute_full.log"
    shell:
        """
        Rscript /pipeline/scripts/Impute.R {input.snmf} {input.lfmm} {params.k} {output} > {log} 2>&1
        """

rule lfmm2vcf_full:
    """Convert imputed full LFMM to VCF format."""
    input:  vcf = W['vcf_filt'], lfmm_imp = W['lfmm_imp_full']
    output: W['vcf_imp_full']
    params: blocksize = 20000
    log:    f"{LOGDIR}lfmm2vcf_full.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm2vcf.R {input.lfmm_imp} {input.vcf} {params.blocksize} > {log} 2>&1
        """

# EMMAX analysis
rule emmax_analysis:
    """Run EMMAX association analysis for all traits."""
    input:
        vcf = W['vcf_filt'],
        traits = O['climate_site_scaled'],
        covariates = W['pca_projections'],  # LEA PCA projections
        metadata = O['metadata']
    output: assoc_pvalues("EMMAX")
    params:
        k = K_BEST,
        predictors = PREDICTORS_SELECTED,
        inter_dir = INTER,
        tables_dir = f"{TABLES}association/EMMAX/",
        emmax_work = W['emmax_work']
    log: f"{LOGDIR}emmax_analysis.log"
    shell:
        """
        Rscript /pipeline/scripts/emmax.R \
            {input.vcf} {params.k} {input.traits} {input.covariates} \
            {params.predictors} {params.inter_dir} {input.metadata} \
            {params.tables_dir} {params.emmax_work} > {log} 2>&1
        """

# LFMM analysis
rule lfmm_analysis:
    """Run LFMM association analysis for all traits."""
    input:
        lfmm_ld = W['lfmm_imp'],       # LD-pruned imputed (for model training)
        lfmm_full = W['lfmm_imp_full'], # Full imputed (for testing)
        climate = O['climate_site_scaled'],
        vcfsnp = W['vcfsnp_full']
    output: assoc_pvalues("LFMM")
    params:
        k = K_BEST,
        predictors = PREDICTORS_SELECTED,
        tables_dir = f"{TABLES}association/LFMM/"
    log: f"{LOGDIR}lfmm_analysis.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm.R \
            {input.lfmm_ld} {input.lfmm_full} {input.climate} \
            {params.k} {params.predictors} {input.vcfsnp} \
            {params.tables_dir} > {log} 2>&1
        """

# Manhattan plots - simple version (runs early, no region dependency)
# Produces both PNG and SVG in a single run
rule manhattan_plot:
    """Generate simple Manhattan plot for a specific trait and method."""
    input: assoc = lambda wc: assoc_pvalues(wc.method)
    output:
        png = f"{PLOTS}{{method}}/Manhattan_{{trait}}_K{K_BEST}_{{adjust}}.png",
        svg = f"{PLOTS}{{method}}/Manhattan_{{trait}}_K{K_BEST}_{{adjust}}.svg"
    wildcard_constraints:
        method = r"EMMAX|LFMM",
        trait = r"bio_\d+",
        adjust = r"\w+_[\d.]+"
    params:
        k = K_BEST,
        plot_dir = lambda wc: f"{PLOTS}{wc.method}/",
        regions = "NULL",  # No regions for simple plot
        selected_snps = "NULL"  # No selected SNPs for simple plot
    log: f"{LOGDIR}manhattan_{{method}}_{{trait}}_{{adjust}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_manhattan.R \
            {input.assoc} {wildcards.adjust} {params.k} {wildcards.method} \
            {wildcards.trait} {params.plot_dir} {params.regions} {params.selected_snps} > {log} 2>&1
        """

# Manhattan plots with regions highlighted (runs after regions are created)
# Produces both PNG and SVG in a single run
# For combined methods (Sum/Overlap/PairOverlap), shows all selected SNPs
# with different shapes: circle=current method, triangle=other method, diamond=both
rule manhattan_plot_regions:
    """Generate Manhattan plot with significant regions highlighted."""
    input:
        assoc = lambda wc: assoc_pvalues(wc.method),
        regions = O['regions'],
        # All sigSNPs files from all methods for per-trait method attribution
        sigsnps = lambda wc: [assoc_sigsnps(method, adjust) for method, adjust in ASSOC_CONFIGS.items()]
    output:
        png = f"{PLOTS}{{method}}/Manhattan_{{trait}}_K{K_BEST}_{{adjust}}_regions.png",
        svg = f"{PLOTS}{{method}}/Manhattan_{{trait}}_K{K_BEST}_{{adjust}}_regions.svg"
    wildcard_constraints:
        method = r"EMMAX|LFMM",
        trait = r"bio_\d+",
        adjust = r"\w+_[\d.]+"
    params:
        k = K_BEST,
        plot_dir = lambda wc: f"{PLOTS}{wc.method}/",
        sigsnps_str = lambda wc, input: ','.join(input.sigsnps)
    log: f"{LOGDIR}manhattan_regions_{{method}}_{{trait}}_{{adjust}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_manhattan.R \
            {input.assoc} {wildcards.adjust} {params.k} {wildcards.method} \
            {wildcards.trait} {params.plot_dir} {input.regions} "{params.sigsnps_str}" > {log} 2>&1
        """

# Find significant SNPs - one rule per method
rule find_sig_snps:
    """Find significant SNPs for a specific method."""
    input: assoc = lambda wc: assoc_pvalues(wc.method)
    output: f"{TABLES}association/{{method}}/{{method}}_pvalues_K{K_BEST}_sigSNPs_{{adjust}}.tsv"
    wildcard_constraints:
        method = r"EMMAX|LFMM",
        adjust = r"\w+_[\d.]+"
    params: snp_dist = SNP_DISTANCE
    log: f"{LOGDIR}find_sig_snps_{{method}}_{{adjust}}.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/find_sig_snps.R \
            {input.assoc} {wildcards.adjust} {params.snp_dist} \
            {wildcards.method} {threads} {output} > {log} 2>&1
        """

# Combine significant SNPs from different methods into Selected_SNPs table
rule combine_selected_snps:
    """Combine significant SNPs from all association methods."""
    input:
        sigsnps = lambda wc: [assoc_sigsnps(method, adjust) for method, adjust in ASSOC_CONFIGS.items()]
    output: O['selected_snps']
    params:
        sigsnps_str = lambda wc, input: ' '.join(input.sigsnps),
        method = SIGSNPS_METHOD,
        gap = SIGSNPS_GAP,
        predictors = PREDICTORS_SELECTED
    log: f"{LOGDIR}combine_selected_snps.log"
    shell:
        """
        Rscript /pipeline/scripts/combine_selected_snps.R \
            "{params.sigsnps_str}" {params.method} {params.gap} \
            {params.predictors} {output} > {log} 2>&1
        """

# Merge SNPs into regions based on distance
rule create_regions:
    """Merge nearby significant SNPs into regions."""
    input: selected_snps = O['selected_snps']
    output: O['regions']
    params: region_dist = REGION_DISTANCE
    log: f"{LOGDIR}create_regions.log"
    shell:
        """
        Rscript /pipeline/scripts/create_regions.R \
            {input.selected_snps} {params.region_dist} {output} > {log} 2>&1
        """

# Find genes around regions
rule find_genes_around_regions:
    """Find genes within GENE_DISTANCE of significant regions."""
    input:
        regions = O['regions'],
        gff = f"{INDIR}{GFF}",
        vcfsnp = W['vcfsnp_full']
    output:
        genes = O['genes_per_region'],
        collapsed = O['genes_per_region_collapsed']
    params:
        feature = GFF_FEATURE,
        distance = GENE_DISTANCE,
        promoter_len = PROMOTER_LENGTH,
        top_regions = TOP_REGIONS
    log: f"{LOGDIR}find_genes_around_regions.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/find_genes_around_regions.R \
            {input.gff} {input.regions} {params.feature} {params.distance} \
            {params.promoter_len} {input.vcfsnp} {threads} {params.top_regions} \
            {output.genes} {output.collapsed} > {log} 2>&1
        """

# GO enrichment analysis (only if GO_FIELD is specified)
rule run_enrichment:
    """Run GO enrichment analysis for genes around significant regions."""
    input:
        genes = O['genes_per_region_collapsed'],
        gff = f"{INDIR}{GFF}"
    output: O['enrichment']
    params:
        go_field = GO_FIELD,
        feature = GFF_FEATURE,
        tables_dir = f"{TABLES}association/enrichment/"
    log: f"{LOGDIR}run_enrichment.log"
    shell:
        """
        Rscript /pipeline/scripts/run_enrichment.R \
            {input.genes} {input.gff} {params.go_field} {params.feature} \
            {params.tables_dir} {output} > {log} 2>&1
        """

#=============================================================================
# MODULE 5: REGIONPLOT
#=============================================================================

rule gff2topr:
    """Convert GFF to topr-compatible gene annotation format."""
    input: gff = f"{INDIR}{GFF}"
    output: O['gff_topr']
    params:
        feature = GFF_FEATURE,
        genename = GFF_GENENAME,
        biotype = GFF_BIOTYPE
    log: f"{LOGDIR}gff2topr.log"
    shell:
        """
        python3 /pipeline/scripts/gff2topr.py \
            {input.gff} {params.feature} {params.genename} {params.biotype} \
            {output} > {log} 2>&1
        """

rule regionplot:
    """Generate regional Manhattan plots for top regions with all methods overlaid."""
    input:
        regions = O['regions'],
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
        plot_dir = f"{PLOTS}regionplot/",
        custom_region = REGIONPLOT_REGION,
        custom_traits = REGIONPLOT_TRAITS,
        custom_methods = REGIONPLOT_ASSOCMETHOD
    log: f"{LOGDIR}regionplot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_regionplot.R \
            {input.regions} {input.gff_topr} {params.assoc_str} \
            {params.top_regions} {params.genes} {params.plot_dir} \
            {params.custom_region} {params.custom_traits} \
            {params.custom_methods} > {log} 2>&1
        """

#=============================================================================
# MODULE 6: MALADAPTATION
#=============================================================================

# Per-model future climate download (runs in parallel via Snakemake)
rule download_climate_future_model:
    """Download CMIP6 future climate data for a single model."""
    input: samples = O['metadata']
    output: f"{INTER}Climate_future_year{YEAR}_ssp{SSP}_{{model}}.grd"
    wildcard_constraints: model = r"[A-Za-z0-9_-]+"
    params:
        crop = CROP_REGION,
        gap = GAP,
        resolution = RESOLUTION,
        data_dir = INDIR
    log: f"{LOGDIR}download_climate_future_{{model}}.log"
    shell:
        """
        Rscript /pipeline/scripts/download_climate_future_model.R \
            {input.samples} {params.crop} {params.gap} {params.data_dir} \
            {SSP} {YEAR} {wildcards.model} {params.resolution} \
            {output} > {log} 2>&1
        """

# Merge per-model rasters into averaged future climate
rule merge_climate_future:
    """Average future climate across models and extract site values."""
    input:
        samples = O['metadata'],
        model_rasters = [f"{INTER}Climate_future_year{YEAR}_ssp{SSP}_{model}.grd" for model in MODELS_LIST],
        present_raster = W['climate_raster'],
        present_all = O['climate_all']
    output:
        raster = W['climate_future_raster'],
        all_vals = O['climate_future_all'],
        site_vals = O['climate_future_site']
    params:
        raster_str = lambda wc, input: ','.join(input.model_rasters),
        n_models = len(MODELS_LIST)
    log: f"{LOGDIR}merge_climate_future.log"
    shell:
        """
        Rscript /pipeline/scripts/merge_climate_future.R \
            {input.samples} {params.raster_str} {params.n_models} \
            {input.present_raster} {input.present_all} \
            {output.raster} {output.all_vals} {output.site_vals} > {log} 2>&1
        """

rule density_plot_future:
    """Generate combined density plot for future climate predictors."""
    input: climate = O['climate_future_site']
    output: O['density_future']
    params:
        predictors = PREDICTORS_SELECTED,
        plot_dir = f"{PLOTS}{CLIMATE_TAG}/",
        inter_dir = INTER,
        prefix = f"DensityPlot_future_ssp{SSP}_{YEAR}"
    log: f"{LOGDIR}density_plot_future.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_density.R \
            {input.climate} {params.predictors} {params.plot_dir} {params.inter_dir} \
            {params.prefix} > {log} 2>&1
        """

# Gradient Forest - adaptive model
rule gradient_forest_adaptive:
    """Build adaptive Gradient Forest model using significant SNPs."""
    input:
        lfmm = W['lfmm_full'],
        sigsnps = O['selected_snps'],
        vcfsnp = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = O['metadata'],
        climate = O['climate_site']
    output: W['gf_adaptive']
    params:
        predictors = PREDICTORS_SELECTED,
        ntree = NTREE,
        cor_threshold = COR_THRESHOLD,
        pcnm = PCNM
    log: f"{LOGDIR}gradient_forest_adaptive.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_model.R \
            {input.lfmm} {input.sigsnps} {input.vcfsnp} {input.removed} \
            {input.samples} {input.climate} {params.predictors} \
            {params.ntree} {params.cor_threshold} {params.pcnm} \
            adaptive {output} > {log} 2>&1
        """

# Gradient Forest - random/neutral model (optional)
rule gradient_forest_random:
    """Build neutral Gradient Forest model using random SNPs."""
    input:
        lfmm = W['lfmm_full'],
        sigsnps = O['selected_snps'],
        vcfsnp = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = O['metadata'],
        climate = O['climate_site']
    output: W['gf_random']
    params:
        predictors = PREDICTORS_SELECTED,
        ntree = NTREE,
        cor_threshold = COR_THRESHOLD,
        pcnm = PCNM
    log: f"{LOGDIR}gradient_forest_random.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_model.R \
            {input.lfmm} {input.sigsnps} {input.vcfsnp} {input.removed} \
            {input.samples} {input.climate} {params.predictors} \
            {params.ntree} {params.cor_threshold} {params.pcnm} \
            random {output} > {log} 2>&1
        """

# Genetic offset calculation
rule gradient_forest_offset:
    """Calculate genetic offset between present and future climate."""
    input:
        gf = W['gf_adaptive'],
        future_all = O['climate_future_all'],
        present_all = O['climate_all'],
        present_raster = W['climate_raster'],
        samples = O['metadata']
    output:
        raster = W['gf_offset_raster'],
        map_values = O['gf_offset_map_values'],
        site_values = O['gf_offset_site_values']
    params:
        predictors = PREDICTORS_SELECTED
    log: f"{LOGDIR}gradient_forest_offset.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_offset.R \
            {input.gf} {params.predictors} {input.future_all} {input.present_all} \
            {input.present_raster} {input.samples} \
            {output.raster} {output.map_values} {output.site_values} > {log} 2>&1
        """

# Cumulative importance plot
rule plot_gf_cumimp:
    """Plot cumulative importance curves for adaptive (and optionally neutral) GF model."""
    input:
        gf = W['gf_adaptive'],
        gf_random = W['gf_random'] if GF_RANDOM_MODEL else []
    output: O['gf_cumimp']
    params:
        gf_random_path = W['gf_random'] if GF_RANDOM_MODEL else 'NULL',
        predictors = PREDICTORS_SELECTED,
        plot_dir = f"{PLOTS}gradientForest/",
        inter_dir = INTER,
        suffix = f"{GF_SUFFIX}_{PCNM}PCNM"
    log: f"{LOGDIR}plot_gf_cumimp.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_gf_cumimp.R \
            {input.gf} {params.gf_random_path} {params.predictors} \
            {params.plot_dir} {params.inter_dir} {params.suffix} > {log} 2>&1
        """

# Overall importance plot
rule plot_gf_importance:
    """Plot R2-weighted importance for adaptive (and optionally neutral) GF model."""
    input:
        gf = W['gf_adaptive'],
        gf_random = W['gf_random'] if GF_RANDOM_MODEL else []
    output: O['gf_importance']
    params:
        gf_random_path = W['gf_random'] if GF_RANDOM_MODEL else 'NULL',
        plot_dir = f"{PLOTS}gradientForest/",
        inter_dir = INTER,
        suffix = f"{GF_SUFFIX}_{PCNM}PCNM"
    log: f"{LOGDIR}plot_gf_importance.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_gf_importance.R \
            {input.gf} {params.gf_random_path} \
            {params.plot_dir} {params.inter_dir} {params.suffix} > {log} 2>&1
        """

# Genetic offset PieMap
rule plot_gf_offset_piemap:
    """Plot genetic offset on map with population structure pie charts (uniform pie size)."""
    input:
        offset_raster = W['gf_offset_raster'],
        samples = O['metadata'],
        clusters = clusters_table(K_BEST)
    output: O['gf_offset_piemap']
    params:
        pie_alpha = PIEMAP_PIE_ALPHA,
        pop_label = PIEMAP_POP_LABEL,
        pop_label_size = PIEMAP_POP_LABEL_SIZE,
        plot_dir = f"{PLOTS}gradientForest/",
        inter_dir = INTER,
        suffix = f"{GF_SUFFIX}_{PCNM}PCNM"
    log: f"{LOGDIR}plot_gf_offset_piemap.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.offset_raster} 1 "Genetic Offset" \
            {input.samples} {input.clusters} \
            NULL NULL \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            GeneticOffsetPieMap_{params.suffix} > {log} 2>&1
        """

rule plot_gf_offset_piemap_tajima:
    """Plot genetic offset with Tajima's D-scaled pie sizes (requires CALC_POP_STATS=TRUE)."""
    input:
        offset_raster = W['gf_offset_raster'],
        samples = O['metadata'],
        clusters = clusters_table(K_BEST),
        tajima = O['tajima']
    output: O['gf_offset_piemap_tajima']
    params:
        pie_alpha = PIEMAP_PIE_ALPHA,
        pop_label = PIEMAP_POP_LABEL,
        pop_label_size = PIEMAP_POP_LABEL_SIZE,
        plot_dir = f"{PLOTS}gradientForest/",
        inter_dir = INTER,
        suffix = f"{GF_SUFFIX}_{PCNM}PCNM"
    log: f"{LOGDIR}plot_gf_offset_piemap_tajima.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.offset_raster} 1 "Genetic Offset" \
            {input.samples} {input.clusters} \
            {input.tajima} "Tajima's D" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            GeneticOffsetPieMap_{params.suffix}_TajimaD > {log} 2>&1
        """

rule plot_gf_offset_piemap_diversity:
    """Plot genetic offset with Pi Diversity-scaled pie sizes (requires CALC_POP_STATS=TRUE)."""
    input:
        offset_raster = W['gf_offset_raster'],
        samples = O['metadata'],
        clusters = clusters_table(K_BEST),
        diversity = O['pi_div']
    output: O['gf_offset_piemap_diversity']
    params:
        pie_alpha = PIEMAP_PIE_ALPHA,
        pop_label = PIEMAP_POP_LABEL,
        pop_label_size = PIEMAP_POP_LABEL_SIZE,
        plot_dir = f"{PLOTS}gradientForest/",
        inter_dir = INTER,
        suffix = f"{GF_SUFFIX}_{PCNM}PCNM"
    log: f"{LOGDIR}plot_gf_offset_piemap_diversity.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.offset_raster} 1 "Genetic Offset" \
            {input.samples} {input.clusters} \
            {input.diversity} "Pi Diversity" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            GeneticOffsetPieMap_{params.suffix}_PiDiversity > {log} 2>&1
        """
