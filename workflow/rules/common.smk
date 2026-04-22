# ADAPTOGENE Pipeline - Refactored
# vim: filetype=python
import os
import sys
import subprocess
from pathlib import Path

# Method registries — import before config parsing so GAPIT_MODELS can be derived
sys.path.insert(0, os.path.join(workflow.basedir, "workflow"))
from methods.gea import GEA_METHODS
from methods.gwas import GWAS_METHODS
from methods.maladaptation import MALADAPTATION_METHODS

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
# Helper to read nested config with defaults
def _cfg(group, key, default=None):
    """Read config[group][key] with fallback to default."""
    return config.get(group, {}).get(key, default)

def _cfg_bool(group, key, default=False):
    """Read a boolean config value, handling YAML string representations ('false'/'true')."""
    val = _cfg(group, key, default)
    if isinstance(val, bool): return val
    if isinstance(val, str): return val.strip().lower() not in ('false', 'f', '0', 'no', '')
    return bool(val)

# INPUT parameters
INDIR = '/pipeline/' + config['Input']['dir']
PROJECT = config['project_name']
OUTDIR = f'/pipeline/{PROJECT}_results/'
LOGDIR = f'/pipeline/{PROJECT}_logs/'

CPU = config.get('cpu', max(1, os.cpu_count() - 2)); check_numeric(CPU, 'cpu')
VCF_RAW = config['Input']['vcf']; check_file_exists(INDIR, VCF_RAW, 'input.vcf')
SAMPLES = config['Input']['metadata']; check_file_exists(INDIR, SAMPLES, 'input.metadata')
VCF_BASE = get_vcf_basename(VCF_RAW)

# FILTER parameters
MAF = config['Filter']['maf']; check_float(MAF, 'filter.maf')
MISS = config['Filter']['snp_miss']; check_float(MISS, 'filter.snp_miss')
SAMPLE_MISS = _cfg('Filter', 'sample_miss', 0.5); check_float(SAMPLE_MISS, 'filter.sample_miss')
MIN_DEPTH    = _cfg('Filter', 'min_depth', None)
MAX_DEPTH    = _cfg('Filter', 'max_depth', None)
HET_OUTLIER_SD = _cfg('Filter', 'het_outlier_sd', None)

# Detect FORMAT/DP field in raw VCF header at parse time (reads header only — fast)
_raw_vcf_path = os.path.join(INDIR, VCF_RAW)
try:
    _dp_check = subprocess.run(
        ['bcftools', 'view', '-h', _raw_vcf_path],
        capture_output=True, text=True, timeout=30
    )
    HAS_FORMAT_DP = '##FORMAT=<ID=DP' in _dp_check.stdout
except Exception:
    HAS_FORMAT_DP = False

# LD parameters
LD_WIN = config['LD']['window']; check_numeric(LD_WIN, 'ld.window')
LD_STEP = config['LD']['step']; check_numeric(LD_STEP, 'ld.step')
LD_R2 = config['LD']['r2']; check_float(LD_R2, 'ld.r2')

# SNMF parameters
K_START = config['sNMF']['k_start']; check_numeric(K_START, 'snmf.k_start')
K_END = config['sNMF']['k_end']; check_numeric(K_END, 'snmf.k_end')
PLOIDY = 2  # diploid only
REPEAT = config['sNMF']['repeats']; check_numeric(REPEAT, 'snmf.repeats')
K_BEST = int(_cfg('sNMF', 'k_best', None)) if _cfg('sNMF', 'k_best', None) is not None else None
SNMF_PROJECT_MODE = 'new'  # LEA project mode: 'new' for fresh runs, 'continue' to resume

# MAP parameters
CLIMATE_EXTENT = _cfg('Map', 'climate_extent', 'auto')
GAP = _cfg('Map', 'gap', 0.5)
RESOLUTION = _cfg('Map', 'resolution', 0.5)
REGIONMAP_EXTENT = _cfg('Map', 'zoom_extent', 'NULL')
HAS_REGIONMAP = REGIONMAP_EXTENT not in ['NULL', '', None]

def get_zoom_suffix():
    """Convert regionmap coordinates to filename-safe suffix."""
    if not HAS_REGIONMAP:
        return ''
    # Replace periods with 'p' and commas with underscores
    # Example: "34.8,35.8,32.0,33.3" -> "_zoom_34p8_35p8_32p0_33p3"
    coords = str(REGIONMAP_EXTENT).replace('.', 'p').replace(',', '_')
    return f"_zoom_{coords}"

# CLIMATE parameters
CLIMATE_ENABLED = _cfg_bool('Climate', 'enabled', True)
PREDICTORS_SELECTED = _cfg('Climate', 'predictors', '') if CLIMATE_ENABLED else ''

# All 19 WorldClim bioclimatic variables — always available in the raster stack
ALL_BIO     = [f"bio_{i}" for i in range(1, 20)]
ALL_BIO_STR = ",".join(ALL_BIO)

# POP parameters
CALC_POP_STATS = _cfg_bool('Population', 'calc_stats', False)
METRICS_WINSIZE = _cfg('Population', 'window_size', 10000)
CUSTOM_TRAIT = _cfg('Population', 'custom_trait_file', 'NULL')

# PIEMAP parameters (palette fixed to viridis plasma)
PIEMAP_ALPHA = _cfg('Piemap', 'alpha', 0.6)
PIEMAP_SHOW_LABELS = 'T' if _cfg_bool('Piemap', 'show_labels', False) else 'F'
PIEMAP_LABEL_SIZE = _cfg('Piemap', 'label_size', 10)
PIEMAP_PIE_SCALE = _cfg('Piemap', 'pie_scale', 1.0)
PIEMAP_USE_POINTS = 'T' if _cfg_bool('Piemap', 'use_points', False) else 'F'

# LD DECAY parameters
_ld_decay = config.get('LDdecay', {})
LD_DECAY_GROUP_BY = _ld_decay.get('group_by', 'site')
check_in_list(LD_DECAY_GROUP_BY, ['site', 'cluster'], 'LDdecay.group_by')
LD_DECAY_MIN_SAMPLES = int(_ld_decay.get('min_samples', 10))
LD_DECAY_MAX_DISTANCE = int(_ld_decay.get('max_distance', 500))
LD_DECAY_SCOPE = _ld_decay.get('scope', 'both')
check_in_list(LD_DECAY_SCOPE, ['genome_wide', 'per_chromosome', 'both'], 'LDdecay.scope')
LD_DECAY_TAG = f"grp{LD_DECAY_GROUP_BY}_min{LD_DECAY_MIN_SAMPLES}_dist{LD_DECAY_MAX_DISTANCE}"

# GFF parameters
GFF = _cfg('Input', 'gff', '')
GFF_FEATURE = _cfg('GFF', 'feature', 'mRNA')
GFF_GENE_NAME = _cfg('GFF', 'gene_name', 'description')
GFF_BIOTYPE = _cfg('GFF', 'biotype', 'biotype')

# ASSOCIATION parameters
_assoc = config.get('GEA', {})
_gff = config.get('GFF', {})
SIGSNPS_METHOD = 'All'  # pipeline always uses All; combine strategy moved to gradient_forest config

def resolve_region_params(group_dict, defaults=None):
    """Parse region_distance, combine_gap, promoter_length from a config group dict."""
    defaults = defaults or {}
    raw_rdist = group_dict.get('region_distance', defaults.get('region_distance', 1000000))
    auto = (str(raw_rdist).lower() == 'auto')
    return {
        'region_distance':      'auto' if auto else int(raw_rdist),
        'region_distance_auto': auto,
        'combine_gap':          int(group_dict.get('combine_gap',      defaults.get('combine_gap', 100000))),
        'promoter_length':      int(group_dict.get('promoter_length',  defaults.get('promoter_length', 10000))),
    }

_gea_rdp = resolve_region_params(_assoc)
PROMOTER_LENGTH    = _gea_rdp['promoter_length']
SIGSNPS_GAP        = _gea_rdp['combine_gap']
SNP_DISTANCE       = SIGSNPS_GAP
REGION_DISTANCE    = _gea_rdp['region_distance']
REGION_DISTANCE_AUTO = _gea_rdp['region_distance_auto']
GO_FIELD           = _gff.get('go_field', 'NULL')

# ENRICHMENT parameters
_enrich = config.get('Enrichment', {})
ENRICHMENT_TOP_TERMS = _enrich.get('top_terms', 20)
ENRICHMENT_PLOT_WIDTH = _enrich.get('plot_width', 12)
ENRICHMENT_PLOT_HEIGHT = _enrich.get('plot_height', 10)
ENRICHMENT_CNET_LABEL = _enrich.get('cnet_label', 'gene_id')

# FUTURE parameters
_future = config.get('Future', {})
SSP = _future.get('ssp', '585')
YEAR = _future.get('year', '2061-2080')
MODELS_STR = _future.get('models', '')
MODELS_LIST = [m.strip() for m in MODELS_STR.split(',') if m.strip()] if MODELS_STR else []

# GRADIENT FOREST parameters
_gf = config.get('Maladaptation', {}).get('methods', {}).get('gradient_forest', {})
NTREE = _gf.get('ntree', '500')
COR_THRESHOLD = _gf.get('cor_threshold', '0.5')
SPATIAL_CORRECTION = _gf.get('spatial_correction', 'with')
GF_RUN_LABEL = _gf.get('run_label', '')
GF_RANDOM_MODEL = _gf.get('random_model', True)
GF_COMBINE_METHOD = _gf.get('combine_method', 'All')
GF_COMBINE_GAP = int(_gf.get('combine_gap', 100000))
SPATIAL_TAG = 'spatial' if SPATIAL_CORRECTION == 'with' else 'nospatial'

# Validate and derive active maladaptation methods from config
_mala_methods_cfg = config.get('Maladaptation', {}).get('methods', {})
_unknown_mala = set(_mala_methods_cfg.keys()) - set(MALADAPTATION_METHODS.keys())
if _unknown_mala:
    raise ValueError(f"Unknown maladaptation methods in config: {_unknown_mala}. "
                     f"Registered: {list(MALADAPTATION_METHODS.keys())}")
ACTIVE_MALA_METHODS = [m for m in _mala_methods_cfg if m in MALADAPTATION_METHODS] or ['gradient_forest']

def parse_association_configs(configs_list):
    """Parse association configs list into method -> adjust_threshold dict."""
    configs = {}
    for cfg in (configs_list or []):
        method = cfg['method']
        adjust = f"{cfg['adjust']}_{cfg['threshold']}"
        if method in configs:
            raise ValueError(f"Method '{method}' appears multiple times in configs")
        configs[method] = adjust
    return configs

GEA_CONFIGS = parse_association_configs(_assoc.get('configs', []))

# GAPIT model detection — derived from registry (single source of truth)
GAPIT_MODELS = {name for name, cfg in GEA_METHODS.items() if cfg["engine"] == "gapit"}

def split_configs_by_engine(configs):
    """Split configs into GAPIT and non-GAPIT (legacy) engines."""
    gapit = {m: a for m, a in configs.items() if m in GAPIT_MODELS}
    other = {m: a for m, a in configs.items() if m not in GAPIT_MODELS}
    return gapit, other

GEA_GAPIT_CONFIGS, GEA_OTHER_CONFIGS = split_configs_by_engine(GEA_CONFIGS)

# Dynamic wildcard constraint regex for association methods
GEA_METHOD_REGEX = '|'.join(GEA_CONFIGS.keys()) if GEA_CONFIGS else 'EMMAX'

# PHENOTYPE ASSOCIATION parameters (inherits from GEA.* by default)
_pheno = config.get('GWAS', {})
GWAS_CONFIGS = parse_association_configs(_pheno.get('configs', []))
GWAS_GAPIT_CONFIGS, GWAS_OTHER_CONFIGS = split_configs_by_engine(GWAS_CONFIGS)
GWAS_METHOD_REGEX = '|'.join(GWAS_CONFIGS.keys()) if GWAS_CONFIGS else 'EMMAX'
PHENO_MISSING = _pheno.get('missing_strategy', 'DROP')
if GWAS_CONFIGS:
    check_in_list(PHENO_MISSING, ['MEAN', 'MEDIAN', 'DROP'], 'GWAS.missing_strategy')

def _validate_methods_in_registry(configs, registry, context):
    """Fail early with a clear error when a config method is not in the registry."""
    unknown = [m for m in configs if m not in registry]
    if unknown:
        raise ValueError(
            f"Unknown method(s) in {context}: {unknown}. "
            f"Known methods: {sorted(registry.keys())}"
        )

if GEA_CONFIGS:
    _validate_methods_in_registry(GEA_CONFIGS, GEA_METHODS, "GEA.configs")
if GWAS_CONFIGS:
    _validate_methods_in_registry(GWAS_CONFIGS, GWAS_METHODS, "GWAS.configs")

# Inherit from GEA.* with optional override
PHENO_COMBINE_METHOD = 'All'  # pipeline always uses All; combine strategy moved to gradient_forest config
_gwas_rdp = resolve_region_params(_pheno, defaults=_gea_rdp)
PHENO_COMBINE_GAP          = _gwas_rdp['combine_gap']
PHENO_SNP_DISTANCE         = PHENO_COMBINE_GAP
PHENO_REGION_DISTANCE      = _gwas_rdp['region_distance']
PHENO_REGION_DISTANCE_AUTO = _gwas_rdp['region_distance_auto']
PHENO_PROMOTER_LENGTH      = _gwas_rdp['promoter_length']

# BLOCK MODE parameters (LD block partitioning + WZA p-value aggregation)
GEA_BLOCK_MODE = _assoc.get('block_mode', 'snp')
check_in_list(GEA_BLOCK_MODE, ['snp', 'block'], 'GEA.block_mode')
GWAS_BLOCK_MODE = _pheno.get('block_mode', 'snp') if config.get('GWAS') else 'snp'
if config.get('GWAS'):
    check_in_list(GWAS_BLOCK_MODE, ['snp', 'block'], 'GWAS.block_mode')
_ANY_BLOCK_MODE = GEA_BLOCK_MODE == 'block' or GWAS_BLOCK_MODE == 'block'

# Shared block parameters — prefer GEA.block when GEA is in block mode, else GWAS.block
_block_cfg = _assoc.get('block', {}) if GEA_BLOCK_MODE == 'block' else _pheno.get('block', {})
BLOCK_BLOCKS_MAX_KB     = int(_block_cfg.get('blocks_max_kb', 200))
_block_wza              = _block_cfg.get('wza', {})
BLOCK_WZA_MAF_FILTER    = float(_block_wza.get('maf_filter', 0.05))
BLOCK_WZA_MIN_SNPS      = int(_block_wza.get('min_snps_per_block', 3))
_block_within           = _block_cfg.get('within_block', {})
BLOCK_WITHIN_METHOD     = _block_within.get('method', 'within_order')
check_in_list(BLOCK_WITHIN_METHOD,
    ['threshold', 'top_n', 'lead', 'within_order', 'all'], 'block.within_block.method')
BLOCK_WITHIN_MULTIPLIER = float(_block_within.get('multiplier', 10.0))
BLOCK_WITHIN_THRESHOLD  = float(_block_within.get('threshold', 0.05))
BLOCK_WITHIN_TOP_N      = int(_block_within.get('top_n', 5))
BLOCK_ADJUST_STR        = str(_block_cfg.get('adjust', 'bonf_0.05'))

# OVERLAP parameters (GEA + GWAS combined analysis)
_overlap = config.get('GEAxGWAS', {})
_overlap_rdist_raw = _overlap.get('region_distance', None)
if _overlap_rdist_raw is None:
    # Default: inherit auto if either source is auto, otherwise use max of both
    if REGION_DISTANCE_AUTO or PHENO_REGION_DISTANCE_AUTO:
        OVERLAP_REGION_DISTANCE_AUTO = True
        OVERLAP_REGION_DISTANCE = 'auto'
    else:
        OVERLAP_REGION_DISTANCE_AUTO = False
        OVERLAP_REGION_DISTANCE = max(REGION_DISTANCE, PHENO_REGION_DISTANCE)
else:
    OVERLAP_REGION_DISTANCE_AUTO = (str(_overlap_rdist_raw).lower() == 'auto')
    OVERLAP_REGION_DISTANCE = 'auto' if OVERLAP_REGION_DISTANCE_AUTO else int(_overlap_rdist_raw)
# PAIRWISE OVERLAP parameters (trait-vs-trait comparison across all sources)
_pairwise = _overlap.get('pairwise', {})
PAIRWISE_WINDOW_SIZE = int(_pairwise.get('window_size', 500000))
PAIRWISE_MIN_SNPS    = int(_pairwise.get('min_snps', 2))

# Discover phenotype traits from metadata columns 5+ (for directory creation and wildcards)
PHENO_TRAITS = []
PHENO_PREDICTORS = ''
_meta_path = os.path.join(INDIR, SAMPLES)
try:
    with open(_meta_path) as _f:
        _meta_header = _f.readline().strip().split('\t')
    META_HAS_PHENO = len(_meta_header) > 4
except Exception:
    META_HAS_PHENO = False
if GWAS_CONFIGS:
    PHENO_TRAITS = _meta_header[4:] if META_HAS_PHENO else []
    # Filter by GWAS.traits if specified in config
    _pheno_trait_filter = _pheno.get('traits', '')
    if _pheno_trait_filter:
        _allowed = [t.strip() for t in str(_pheno_trait_filter).split(',') if t.strip()]
        if _allowed:
            PHENO_TRAITS = [t for t in PHENO_TRAITS if t in _allowed]
    PHENO_PREDICTORS = ','.join(PHENO_TRAITS)

#=============================================================================
# PATH DEFINITIONS
#=============================================================================
# Directory tags (easy to modify if adding new parameters)
_depth_tag = ""
if MIN_DEPTH is not None or MAX_DEPTH is not None:
    _dp_parts = []
    if MIN_DEPTH is not None: _dp_parts.append(f"dpmin{MIN_DEPTH}")
    if MAX_DEPTH is not None: _dp_parts.append(f"dpmax{MAX_DEPTH}")
    _depth_tag = "_" + "_".join(_dp_parts)
_het_tag = f"_hetsd{HET_OUTLIER_SD}" if HET_OUTLIER_SD is not None else ""
FILT_TAG = f"maf{MAF}_miss{MISS}_smiss{SAMPLE_MISS}{_depth_tag}{_het_tag}"
LD_TAG = f"ld{LD_R2}_win{LD_WIN}_step{LD_STEP}"

# Base directories (new module-based layout)
WORK = f"{OUTDIR}_work/"
INTER = f"{OUTDIR}_intermediate/"

# Working subdirectories (parameter-based structure)
WORK_FILT = f"{WORK}{FILT_TAG}/"
WORK_LD = f"{WORK_FILT}{LD_TAG}/"

# Module output directories
MOD_PROCESSING = f"{OUTDIR}Processing/"
MOD_PRESTRUCT  = f"{OUTDIR}PreStructure/"
MOD_CLIMATE    = f"{OUTDIR}climate/"
MOD_STRUCT     = f"{OUTDIR}Structure/"
MOD_GEA        = f"{OUTDIR}GEA/"
MOD_GWAS       = f"{OUTDIR}GWAS/"
MOD_GEAXGWAS   = f"{OUTDIR}GEAxGWAS/"
MOD_MALAD      = f"{OUTDIR}Maladaptation/"

# Working paths
W = {
    # Samples (in intermediate - shared across all filtering)
    'samples_list': f"{INTER}samples/samples.list",
    'samples_filtered': f"{INTER}samples/samples_filtered.list",
    'samples_removed': f"{INTER}samples/samples_removed.list",
    'samples_missing_stats': f"{INTER}samples/samples_missing_stats.tsv",
    'samples_order': f"{INTER}samples/samples_order.list",
    # Optional: samples list after het outlier removal (used by filter_vcf when HET_OUTLIER_SD set)
    'samples_het_filtered': f"{INTER}samples/samples_het_filtered.list",
    # Normalized GFF (chromosome names match VCF - 'chr' prefix stripped)
    'gff_normalized': f"{INTER}annotation/normalized.gff3",
    # Optional: depth-masked VCF (before MAF/missingness filter, when min/max_depth set)
    'vcf_depth_filtered': f"{WORK}{VCF_BASE}_depth_filtered.vcf",
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
    'summary_done': f"{INTER}flags/summary_done.flag",
    'removed': f"{WORK_LD}{VCF_BASE}.removed",
    'snmf': f"{WORK_LD}{VCF_BASE}.snmfProject",
    # QC intermediate stats (used by plot_qc_processing)
    'het_stats_raw': f"{INTER}qc/het_raw.het",
    'het_stats_filtered': f"{INTER}qc/het_filtered.het",
}

# Output paths (organized by module)
O = {
    'summary': f"{OUTDIR}pipeline_summary.tsv",
    'metadata': f"{MOD_PROCESSING}tables/metadata.tsv",
    'sample_missing_stats': f"{MOD_PROCESSING}tables/sample_missing_stats.tsv",
    'pca': f"{MOD_PRESTRUCT}plots/pca.png",
    'pca_svg': f"{MOD_PRESTRUCT}plots/pca.svg",
    'tracy': f"{MOD_PRESTRUCT}plots/tracy_widom.png",
    'cross_entropy': f"{MOD_PRESTRUCT}plots/cross_entropy_K{K_START}-{K_END}.png",
    # QC intermediate stats
    'qc_raw_summary':      f"{MOD_PROCESSING}tables/vcf_raw_summary.tsv",
    'qc_filtering_summary':f"{MOD_PROCESSING}tables/filtering_summary.tsv",
    'qc_sample_het':       f"{MOD_PROCESSING}tables/sample_heterozygosity.tsv",
    'qc_maf_raw':          f"{INTER}qc/maf_raw.frq",
    'qc_maf_filtered':     f"{INTER}qc/maf_filtered.frq",
    'qc_maf_pos':          f"{INTER}qc/maf_pos.tsv",
    'qc_snp_miss_raw':     f"{INTER}qc/snp_missingness_raw.lmiss",
    'qc_snp_density_raw':  f"{INTER}qc/snp_density_raw.snpden",
    'qc_snp_density_filt': f"{INTER}qc/snp_density_filtered.snpden",
    'qc_depth_sample':     f"{INTER}qc/depth_per_sample.idepth",
    'qc_depth_site':       f"{INTER}qc/depth_per_site.ldepth.mean",
    # QC plots
    'qc_plot_sample_miss': f"{MOD_PROCESSING}plots/sample_missingness_distribution.png",
    'qc_plot_het_miss':    f"{MOD_PROCESSING}plots/het_vs_missingness.png",
    'qc_plot_maf':         f"{MOD_PROCESSING}plots/maf_distribution.png",
    'qc_plot_snp_miss':    f"{MOD_PROCESSING}plots/snp_missingness_distribution.png",
    'qc_plot_attrition':   f"{MOD_PROCESSING}plots/filtering_attrition.png",
    'qc_plot_snp_density': f"{MOD_PROCESSING}plots/snp_density_by_chr.png",
    'qc_plot_depth':       f"{MOD_PROCESSING}plots/depth_distribution.png",
}

# K_BEST dependent paths — pre-populated with a sentinel so .smk files can be parsed
# unconditionally even when K_BEST is None (e.g. processing mode with no k_best set).
# Snakemake rejects empty strings but accepts any non-empty path at parse time.
# add_kbest_paths() / add_association_paths() / add_maladaptation_paths() overwrite
# these with real values when K_BEST is set; rules using them only enter the DAG
# when K_BEST is available, so the sentinel is never actually used as a file path.
def _ph(key):
    """Return a unique non-empty placeholder path for a K_BEST-dependent dict key.
    Snakemake rejects empty strings but accepts any non-empty path at parse time.
    These rules only enter the DAG when K_BEST is set, so the placeholder is never
    used as an actual file path."""
    return f"__placeholder_{key}__"
# --- from add_kbest_paths() ---
W['lfmm_imp']             = _ph('lfmm_imp')
W['vcf_imp']              = _ph('vcf_imp')
W['climate_raster']       = _ph('climate_raster')
W['ld_decay_work']        = _ph('ld_decay_work')
W['ld_decay_sample_lists']= _ph('ld_decay_sample_lists')
W['ld_decay_chr_vcfs']    = _ph('ld_decay_chr_vcfs')
W['ld_decay_stat_gw']     = _ph('ld_decay_stat_gw')
W['ld_decay_stat_chr']    = _ph('ld_decay_stat_chr')
W['ld_decay_prep_done']   = _ph('ld_decay_prep_done')
W['ld_decay_gw_done']     = _ph('ld_decay_gw_done')
W['ld_decay_chr_done']    = _ph('ld_decay_chr_done')
O['climate_site']         = _ph('climate_site')
O['climate_site_scaled']  = _ph('climate_site_scaled')
O['climate_all']          = _ph('climate_all')
O['climate_invariant']    = _ph('climate_invariant')
O['tajima']               = _ph('tajima')
O['pi_div']               = _ph('pi_div')
O['ibd_raw']              = _ph('ibd_raw')
O['ibd_pairs']            = _ph('ibd_pairs')
O['amova']                = _ph('amova')
O['corr_heatmap']         = _ph('corr_heatmap')
O['mantel']               = _ph('mantel')
O['amova_plot']           = _ph('amova_plot')
O['ld_decay_table']       = _ph('ld_decay_table')
O['ld_decay_plot_gw']     = _ph('ld_decay_plot_gw')
O['ld_decay_plot_gw_svg'] = _ph('ld_decay_plot_gw_svg')
O['ld_decay_plot_chr']    = _ph('ld_decay_plot_chr')
O['ld_decay_plot_chr_svg']= _ph('ld_decay_plot_chr_svg')
# --- from add_association_paths() ---
W['blocks_det']           = _ph('blocks_det')
W['geno_full']            = _ph('geno_full')
W['lfmm_full']            = _ph('lfmm_full')
W['vcfsnp_full']          = _ph('vcfsnp_full')
W['removed_full']         = _ph('removed_full')
W['snmf_full']            = _ph('snmf_full')
W['lfmm_imp_full']        = _ph('lfmm_imp_full')
W['vcf_imp_full']         = _ph('vcf_imp_full')
W['emmax_work']           = _ph('emmax_work')
W['assoc_tped']           = _ph('assoc_tped')
W['assoc_tfam']           = _ph('assoc_tfam')
W['assoc_kinship']        = _ph('assoc_kinship')
W['gapit_gd']             = _ph('gapit_gd')
W['gapit_gm']             = _ph('gapit_gm')
W['gapit_work']           = _ph('gapit_work')
O['selected_snps']        = _ph('selected_snps')
O['regions_per_trait']    = _ph('regions_per_trait')
O['regions_combined']     = _ph('regions_combined')
O['genes_per_region']     = _ph('genes_per_region')
O['genes_per_region_collapsed'] = _ph('genes_per_region_collapsed')
O['genes_combined_regions']     = _ph('genes_combined_regions')
O['manhattan_combined']   = _ph('manhattan_combined')
O['qq_combined']          = _ph('qq_combined')
# --- from add_maladaptation_paths() ---
W['climate_future_raster']= _ph('climate_future_raster')
O['climate_future_all']   = _ph('climate_future_all')
O['climate_future_site']  = _ph('climate_future_site')
O['density_future']       = _ph('density_future')

def add_kbest_paths():
    """Add K_BEST dependent paths to W and O dictionaries."""
    if K_BEST is None:
        return

    # Imputed files in LD directory
    W['lfmm_imp'] = f"{WORK_LD}{VCF_BASE}_K{K_BEST}imp.lfmm"
    W['vcf_imp'] = f"{WORK_LD}{VCF_BASE}_K{K_BEST}imp.vcf"

    # Climate rasters
    W['climate_raster'] = f"{MOD_CLIMATE}rasters/present/climate_present_rasterstack.tif"

    # Tables - climate
    O['climate_site'] = f"{MOD_CLIMATE}tables/present/climate_present_site.tsv"
    O['climate_site_scaled'] = f"{MOD_CLIMATE}tables/present/climate_present_site_scaled.tsv"
    O['climate_all'] = f"{MOD_CLIMATE}tables/present/climate_present_all.tsv"
    O['climate_invariant'] = f"{MOD_CLIMATE}tables/present/climate_invariant_predictors.tsv"
    # Tables - structure_k/population stats
    O['tajima'] = f"{MOD_STRUCT}tables/pop_stats/tajima_d_by_pop.tsv"
    O['pi_div'] = f"{MOD_STRUCT}tables/pop_stats/pi_diversity_by_pop.tsv"
    O['ibd_raw'] = f"{MOD_STRUCT}tables/pop_stats/ibd_raw.tsv"
    O['ibd_pairs'] = f"{MOD_STRUCT}tables/pop_stats/ibd_pairs.tsv"
    O['amova'] = f"{MOD_STRUCT}tables/pop_stats/amova.tsv"

    # Plots - climate
    O['corr_heatmap'] = f"{MOD_CLIMATE}plots/correlation_heatmap.png"
    # Plots - structure_k
    O['mantel'] = f"{MOD_STRUCT}plots/pop_stats/mantel_test.png"
    O['amova_plot'] = f"{MOD_STRUCT}plots/pop_stats/amova.png"

    # LD decay paths (parameterized by group_by, min_samples, max_distance)
    _ld_base = f"{INTER}ld_decay/{LD_DECAY_TAG}/"
    W['ld_decay_work']         = _ld_base
    W['ld_decay_sample_lists'] = f"{_ld_base}sample_lists/"
    W['ld_decay_chr_vcfs']     = f"{_ld_base}chr_vcfs/"
    W['ld_decay_stat_gw']      = f"{_ld_base}stat_gw/"
    W['ld_decay_stat_chr']     = f"{_ld_base}stat_chr/"
    W['ld_decay_prep_done']    = f"{_ld_base}prep_done.flag"
    W['ld_decay_gw_done']      = f"{_ld_base}gw_done.flag"
    W['ld_decay_chr_done']     = f"{_ld_base}chr_done.flag"
    O['ld_decay_table'] = f"{MOD_STRUCT}tables/ld_decay_half_distances.tsv"
    O['ld_decay_plot_gw'] = f"{MOD_STRUCT}plots/ld_decay/ld_decay_genome_wide.png"
    O['ld_decay_plot_gw_svg'] = f"{MOD_STRUCT}plots/ld_decay/ld_decay_genome_wide.svg"
    O['ld_decay_plot_chr'] = f"{MOD_STRUCT}plots/ld_decay/ld_decay_per_chromosome.png"
    O['ld_decay_plot_chr_svg'] = f"{MOD_STRUCT}plots/ld_decay/ld_decay_per_chromosome.svg"

add_kbest_paths()

# Association paths (added when K_BEST is set and association mode is used)
def add_association_paths():
    """Add association-specific paths to W and O dictionaries."""
    if K_BEST is None or not GEA_CONFIGS:
        return

    # Full (non-LD pruned) files for association analysis
    W['geno_full'] = f"{WORK_FILT}{VCF_BASE}.geno"
    W['lfmm_full'] = f"{WORK_FILT}{VCF_BASE}.lfmm"
    W['vcfsnp_full'] = f"{WORK_FILT}{VCF_BASE}.vcfsnp"
    W['removed_full'] = f"{WORK_FILT}{VCF_BASE}.removed"
    W['snmf_full'] = f"{WORK_FILT}{VCF_BASE}.snmfProject"
    W['lfmm_imp_full'] = f"{WORK_FILT}{VCF_BASE}_K{K_BEST}imp.lfmm"
    W['vcf_imp_full'] = f"{WORK_FILT}{VCF_BASE}_K{K_BEST}imp.vcf"

    # Engine-specific working paths — only populate for engines actually in use
    _assoc_engines = {GEA_METHODS[m]["engine"] for m in GEA_CONFIGS}
    if "emmax" in _assoc_engines:
        W['emmax_work']    = f"{WORK_FILT}emmax/"
        W['assoc_tped']    = f"{WORK_FILT}emmax/{VCF_BASE}.tped"
        W['assoc_tfam']    = f"{WORK_FILT}emmax/{VCF_BASE}.tfam"
        W['assoc_kinship'] = f"{WORK_FILT}emmax/{VCF_BASE}.aBN.kinf"
    if "gapit" in _assoc_engines:
        W['gapit_gd']   = f"{WORK_FILT}gapit/{VCF_BASE}_GD.tsv"
        W['gapit_gm']   = f"{WORK_FILT}gapit/{VCF_BASE}_GM.tsv"
        W['gapit_work'] = f"{INTER}gapit/gea/"
    # lfmm uses W['lfmm_imp'] / W['lfmm_imp_full'] from add_kbest_paths()

    # LD blocks paths (shared across GEA/GWAS when any source uses block_mode=block)
    if _ANY_BLOCK_MODE:
        W['blocks_det'] = f"{INTER}ld_blocks/blocks.det"
        os.makedirs(f"{INTER}ld_blocks/", exist_ok=True)

    # Combined outputs - association
    O['selected_snps'] = f"{MOD_GEA}tables/selected_snps.tsv"
    O['regions_per_trait'] = f"{MOD_GEA}tables/regions_per_trait.tsv"
    O['regions_combined'] = f"{MOD_GEA}tables/regions_combined.tsv"
    O['genes_per_region'] = f"{MOD_GEA}tables/genes_per_region.tsv"
    O['genes_per_region_collapsed'] = f"{MOD_GEA}tables/genes_per_region_collapsed.tsv"
    O['genes_combined_regions'] = f"{MOD_GEA}tables/genes_combined.tsv"
    # Combined Manhattan plots (all traits and methods)
    O['manhattan_combined'] = f"{MOD_GEA}plots/manhattan/combined/manhattan_combined_K{K_BEST}.png"
    O['qq_combined'] = f"{MOD_GEA}plots/manhattan/combined/qq_combined_K{K_BEST}.png"

add_association_paths()

# Phenotype association paths
def add_pheno_association_paths():
    """Add phenotype association paths to W and O dictionaries."""
    if K_BEST is None or not GWAS_CONFIGS:
        return

    PHENO_WORK = f"{WORK_FILT}phenotypes/"

    # Working paths
    W['pheno_work'] = PHENO_WORK
    W['pheno_emmax_work'] = f"{PHENO_WORK}emmax/"
    W['pheno_gapit_work'] = f"{INTER}gapit/gwas/"

    # Ensure full-VCF derived paths exist (needed by unified downstream gene-finding rules)
    if 'vcfsnp_full' not in W:
        W['geno_full']    = f"{WORK_FILT}{VCF_BASE}.geno"
        W['lfmm_full']    = f"{WORK_FILT}{VCF_BASE}.lfmm"
        W['vcfsnp_full']  = f"{WORK_FILT}{VCF_BASE}.vcfsnp"
        W['removed_full'] = f"{WORK_FILT}{VCF_BASE}.removed"

    # Engine-specific paths for phenotype methods — only populate if not already set by GEA
    _pheno_engines = {GWAS_METHODS[m]["engine"] for m in GWAS_CONFIGS}
    if "gapit" in _pheno_engines and 'gapit_gd' not in W:
        W['gapit_gd']   = f"{WORK_FILT}gapit/{VCF_BASE}_GD.tsv"
        W['gapit_gm']   = f"{WORK_FILT}gapit/{VCF_BASE}_GM.tsv"
        W['gapit_work'] = f"{INTER}gapit/gea/"
    if "emmax" in _pheno_engines and 'assoc_kinship' not in W:
        W['emmax_work']    = f"{WORK_FILT}emmax/"
        W['assoc_tped']    = f"{WORK_FILT}emmax/{VCF_BASE}.tped"
        W['assoc_tfam']    = f"{WORK_FILT}emmax/{VCF_BASE}.tfam"
        W['assoc_kinship'] = f"{WORK_FILT}emmax/{VCF_BASE}.aBN.kinf"

    # Path A (MEAN/MEDIAN) — single set of TPED/kinship
    W['pheno_tped'] = f"{PHENO_WORK}emmax/{VCF_BASE}.tped"
    W['pheno_tfam'] = f"{PHENO_WORK}emmax/{VCF_BASE}.tfam"
    W['pheno_kinship'] = f"{PHENO_WORK}emmax/{VCF_BASE}.aBN.kinf"
    W['pheno_all_phenotypes'] = f"{PHENO_WORK}all_phenotypes.tsv"

    # Prep outputs
    O['pheno_missing_summary'] = f"{MOD_GWAS}tables/phenotype_missing_summary.tsv"

    # Combined analysis outputs
    O['pheno_selected_snps'] = f"{MOD_GWAS}tables/selected_snps.tsv"
    O['pheno_regions_per_trait'] = f"{MOD_GWAS}tables/regions_per_trait.tsv"
    O['pheno_regions_combined'] = f"{MOD_GWAS}tables/regions_combined.tsv"
    O['pheno_genes_per_region'] = f"{MOD_GWAS}tables/genes_per_region.tsv"
    O['pheno_genes_collapsed'] = f"{MOD_GWAS}tables/genes_per_region_collapsed.tsv"
    O['pheno_genes_combined'] = f"{MOD_GWAS}tables/genes_combined.tsv"
    # Manhattan combined
    O['pheno_manhattan_combined'] = f"{MOD_GWAS}plots/manhattan/combined/manhattan_combined_K{K_BEST}.png"
    O['pheno_qq_combined'] = f"{MOD_GWAS}plots/manhattan/combined/qq_combined_K{K_BEST}.png"

add_pheno_association_paths()

# Overlap paths (GEA + GWAS combined analysis)
def add_overlap_paths():
    """Add overlap analysis paths to W and O dictionaries."""
    if K_BEST is None:
        return

    has_gea  = bool(GEA_CONFIGS)
    has_gwas = bool(GWAS_CONFIGS)

    # Pairwise trait overlap — available when at least one source exists
    if has_gea or has_gwas:
        O['pairwise_collapsed_snps'] = f"{MOD_GEAXGWAS}tables/pairwise_collapsed_snps.tsv"
        O['pairwise_overlap_table']  = f"{MOD_GEAXGWAS}tables/pairwise_overlap_table.tsv"

    # GEA × GWAS combined analysis — requires both sources
    if has_gea and has_gwas:
        # Miami plots (static GEA vs GWAS — background PNG + coords JSON for Shiny overlay)
        O['overlap_miami'] = f"{MOD_GEAXGWAS}plots/miami_combined_K{K_BEST}.png"
        O['overlap_miami_svg'] = f"{MOD_GEAXGWAS}plots/miami_combined_K{K_BEST}.svg"

add_overlap_paths()

# Templates for maladaptation outputs (method-parameterized)
def _mala_suffix(run_label, spatial_tag):
    return f"{run_label}_{spatial_tag}"

def mala_plot_dir(method, run_label, spatial_tag):
    return f"{MOD_MALAD}plots/{method}/{_mala_suffix(run_label, spatial_tag)}/"

def mala_table_dir(method, run_label, spatial_tag):
    return f"{MOD_MALAD}tables/{method}/{_mala_suffix(run_label, spatial_tag)}/"

def mala_inter_dir(method, run_label, spatial_tag):
    return f"{INTER}{method}/{_mala_suffix(run_label, spatial_tag)}/"

def mala_model(method, run_label, spatial_tag, kind):
    """kind: 'adaptive' or 'random'"""
    return f"{mala_inter_dir(method, run_label, spatial_tag)}{kind}_model.qs"

def mala_selected_snps(method, run_label, spatial_tag):
    return f"{mala_inter_dir(method, run_label, spatial_tag)}selected_snps.tsv"

def mala_offset_raster(method, run_label, spatial_tag):
    return f"{mala_inter_dir(method, run_label, spatial_tag)}offset_raster.tif"

def mala_offset_map_values(method, run_label, spatial_tag):
    return f"{mala_table_dir(method, run_label, spatial_tag)}genetic_offset_map.tsv"

def mala_offset_site_values(method, run_label, spatial_tag):
    return f"{mala_table_dir(method, run_label, spatial_tag)}genetic_offset_site.tsv"

def mala_cumimp(method, run_label, spatial_tag):
    return f"{mala_plot_dir(method, run_label, spatial_tag)}cumulative_importance.png"

def mala_importance(method, run_label, spatial_tag):
    return f"{mala_plot_dir(method, run_label, spatial_tag)}overall_importance.png"

def mala_offset_piemap(method, run_label, spatial_tag, size_trait):
    """size_trait: 'notrait' | 'tajima_d' | 'pi_diversity'"""
    return f"{mala_plot_dir(method, run_label, spatial_tag)}genetic_offset_piemap_{size_trait}.png"

# Maladaptation paths
def add_maladaptation_paths():
    """Add maladaptation-specific paths to W and O dictionaries."""
    if K_BEST is None or not MODELS_LIST:
        return

    # Per-model future climate rasters (method-agnostic)
    for model in MODELS_LIST:
        W[f'climate_future_{model}'] = f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_{model}.tif"

    # Merged future climate (method-agnostic)
    W['climate_future_raster'] = f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_rasterstack.tif"
    O['climate_future_all'] = f"{MOD_CLIMATE}tables/future/climate_future_year{YEAR}_ssp{SSP}_all.tsv"
    O['climate_future_site'] = f"{MOD_CLIMATE}tables/future/climate_future_year{YEAR}_ssp{SSP}_site.tsv"

    # Future climate density plot (method-agnostic)
    O['density_future'] = f"{MOD_CLIMATE}plots/density_plot_future_ssp{SSP}_{YEAR}.png"

    # Create method-level directories
    for method in ACTIVE_MALA_METHODS:
        os.makedirs(mala_plot_dir(method, GF_RUN_LABEL, SPATIAL_TAG), exist_ok=True)
        os.makedirs(mala_table_dir(method, GF_RUN_LABEL, SPATIAL_TAG), exist_ok=True)
        os.makedirs(mala_inter_dir(method, GF_RUN_LABEL, SPATIAL_TAG), exist_ok=True)

add_maladaptation_paths()

# Templates for K-dependent outputs
def clusters_table(k): return f"{MOD_PRESTRUCT}tables/K{k}/clusters_K{k}.tsv"
def structure_plot(k): return f"{MOD_PRESTRUCT}plots/K{k}/structure_K{k}.png"
def pca_struct_plot(k): return f"{MOD_PRESTRUCT}plots/K{k}/pca_structure_K{k}.png"
def pop_diff_plot(k): return f"{MOD_PRESTRUCT}plots/K{k}/pop_diff_K{k}.png"

# Templates for climate/trait-dependent outputs
DENSITY_PLOT_COMBINED    = f"{MOD_CLIMATE}plots/density_plot_present.png"
DENSITY_PLOT_PHENOTYPES  = f"{MOD_CLIMATE}plots/density_plot_phenotypes.png"
def piemap_tajima(bio): return f"{MOD_STRUCT}plots/piemap/piemap_{bio}_tajima_d.png"
def piemap_diversity(bio): return f"{MOD_STRUCT}plots/piemap/piemap_{bio}_pi_diversity.png"
def piemap_notrait(bio): return f"{MOD_STRUCT}plots/piemap/piemap_{bio}.png"

# Templates for association outputs
def assoc_pvalues(method): return f"{MOD_GEA}tables/methods/{method}/{method}_pvalues_K{K_BEST}.tsv"
def assoc_sigsnps(method, adjust): return f"{MOD_GEA}tables/methods/{method}/{method}_pvalues_K{K_BEST}_sig_snps_{adjust}.tsv"
def manhattan_plot(method, trait, adjust): return f"{MOD_GEA}plots/manhattan/{method}/manhattan_{trait}_K{K_BEST}_{adjust}.png"
def qq_plot(method, trait, adjust): return f"{MOD_GEA}plots/manhattan/{method}/qq_{trait}_K{K_BEST}_{adjust}.png"

# Templates for phenotype association outputs
def pheno_pvalues(method): return f"{MOD_GWAS}tables/methods/{method}/{method}_pvalues_K{K_BEST}.tsv"
def pheno_qvalues(method): return f"{MOD_GWAS}tables/methods/{method}/{method}_qvalues_K{K_BEST}.tsv"
def pheno_sigsnps(method, adjust): return f"{MOD_GWAS}tables/methods/{method}/{method}_pvalues_K{K_BEST}_sig_snps_{adjust}.tsv"
def pheno_manhattan(method, trait, adjust): return f"{MOD_GWAS}plots/manhattan/{method}/manhattan_{trait}_K{K_BEST}_{adjust}.png"
def pheno_qq(method, trait, adjust): return f"{MOD_GWAS}plots/manhattan/{method}/qq_{trait}_K{K_BEST}_{adjust}.png"

# Per-trait DROP mode paths
def pheno_trait_vcf(trait): return f"{WORK_FILT}phenotypes/{trait}/{VCF_BASE}.vcf"
def pheno_trait_tped(trait): return f"{WORK_FILT}phenotypes/{trait}/emmax/{VCF_BASE}.tped"
def pheno_trait_tfam(trait): return f"{WORK_FILT}phenotypes/{trait}/emmax/{VCF_BASE}.tfam"
def pheno_trait_kinship(trait): return f"{WORK_FILT}phenotypes/{trait}/emmax/{VCF_BASE}.aBN.kinf"
def pheno_trait_pvalues(method, trait): return f"{MOD_GWAS}tables/methods/{method}/{trait}_pvalues_K{K_BEST}.tsv"

#=============================================================================
# PHASE 4: Source-indexed config for unified GEA/GWAS downstream rules
#=============================================================================
# Keys match output directory names. "gea" and "gwas" will be set in Sub-step C.
ASSOC_SOURCES = {}

if K_BEST is not None and GEA_CONFIGS:
    ASSOC_SOURCES["GEA"] = {
        "mod":             MOD_GEA,
        "logdir":          f"{LOGDIR}GEA/",
        "configs":         GEA_CONFIGS,
        "method_regex":    GEA_METHOD_REGEX,
        "trait_regex":     r"bio_\d+",
        "predictors":      PREDICTORS_SELECTED,
        "snp_distance":    SNP_DISTANCE,
        "region_distance": REGION_DISTANCE,
        "region_auto":     REGION_DISTANCE_AUTO,
        "promoter_length": PROMOTER_LENGTH,
        "combine_method":  SIGSNPS_METHOD,
        "combine_gap":     SIGSNPS_GAP,
        "pvalues_fn":      assoc_pvalues,
        "sigsnps_fn":      assoc_sigsnps,
        "block_mode":      GEA_BLOCK_MODE,
    }

if K_BEST is not None and GWAS_CONFIGS:
    ASSOC_SOURCES["GWAS"] = {
        "mod":             MOD_GWAS,
        "logdir":          f"{LOGDIR}GWAS/",
        "configs":         GWAS_CONFIGS,
        "method_regex":    GWAS_METHOD_REGEX,
        "trait_regex":     r"[a-zA-Z]\w*",
        "predictors":      PHENO_PREDICTORS,
        "snp_distance":    PHENO_SNP_DISTANCE,
        "region_distance": PHENO_REGION_DISTANCE,
        "region_auto":     PHENO_REGION_DISTANCE_AUTO,
        "promoter_length": PHENO_PROMOTER_LENGTH,
        "combine_method":  PHENO_COMBINE_METHOD,
        "combine_gap":     PHENO_COMBINE_GAP,
        "pvalues_fn":      pheno_pvalues,
        "sigsnps_fn":      pheno_sigsnps,
        "block_mode":      GWAS_BLOCK_MODE,
    }

def _src(source, key):
    """Look up a source-indexed config value in unified downstream rule bodies."""
    return ASSOC_SOURCES[source][key]

# Wildcard regex matching all configured sources (used in _assoc_downstream.smk)
SOURCE_REGEX = "|".join(ASSOC_SOURCES.keys()) if ASSOC_SOURCES else "gea"
# Trait wildcard: union of GEA (bio_\d+) and GWAS (identifier) patterns.
# {source} in the output path disambiguates which branch a match belongs to.
TRAIT_REGEX_ANY = r"(bio_\d+|[a-zA-Z]\w*)"

# Source classification by block_mode — used in both _assoc_downstream.smk and ld_blocks.smk
_BLOCK_SOURCES = sorted([s for s, cfg in ASSOC_SOURCES.items() if cfg.get("block_mode", "snp") == "block"])
_SNP_SOURCES   = sorted([s for s, cfg in ASSOC_SOURCES.items() if cfg.get("block_mode", "snp") == "snp"])
if _BLOCK_SOURCES:
    _BLOCK_SOURCE_REGEX      = "|".join(_BLOCK_SOURCES)
    _ALL_BLOCK_METHODS_REGEX = "|".join(sorted({
        m for s in _BLOCK_SOURCES for m in ASSOC_SOURCES[s]["configs"].keys()
    }))
    _SNP_SOURCE_REGEX = "|".join(_SNP_SOURCES) if _SNP_SOURCES else "NOSOURCE_snp"
else:
    _SNP_SOURCE_REGEX = SOURCE_REGEX  # all sources are snp-mode, no partitioning needed

def assoc_out(source, key):
    """Return the per-source output path for a logical downstream key."""
    _templates = {
        "selected_snps":              "tables/selected_snps.tsv",
        "selected_blocks":            "tables/selected_blocks.tsv",
        "regions_per_trait":          "tables/regions_per_trait.tsv",
        "regions_combined":           "tables/regions_combined.tsv",
        "genes_per_region":           "tables/genes_per_region.tsv",
        "genes_per_region_collapsed": "tables/genes_per_region_collapsed.tsv",
        "genes_combined":             "tables/genes_combined.tsv",
        "manhattan_combined_png":     f"plots/manhattan/combined/manhattan_combined_K{K_BEST}.png",
        "manhattan_combined_svg":     f"plots/manhattan/combined/manhattan_combined_K{K_BEST}.svg",
        "qq_combined_png":            f"plots/manhattan/combined/qq_combined_K{K_BEST}.png",
        "qq_combined_svg":            f"plots/manhattan/combined/qq_combined_K{K_BEST}.svg",
        "manhattan_combined_bg":      f"plots/manhattan/combined/manhattan_combined_K{K_BEST}_background.png",
        "manhattan_combined_coords":  f"plots/manhattan/combined/manhattan_combined_K{K_BEST}_coords.json",
    }
    return f"{ASSOC_SOURCES[source]['mod']}{_templates[key]}"

def assoc_block_pvalues(source, method):
    """Per-method WZA block-level p-value table."""
    return f"{ASSOC_SOURCES[source]['mod']}tables/methods/{method}/{method}_block_pvalues_K{K_BEST}.tsv"

def assoc_sig_blocks(source, method, adjust):
    """Significant blocks for a (source, method, adjust) combination."""
    return f"{ASSOC_SOURCES[source]['mod']}tables/methods/{method}/{method}_block_pvalues_K{K_BEST}_sig_blocks_{adjust}.tsv"

#=============================================================================
# RULE FACTORIES — per-engine helpers for dynamic rule declaration
#=============================================================================
# These are called from the top-level for-loops in association.smk /
# phenotype_assoc.smk (same loop idiom). Each helper returns a
# plain Python dict or string; the rule: blocks live in the .smk files.

def _gea_inputs(engine):
    """Input dict for a GEA non-GAPIT rule of a given engine."""
    if engine == "emmax":
        return dict(
            vcf      = W['vcf_filt'],
            tped     = W['assoc_tped'],
            kinship  = W['assoc_kinship'],
            traits   = O['climate_site_scaled'],
            covariates = W['pca_projections'],
            metadata = O['metadata'],
        )
    elif engine == "lfmm":
        return dict(
            lfmm_ld  = W['lfmm_imp'],
            lfmm_full = W['lfmm_imp_full'],
            climate  = O['climate_site_scaled'],
            vcfsnp   = W['vcfsnp_full'],
        )
    raise ValueError(f"_gea_inputs: unknown engine '{engine}'")

def _gea_params(engine, method):
    """Params dict for a GEA non-GAPIT rule of a given engine."""
    if engine == "emmax":
        return dict(
            k           = K_BEST,
            predictors  = PREDICTORS_SELECTED,
            inter_dir   = INTER,
            tables_dir  = f"{MOD_GEA}tables/methods/EMMAX/",
            tped_prefix = f"{WORK_FILT}emmax/{VCF_BASE}",
        )
    elif engine == "lfmm":
        return dict(
            k          = K_BEST,
            predictors = PREDICTORS_SELECTED,
            tables_dir = f"{MOD_GEA}tables/methods/LFMM/",
        )
    raise ValueError(f"_gea_params: unknown engine '{engine}'")


def _pheno_a_inputs(engine):
    """Input dict for phenotype Path A (MEAN/MEDIAN) rule of a given engine."""
    if engine == "emmax":
        return dict(
            vcf        = W['vcf_filt'],
            tped       = W['pheno_tped'],
            kinship    = W['pheno_kinship'],
            pca        = W['pca_projections'],
            phenotypes = W['pheno_all_phenotypes'],
        )
    raise ValueError(f"_pheno_a_inputs: unknown engine '{engine}'")

def _pheno_a_params(engine, method):
    """Params dict for phenotype Path A (MEAN/MEDIAN) rule of a given engine."""
    if engine == "emmax":
        return dict(
            tped_prefix = f"{WORK_FILT}phenotypes/emmax/{VCF_BASE}",
            k           = K_BEST,
            tables_dir  = f"{MOD_GWAS}tables/methods/",
        )
    raise ValueError(f"_pheno_a_params: unknown engine '{engine}'")


#=============================================================================
# CREATE DIRECTORIES
#=============================================================================
dirs_to_create = [
    WORK, WORK_FILT, WORK_LD, INTER,
    f"{INTER}samples/", f"{INTER}annotation/", f"{INTER}flags/", f"{INTER}qc/",
    LOGDIR,
    # Module directories
    f"{MOD_PROCESSING}tables/", f"{MOD_PROCESSING}plots/",
    f"{MOD_PRESTRUCT}plots/", f"{MOD_PRESTRUCT}tables/",
    *[f"{MOD_PRESTRUCT}plots/K{k}/" for k in k_range(K_START, K_END)],
    *[f"{MOD_PRESTRUCT}tables/K{k}/" for k in k_range(K_START, K_END)],
    *([ f"{MOD_CLIMATE}plots/", f"{MOD_CLIMATE}tables/present/", f"{MOD_CLIMATE}rasters/present/" ] if CLIMATE_ENABLED else []),
    f"{MOD_STRUCT}plots/piemap/", f"{MOD_STRUCT}tables/",
    # Log subdirectories
    f"{LOGDIR}processing/", f"{LOGDIR}prestructure/", f"{LOGDIR}structure/",
    f"{LOGDIR}gea/", f"{LOGDIR}GEA/", f"{LOGDIR}maladaptation/",
]

# Add association directories for each method
for method in GEA_CONFIGS:
    dirs_to_create.append(f"{MOD_GEA}plots/manhattan/{method}/")
    dirs_to_create.append(f"{MOD_GEA}tables/methods/{method}/")

# Add enrichment directories
dirs_to_create.append(f"{MOD_GEA}tables/")
dirs_to_create.append(f"{MOD_GEA}tables/enrichment/")
dirs_to_create.append(f"{MOD_GEA}plots/")
dirs_to_create.append(f"{MOD_GEA}plots/enrichment/")
dirs_to_create.append(f"{MOD_GEA}plots/manhattan/combined/")
dirs_to_create.append(f"{INTER}enrichment/gea/")

# Add GAPIT directories
if GEA_GAPIT_CONFIGS:
    dirs_to_create.append(f"{WORK_FILT}gapit/")
    dirs_to_create.append(f"{INTER}gapit/gea/")
    dirs_to_create.append(f"{MOD_GEA}GAPIT_native_output/")
    for model in GEA_GAPIT_CONFIGS:
        dirs_to_create.append(f"{MOD_GEA}GAPIT_native_output/{model}/")

# Add maladaptation directories (require climate)
if CLIMATE_ENABLED:
    for _m in ACTIVE_MALA_METHODS:
        dirs_to_create.append(f"{INTER}{_m}/")
    dirs_to_create.append(f"{MOD_CLIMATE}rasters/future/")
    dirs_to_create.append(f"{MOD_CLIMATE}tables/future/")

# Add phenotype association directories
if GWAS_CONFIGS:
    dirs_to_create.append(f"{MOD_GWAS}tables/")
    dirs_to_create.append(f"{MOD_GWAS}tables/enrichment/")
    dirs_to_create.append(f"{MOD_GWAS}plots/")
    dirs_to_create.append(f"{MOD_GWAS}plots/enrichment/")
    dirs_to_create.append(f"{MOD_GWAS}plots/piemap/")
    dirs_to_create.append(f"{MOD_GWAS}plots/manhattan/combined/")
    dirs_to_create.append(f"{INTER}enrichment_gwas/")
    dirs_to_create.append(f"{LOGDIR}gwas/")
    dirs_to_create.append(f"{LOGDIR}GWAS/")
    for method in GWAS_CONFIGS:
        dirs_to_create.append(f"{MOD_GWAS}plots/manhattan/{method}/")
        dirs_to_create.append(f"{MOD_GWAS}tables/methods/{method}/")
    if PHENO_MISSING == 'DROP':
        for trait in PHENO_TRAITS:
            dirs_to_create.append(f"{WORK_FILT}phenotypes/{trait}/emmax/")
    if GWAS_GAPIT_CONFIGS:
        dirs_to_create.append(f"{INTER}gapit/gwas/")
        dirs_to_create.append(f"{MOD_GWAS}GAPIT_native_output/")
        for model in GWAS_GAPIT_CONFIGS:
            dirs_to_create.append(f"{MOD_GWAS}GAPIT_native_output/{model}/")
        if not GEA_GAPIT_CONFIGS:
            # GAPIT numeric dir needed even when only pheno uses GAPIT
            dirs_to_create.append(f"{WORK_FILT}gapit/")

# Add overlapping analysis directories
if GEA_CONFIGS and GWAS_CONFIGS:
    dirs_to_create.append(f"{MOD_GEAXGWAS}tables/")
    dirs_to_create.append(f"{MOD_GEAXGWAS}tables/enrichment/")
    dirs_to_create.append(f"{MOD_GEAXGWAS}plots/")
    dirs_to_create.append(f"{MOD_GEAXGWAS}plots/enrichment/")
    dirs_to_create.append(f"{INTER}enrichment_gea_x_gwas/")
    dirs_to_create.append(f"{LOGDIR}gea_x_gwas/")
    dirs_to_create.append(f"{LOGDIR}GEAxGWAS/")

# Add LD decay directories
if K_BEST is not None:
    dirs_to_create.append(W['ld_decay_work'])
    dirs_to_create.append(W['ld_decay_sample_lists'])
    dirs_to_create.append(W['ld_decay_chr_vcfs'])
    dirs_to_create.append(W['ld_decay_stat_gw'])
    dirs_to_create.append(W['ld_decay_stat_chr'])

# Add pop stats directories if enabled
if CALC_POP_STATS:
    dirs_to_create.append(f"{MOD_STRUCT}plots/pop_stats/")
    dirs_to_create.append(f"{MOD_STRUCT}tables/pop_stats/")

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

def _targets_for_assoc_source(source):
    """Build the shared downstream target list for a GEA or GWAS source."""
    src = ASSOC_SOURCES[source]
    traits = get_predictors_list() if source == "GEA" else PHENO_TRAITS
    block_mode = src.get("block_mode", "snp")
    targets = []
    for method, adjust in src["configs"].items():
        targets.append(src["pvalues_fn"](method))
        if block_mode == "snp":
            # Standard: per-SNP sig SNPs file
            targets.append(src["sigsnps_fn"](method, adjust))
        else:
            # Block mode: WZA block p-values + sig blocks per method
            targets.append(assoc_block_pvalues(source, method))
            targets.append(assoc_sig_blocks(source, method, adjust))
        for trait in traits:
            targets.append(f"{src['mod']}plots/manhattan/{method}/manhattan_{trait}_K{K_BEST}_{adjust}.png")
            targets.append(f"{src['mod']}plots/manhattan/{method}/qq_{trait}_K{K_BEST}_{adjust}.png")
    # Shared downstream: same output paths regardless of mode
    for key in (
        "selected_snps", "regions_per_trait", "regions_combined",
        "genes_per_region", "genes_per_region_collapsed", "genes_combined",
        "manhattan_combined_png", "qq_combined_png",
    ):
        targets.append(assoc_out(source, key))
    if block_mode == "block":
        targets.append(assoc_out(source, "selected_blocks"))
    return targets


def get_targets(mode):
    if mode == 'processing':
        targets = [
            W['samples_missing_stats'], W['samples_removed'],  # Sample missingness outputs
            W['vcf_filt'], W['vcf_ld'], W['geno'], W['lfmm'], O['metadata'],
            # QC outputs
            O['qc_raw_summary'], O['qc_filtering_summary'], O['qc_sample_het'],
            O['qc_plot_sample_miss'], O['qc_plot_het_miss'],
            O['qc_plot_maf'], O['qc_plot_snp_miss'],
            O['qc_plot_attrition'], O['qc_plot_snp_density'],
            W['summary_done']
        ]
        if GFF:
            targets.append(W['gff_normalized'])
        if HAS_FORMAT_DP:
            targets.append(O['qc_plot_depth'])
        return targets
    
    elif mode == 'prestructure':
        ks = k_range(K_START, K_END)
        return (
            [clusters_table(k) for k in ks] +
            [structure_plot(k) for k in ks] +
            [pca_struct_plot(k) for k in ks] +
            [pop_diff_plot(k) for k in ks] +
            [O['pca'], O['tracy'], O['cross_entropy'], W['summary_done']]
        )
    
    elif mode == 'structure':
        check_numeric(K_BEST, 'K_BEST')

        # Imputed data (always needed)
        targets = [W['lfmm_imp'], W['vcf_imp']]

        # Climate-dependent targets
        if CLIMATE_ENABLED:
            targets += [O['climate_site'], O['climate_site_scaled'], O['climate_all'], W['climate_raster']]
            targets += [O['climate_invariant']]
            targets += [DENSITY_PLOT_COMBINED]
            targets += [O['corr_heatmap']]
            # Simple PieMaps for ALL 19 BIO variables (exploration step)
            targets += [piemap_notrait(bio) for bio in ALL_BIO]
        if META_HAS_PHENO:
            targets += [DENSITY_PLOT_PHENOTYPES]

        # Population statistics (requires >= 3 samples per population)
        if CALC_POP_STATS:
            targets += [O['tajima'], O['pi_div'], O['ibd_raw'], O['ibd_pairs']]
            targets += [O['amova'], O['amova_plot']]
            if CLIMATE_ENABLED:
                targets += [O['mantel']]
                # PieMaps with trait overlays for ALL 19 BIO variables
                targets += [piemap_tajima(bio) for bio in ALL_BIO]
                targets += [piemap_diversity(bio) for bio in ALL_BIO]

        # LD decay — table always produced
        targets.append(O['ld_decay_table'])
        if LD_DECAY_SCOPE in ('genome_wide', 'both'):
            targets += [O['ld_decay_plot_gw'], O['ld_decay_plot_gw_svg']]
        if LD_DECAY_SCOPE in ('per_chromosome', 'both'):
            targets += [O['ld_decay_plot_chr'], O['ld_decay_plot_chr_svg']]

        targets.append(W['summary_done'])
        return targets

    elif mode == 'gea':
        if not CLIMATE_ENABLED:
            raise ValueError("association mode requires climate.enabled: true (GEA uses climate predictors as traits)")
        check_numeric(K_BEST, 'K_BEST')
        if not get_predictors_list():
            raise ValueError("climate.predictors must be set for association mode")
        if not GEA_CONFIGS:
            raise ValueError("GEA.configs must be set for association mode")
        if GFF:
            check_file_exists(INDIR, GFF, 'GFF')

        targets = _targets_for_assoc_source("GEA")
        targets.append(W['summary_done'])
        return targets


    elif mode == 'maladaptation':
        if not CLIMATE_ENABLED:
            raise ValueError("maladaptation mode requires climate.enabled: true")
        check_numeric(K_BEST, 'K_BEST')
        check_numeric(SSP, 'SSP')
        check_numeric(NTREE, 'NTREE')
        check_float(COR_THRESHOLD, 'COR_THRESHOLD')
        if not MODELS_LIST:
            raise ValueError("MODELS must be set for maladaptation mode")

        targets = [
            # Future climate (method-agnostic)
            O['climate_future_site'],
            O['climate_future_all'],
            W['climate_future_raster'],
            O['density_future'],
        ]

        for method in ACTIVE_MALA_METHODS:
            targets += [
                mala_model(method, GF_RUN_LABEL, SPATIAL_TAG, 'adaptive'),
                mala_offset_map_values(method, GF_RUN_LABEL, SPATIAL_TAG),
                mala_offset_site_values(method, GF_RUN_LABEL, SPATIAL_TAG),
                mala_cumimp(method, GF_RUN_LABEL, SPATIAL_TAG),
                mala_importance(method, GF_RUN_LABEL, SPATIAL_TAG),
                mala_offset_piemap(method, GF_RUN_LABEL, SPATIAL_TAG, 'notrait'),
            ]
            if MALADAPTATION_METHODS[method]['supports_random_model'] and GF_RANDOM_MODEL:
                targets.append(mala_model(method, GF_RUN_LABEL, SPATIAL_TAG, 'random'))
            if CALC_POP_STATS:
                targets += [
                    mala_offset_piemap(method, GF_RUN_LABEL, SPATIAL_TAG, 'tajima_d'),
                    mala_offset_piemap(method, GF_RUN_LABEL, SPATIAL_TAG, 'pi_diversity'),
                ]

        targets.append(W['summary_done'])
        return targets

    elif mode == 'gwas':
        check_numeric(K_BEST, 'K_BEST')
        if not PHENO_TRAITS:
            raise ValueError("No phenotype columns found in metadata (expected columns 5+ after site, sample, latitude, longitude)")
        if not GWAS_CONFIGS:
            raise ValueError("GWAS.configs must be set for association_phenotypes mode")
        if GFF:
            check_file_exists(INDIR, GFF, 'GFF')

        targets = [O['pheno_missing_summary']] + _targets_for_assoc_source("GWAS")

        # Phenotype piemaps (require climate raster for background)
        if CLIMATE_ENABLED:
            for trait in PHENO_TRAITS:
                targets.append(f"{MOD_GWAS}plots/piemap/phenomap_{trait}.png")

        targets.append(W['summary_done'])
        return targets

    elif mode == 'gea_x_gwas':
        check_numeric(K_BEST, 'K_BEST')
        has_gea  = bool(GEA_CONFIGS) and CLIMATE_ENABLED
        has_gwas = bool(GWAS_CONFIGS)
        if not has_gea and not has_gwas:
            raise ValueError(
                "overlapping mode requires at least one of: "
                "association configs (GEA) or phenotype_association configs (GWAS)"
            )
        if GFF:
            check_file_exists(INDIR, GFF, 'GFF')

        # Pairwise trait overlap (works with either or both sources)
        targets = [
            O['pairwise_collapsed_snps'],
            O['pairwise_overlap_table'],
        ]

        # Miami plot (only when both GEA and GWAS available)
        if has_gea and has_gwas:
            targets.extend([
                O['overlap_miami'], O['overlap_miami_svg'],
            ])

        targets.append(W['summary_done'])
        return targets

    elif mode is None:
        raise ValueError("Specify mode: --config mode=processing or mode=prestructure or mode=structure or mode=gea or mode=gwas or mode=gea_x_gwas or mode=maladaptation")
    else:
        raise ValueError(f"Unknown mode: {mode}")

# Helper: provide LD decay table as input when region_distance is "auto"
def ld_decay_input(auto_flag):
    """Return LD decay table path as input dependency when auto mode is enabled."""
    return O.get('ld_decay_table', []) if auto_flag else []

#=============================================================================
# MAIN RULE
