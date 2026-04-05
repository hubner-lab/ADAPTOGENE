# ADAPTOGENE Pipeline - Refactored
# vim: filetype=python
import os
import subprocess
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
# CONFIG MIGRATION (backward compatibility with flat YAML keys)
#=============================================================================
def _migrate_config(cfg):
    """Detect old flat config keys and map them to nested structure.
    Allows old configs to keep working with deprecation warnings."""
    # If 'input' key exists as a dict, assume new-style config
    if isinstance(cfg.get('input'), dict):
        return cfg

    import sys
    print("WARNING: Detected old flat-key config format. Migrating automatically.", file=sys.stderr)
    print("WARNING: Please update your config to nested YAML format. See config_SIMDATA.yaml for reference.", file=sys.stderr)

    migrated = dict(cfg)

    # Mapping: old flat key -> (nested group, nested key)
    _FLAT_TO_NESTED = {
        'INPUT_DIR': ('input', 'dir'),
        'INPUT_VCF': ('input', 'vcf'),
        'INPUT_METADATA': ('input', 'metadata'),
        'INPUT_GFF': ('input', 'gff'),
        'PROJECT_NAME': ('project_name', None),
        'CPU': ('cpu', None),
        'FILTER_MAF': ('filter', 'maf'),
        'FILTER_SNP_MISS': ('filter', 'snp_miss'),
        'FILTER_SAMPLE_MISS': ('filter', 'sample_miss'),
        'LD_WINDOW': ('ld', 'window'),
        'LD_STEP': ('ld', 'step'),
        'LD_R2': ('ld', 'r2'),
        'SNMF_K_START': ('snmf', 'k_start'),
        'SNMF_K_END': ('snmf', 'k_end'),
        'SNMF_K_BEST': ('snmf', 'k_best'),
        'SNMF_PLOIDY': ('snmf', 'ploidy'),
        'SNMF_REPEATS': ('snmf', 'repeats'),
        'MAP_CROP_REGION': ('map', 'climate_extent'),
        'MAP_GAP': ('map', 'gap'),
        'MAP_RESOLUTION': ('map', 'resolution'),
        'MAP_REGIONMAP_EXTENT': ('map', 'zoom_extent'),
        'CLIMATE_PREDICTORS': ('climate', 'predictors'),
        'POP_CALC_STATS': ('pop', 'calc_stats'),
        'POP_WINDOW_SIZE': ('pop', 'window_size'),
        'POP_CUSTOM_TRAIT': ('pop', 'custom_trait_file'),
        'PIEMAP_ALPHA': ('piemap', 'alpha'),
        'PIEMAP_SHOW_LABELS': ('piemap', 'show_labels'),
        'PIEMAP_LABEL_SIZE': ('piemap', 'label_size'),
        'PIEMAP_PIE_SCALE': ('piemap', 'pie_scale'),
        'PIEMAP_USE_POINTS': ('piemap', 'use_points'),
        'ASSOC_CONFIGS': ('association', 'configs'),
        'ASSOC_COMBINE_METHOD': ('association', 'combine_method'),
        'ASSOC_COMBINE_GAP': ('association', 'combine_gap'),
        'ASSOC_REGION_DISTANCE': ('association', 'region_distance'),
        'ASSOC_PROMOTER_LENGTH': ('association', 'promoter_length'),
        'ASSOC_GO_FIELD': ('association', 'go_field'),
        'GFF_FEATURE': ('gff', 'feature'),
        'GFF_GENE_NAME': ('gff', 'gene_name'),
        'GFF_BIOTYPE': ('gff', 'biotype'),
        'ENRICHMENT_TOP_TERMS': ('enrichment', 'top_terms'),
        'ENRICHMENT_PLOT_WIDTH': ('enrichment', 'plot_width'),
        'ENRICHMENT_PLOT_HEIGHT': ('enrichment', 'plot_height'),
        'ENRICHMENT_CNET_LABEL': ('enrichment', 'cnet_label'),
        'FUTURE_SSP': ('future', 'ssp'),
        'FUTURE_YEAR': ('future', 'year'),
        'FUTURE_MODELS': ('future', 'models'),
        'GF_NTREE': ('gradient_forest', 'ntree'),
        'GF_COR_THRESHOLD': ('gradient_forest', 'cor_threshold'),
        'GF_PCNM': ('gradient_forest', 'spatial_correction'),
        'GF_SUFFIX': ('gradient_forest', 'run_label'),
        'GF_RANDOM_MODEL': ('gradient_forest', 'random_model'),
        'PHENO_ASSOC_CONFIGS': ('phenotype_association', 'configs'),
        'PHENO_MISSING_STRATEGY': ('phenotype_association', 'missing_strategy'),
        'PHENO_COMBINE_METHOD': ('phenotype_association', 'combine_method'),
        'PHENO_COMBINE_GAP': ('phenotype_association', 'combine_gap'),
        'PHENO_REGION_DISTANCE': ('phenotype_association', 'region_distance'),
        'PHENO_PROMOTER_LENGTH': ('phenotype_association', 'promoter_length'),
        'OVERLAP_REGION_DISTANCE': ('overlap', 'region_distance'),
        'OVERLAP_PROMOTER_LENGTH': ('overlap', 'promoter_length'),
        'HAPLOTYPE_SCAN_REGIONS_SOURCE': ('haplotype', 'scan', 'regions_source'),
        'HAPLOTYPE_SCAN_REGIONS_FILE': ('haplotype', 'scan', 'regions_file'),
        'HAPLOTYPE_SCAN_TOP_REGIONS': ('haplotype', 'scan', 'top_regions'),
        'HAPLOTYPE_SCAN_MIN_SNPS': ('haplotype', 'scan', 'min_snps'),
        'HAPLOTYPE_SCAN_MGMIN': ('haplotype', 'scan', 'min_group_size'),
        'HAPLOTYPE_SCAN_MINHAP': ('haplotype', 'scan', 'min_haplotype_size'),
        'HAPLOTYPE_SCAN_EPSILON_RANGE': ('haplotype', 'scan', 'epsilon_range'),
        'HAPLOTYPE_SCAN_METADATA_TYPE': ('haplotype', 'scan', 'metadata_type'),
        'HAPLOTYPE_EPSILON_SELECTED': ('haplotype', 'epsilon_selected'),
    }

    # Also migrate old METHOD/ADJUST/THRESHOLD keys in configs lists
    def _migrate_configs_list(configs_list):
        """Convert [{METHOD:..., ADJUST:..., THRESHOLD:...}] to [{method:..., adjust:..., threshold:...}]."""
        if not configs_list:
            return configs_list
        migrated_list = []
        for entry in configs_list:
            if 'METHOD' in entry:
                migrated_list.append({
                    'method': entry['METHOD'],
                    'adjust': entry['ADJUST'],
                    'threshold': entry['THRESHOLD']
                })
            else:
                migrated_list.append(entry)
        return migrated_list

    for flat_key, path in _FLAT_TO_NESTED.items():
        if flat_key not in cfg:
            continue
        val = cfg[flat_key]

        if path[1] is None:
            # Top-level key (e.g., project_name, cpu)
            migrated[path[0]] = val
        elif len(path) == 2:
            group, key = path
            if group not in migrated or not isinstance(migrated.get(group), dict):
                migrated[group] = {}
            migrated[group][key] = val
        elif len(path) == 3:
            group, sub, key = path
            if group not in migrated or not isinstance(migrated.get(group), dict):
                migrated[group] = {}
            if sub not in migrated[group] or not isinstance(migrated[group].get(sub), dict):
                migrated[group][sub] = {}
            migrated[group][sub][key] = val

    # Migrate ASSOC_CONFIGS list entries (METHOD -> method)
    assoc_cfg = migrated.get('association', {})
    if 'configs' in assoc_cfg:
        assoc_cfg['configs'] = _migrate_configs_list(assoc_cfg['configs'])
    pheno_cfg = migrated.get('phenotype_association', {})
    if 'configs' in pheno_cfg:
        pheno_cfg['configs'] = _migrate_configs_list(pheno_cfg['configs'])

    return migrated

config = _migrate_config(config)

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
INDIR = '/pipeline/' + config['input']['dir']
PROJECT = config['project_name']
OUTDIR = f'/pipeline/{PROJECT}_results/'
LOGDIR = f'/pipeline/{PROJECT}_logs/'

CPU = config.get('cpu', max(1, os.cpu_count() - 2)); check_numeric(CPU, 'cpu')
VCF_RAW = config['input']['vcf']; check_file_exists(INDIR, VCF_RAW, 'input.vcf')
SAMPLES = config['input']['metadata']; check_file_exists(INDIR, SAMPLES, 'input.metadata')
VCF_BASE = get_vcf_basename(VCF_RAW)

# FILTER parameters
MAF = config['filter']['maf']; check_float(MAF, 'filter.maf')
MISS = config['filter']['snp_miss']; check_float(MISS, 'filter.snp_miss')
SAMPLE_MISS = _cfg('filter', 'sample_miss', 0.5); check_float(SAMPLE_MISS, 'filter.sample_miss')
MIN_DEPTH    = _cfg('filter', 'min_depth', None)
MAX_DEPTH    = _cfg('filter', 'max_depth', None)
HET_OUTLIER_SD = _cfg('filter', 'het_outlier_sd', None)

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
LD_WIN = config['ld']['window']; check_numeric(LD_WIN, 'ld.window')
LD_STEP = config['ld']['step']; check_numeric(LD_STEP, 'ld.step')
LD_R2 = config['ld']['r2']; check_float(LD_R2, 'ld.r2')

# SNMF parameters
K_START = config['snmf']['k_start']; check_numeric(K_START, 'snmf.k_start')
K_END = config['snmf']['k_end']; check_numeric(K_END, 'snmf.k_end')
PLOIDY = 2  # diploid only; haploid support requires changes in lfmm2vcf.R, amova.R, plink
REPEAT = config['snmf']['repeats']; check_numeric(REPEAT, 'snmf.repeats')
K_BEST = int(_cfg('snmf', 'k_best', None)) if _cfg('snmf', 'k_best', None) is not None else None
SNMF_PROJECT_MODE = 'new'  # LEA project mode: 'new' for fresh runs, 'continue' to resume

# MAP parameters
CROP_REGION = _cfg('map', 'climate_extent', 'auto')
GAP = _cfg('map', 'gap', 0.5)
RESOLUTION = _cfg('map', 'resolution', 0.5)
REGIONMAP_EXTENT = _cfg('map', 'zoom_extent', 'NULL')
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
CLIMATE_ENABLED = _cfg_bool('climate', 'enabled', True)
PREDICTORS_SELECTED = _cfg('climate', 'predictors', '') if CLIMATE_ENABLED else ''

# All 19 WorldClim bioclimatic variables — always available in the raster stack
ALL_BIO     = [f"bio_{i}" for i in range(1, 20)]
ALL_BIO_STR = ",".join(ALL_BIO)

# POP parameters
CALC_POP_STATS = _cfg_bool('pop', 'calc_stats', False)
METRICS_WINSIZE = _cfg('pop', 'window_size', 10000)
CUSTOM_TRAIT = _cfg('pop', 'custom_trait_file', 'NULL')

# PIEMAP parameters (palette fixed to viridis plasma)
PIEMAP_ALPHA = _cfg('piemap', 'alpha', 0.6)
PIEMAP_SHOW_LABELS = 'T' if _cfg_bool('piemap', 'show_labels', False) else 'F'
PIEMAP_LABEL_SIZE = _cfg('piemap', 'label_size', 10)
PIEMAP_PIE_SCALE = _cfg('piemap', 'pie_scale', 1.0)
PIEMAP_USE_POINTS = 'T' if _cfg_bool('piemap', 'use_points', False) else 'F'

# LD DECAY parameters
_ld_decay = config.get('ld_decay', {})
LD_DECAY_GROUP_BY = _ld_decay.get('group_by', 'site')
check_in_list(LD_DECAY_GROUP_BY, ['site', 'cluster'], 'ld_decay.group_by')
LD_DECAY_MIN_SAMPLES = int(_ld_decay.get('min_samples', 10))
LD_DECAY_MAX_DISTANCE = int(_ld_decay.get('max_distance', 500))
LD_DECAY_SCOPE = _ld_decay.get('scope', 'both')
check_in_list(LD_DECAY_SCOPE, ['genome_wide', 'per_chromosome', 'both'], 'ld_decay.scope')
LD_DECAY_TAG = f"grp{LD_DECAY_GROUP_BY}_min{LD_DECAY_MIN_SAMPLES}_dist{LD_DECAY_MAX_DISTANCE}"

# GFF parameters
GFF = _cfg('input', 'gff', '')
GFF_FEATURE = _cfg('gff', 'feature', 'mRNA')
GFF_GENE_NAME = _cfg('gff', 'gene_name', 'description')
GFF_BIOTYPE = _cfg('gff', 'biotype', 'biotype')

# ASSOCIATION parameters
_assoc = config.get('association', {})
PROMOTER_LENGTH = _assoc.get('promoter_length', 10000)
SIGSNPS_METHOD = 'All'  # pipeline always uses All; combine strategy moved to gradient_forest config
SIGSNPS_GAP = _assoc.get('combine_gap', 100000)
SNP_DISTANCE = SIGSNPS_GAP  # overlap annotation distance
_region_distance_raw = _assoc.get('region_distance', 1000000)
REGION_DISTANCE_AUTO = (str(_region_distance_raw).lower() == 'auto')
REGION_DISTANCE = 'auto' if REGION_DISTANCE_AUTO else int(_region_distance_raw)
GO_FIELD = _assoc.get('go_field', 'NULL')

# ENRICHMENT parameters
_enrich = config.get('enrichment', {})
ENRICHMENT_TOP_TERMS = _enrich.get('top_terms', 20)
ENRICHMENT_PLOT_WIDTH = _enrich.get('plot_width', 12)
ENRICHMENT_PLOT_HEIGHT = _enrich.get('plot_height', 10)
ENRICHMENT_CNET_LABEL = _enrich.get('cnet_label', 'gene_id')

# Note: regionplot mode deprecated — regionplots generated on-demand in Shiny app
GENES_TO_HIGHLIGHT = config.get('regionplot', {}).get('genes', 'all')

# FUTURE parameters
_future = config.get('future', {})
SSP = _future.get('ssp', '585')
YEAR = _future.get('year', '2061-2080')
MODELS_STR = _future.get('models', '')
MODELS_LIST = [m.strip() for m in MODELS_STR.split(',') if m.strip()] if MODELS_STR else []

# GRADIENT FOREST parameters
_gf = config.get('gradient_forest', {})
NTREE = _gf.get('ntree', '500')
COR_THRESHOLD = _gf.get('cor_threshold', '0.5')
PCNM = _gf.get('spatial_correction', 'with')
GF_SUFFIX = _gf.get('run_label', _gf.get('suffix', ''))
GF_RANDOM_MODEL = _gf.get('random_model', True)
GF_COMBINE_METHOD = _gf.get('combine_method', 'All')
GF_COMBINE_GAP = int(_gf.get('combine_gap', 100000))
PCNM_TAG = 'PCNM' if PCNM == 'with' else 'noPCNM'

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

ASSOC_CONFIGS = parse_association_configs(_assoc.get('configs', []))

# GAPIT model detection — auto-route methods to GAPIT or legacy (EMMAX/LFMM) engine
GAPIT_MODELS = {'GLM', 'MLM', 'CMLM', 'ECMLM', 'SUPER', 'MLMM', 'FarmCPU', 'BLINK'}

def split_configs_by_engine(configs):
    """Split configs into GAPIT and non-GAPIT (legacy) engines."""
    gapit = {m: a for m, a in configs.items() if m in GAPIT_MODELS}
    other = {m: a for m, a in configs.items() if m not in GAPIT_MODELS}
    return gapit, other

ASSOC_GAPIT_CONFIGS, ASSOC_OTHER_CONFIGS = split_configs_by_engine(ASSOC_CONFIGS)

# Dynamic wildcard constraint regex for association methods
ASSOC_METHOD_REGEX = '|'.join(ASSOC_CONFIGS.keys()) if ASSOC_CONFIGS else 'EMMAX'

# PHENOTYPE ASSOCIATION parameters (inherits from association.* by default)
_pheno = config.get('phenotype_association', {})
PHENO_ASSOC_CONFIGS = parse_association_configs(_pheno.get('configs', []))
PHENO_GAPIT_CONFIGS, PHENO_OTHER_CONFIGS = split_configs_by_engine(PHENO_ASSOC_CONFIGS)
PHENO_METHOD_REGEX = '|'.join(PHENO_ASSOC_CONFIGS.keys()) if PHENO_ASSOC_CONFIGS else 'EMMAX'
PHENO_MISSING = _pheno.get('missing_strategy', 'DROP')
if PHENO_ASSOC_CONFIGS:
    check_in_list(PHENO_MISSING, ['MEAN', 'MEDIAN', 'DROP'], 'phenotype_association.missing_strategy')
# Inherit from association.* with optional override
PHENO_COMBINE_METHOD = 'All'  # pipeline always uses All; combine strategy moved to gradient_forest config
PHENO_COMBINE_GAP = _pheno.get('combine_gap', SIGSNPS_GAP)
PHENO_SNP_DISTANCE = PHENO_COMBINE_GAP  # overlap annotation distance
_pheno_rdist_raw = _pheno.get('region_distance', 1000000)
PHENO_REGION_DISTANCE_AUTO = (str(_pheno_rdist_raw).lower() == 'auto')
PHENO_REGION_DISTANCE = 'auto' if PHENO_REGION_DISTANCE_AUTO else int(_pheno_rdist_raw)
PHENO_PROMOTER_LENGTH = _pheno.get('promoter_length', PROMOTER_LENGTH)

# OVERLAP parameters (GEA + GWAS combined analysis)
_overlap = config.get('overlap', {})
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
OVERLAP_PROMOTER_LENGTH = _overlap.get('promoter_length', None)
if OVERLAP_PROMOTER_LENGTH is None:
    OVERLAP_PROMOTER_LENGTH = max(PROMOTER_LENGTH, PHENO_PROMOTER_LENGTH)

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
if PHENO_ASSOC_CONFIGS:
    PHENO_TRAITS = _meta_header[4:] if META_HAS_PHENO else []
    PHENO_PREDICTORS = ','.join(PHENO_TRAITS)

# HAPLOTYPE parameters
_hap = config.get('haplotype', {})
_hap_scan = _hap.get('scan', {})
HAP_SCAN_SOURCE = _hap_scan.get('regions_source', 'association')
HAP_SCAN_REGIONS_FILE = _hap_scan.get('regions_file', 'NULL')
HAP_SCAN_TOP_REGIONS = _hap_scan.get('top_regions', 5)
HAP_SCAN_MIN_SNPS = int(_hap_scan.get('min_snps', 10))
HAP_SCAN_MGMIN = _hap_scan.get('min_group_size', 50)
HAP_SCAN_MINHAP = _hap_scan.get('min_haplotype_size', 15)
HAP_SCAN_EPSILON_RANGE = _hap_scan.get('epsilon_range', [0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0])
HAP_SCAN_META_TYPE = _hap_scan.get('metadata_type', 'site')
HAP_EPSILON_SELECTED = _hap.get('epsilon_selected', 'NULL')

# Expand metadata types: 'both' -> ['site', 'cluster_K{K_BEST}']
def get_hap_meta_types():
    if HAP_SCAN_META_TYPE == 'site':
        return ['site']
    elif HAP_SCAN_META_TYPE == 'cluster':
        return [f'cluster_K{K_BEST}'] if K_BEST else ['cluster']
    elif HAP_SCAN_META_TYPE == 'both':
        types = ['site']
        if K_BEST:
            types.append(f'cluster_K{K_BEST}')
        return types
    return ['site']

HAP_META_TYPES = get_hap_meta_types()

# Expand region sources: 'all' -> every available source
def get_hap_sources():
    if HAP_SCAN_SOURCE == 'all':
        sources = []
        if ASSOC_CONFIGS:
            sources.append('association')
        if PHENO_ASSOC_CONFIGS:
            sources.append('association_phenotypes')
        if ASSOC_CONFIGS and PHENO_ASSOC_CONFIGS:
            sources.append('overlapping')
        return sources if sources else ['association']
    # Validate source against available configs
    src = HAP_SCAN_SOURCE
    if src == 'association_phenotypes' and not PHENO_ASSOC_CONFIGS:
        raise ValueError(f"HAPLOTYPE_SCAN_REGIONS_SOURCE is '{src}' but PHENO_ASSOC_CONFIGS is not set. "
                         f"Run mode=association_phenotypes first or change source to 'association'.")
    if src == 'overlapping' and not (ASSOC_CONFIGS and PHENO_ASSOC_CONFIGS):
        raise ValueError(f"HAPLOTYPE_SCAN_REGIONS_SOURCE is '{src}' but both ASSOC_CONFIGS and PHENO_ASSOC_CONFIGS are required. "
                         f"Run mode=association and mode=association_phenotypes first.")
    return [src]

HAP_SOURCES = get_hap_sources()

# Build list of (meta_tag, source) combos as combined tags
HAP_TAGS = [f"{mt}_{src}" for mt in HAP_META_TYPES for src in HAP_SOURCES]

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
MOD_PROC = f"{OUTDIR}processing/"
MOD_STRUCT = f"{OUTDIR}structure/"
MOD_CLIMATE = f"{OUTDIR}climate/"
MOD_STRUCTK = f"{OUTDIR}structure_k/"
MOD_ASSOC = f"{OUTDIR}association/"
MOD_PHENO = f"{OUTDIR}phenotype_association/"
MOD_OVERLAP = f"{OUTDIR}overlapping/"
MOD_MALAD = f"{OUTDIR}maladaptation/"

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
    'metadata': f"{MOD_PROC}tables/metadata.tsv",
    'sample_missing_stats': f"{MOD_PROC}tables/sample_missing_stats.tsv",
    'pca': f"{MOD_STRUCT}plots/pca.png",
    'pca_svg': f"{MOD_STRUCT}plots/pca.svg",
    'tracy': f"{MOD_STRUCT}plots/tracy_widom.png",
    'cross_entropy': f"{MOD_STRUCT}plots/cross_entropy_K{K_START}-{K_END}.png",
    # QC intermediate stats
    'qc_raw_summary':      f"{MOD_PROC}tables/vcf_raw_summary.tsv",
    'qc_filtering_summary':f"{MOD_PROC}tables/filtering_summary.tsv",
    'qc_sample_het':       f"{MOD_PROC}tables/sample_heterozygosity.tsv",
    'qc_maf_raw':          f"{INTER}qc/maf_raw.frq",
    'qc_maf_filtered':     f"{INTER}qc/maf_filtered.frq",
    'qc_snp_miss_raw':     f"{INTER}qc/snp_missingness_raw.lmiss",
    'qc_snp_density_raw':  f"{INTER}qc/snp_density_raw.snpden",
    'qc_snp_density_filt': f"{INTER}qc/snp_density_filtered.snpden",
    'qc_depth_sample':     f"{INTER}qc/depth_per_sample.idepth",
    'qc_depth_site':       f"{INTER}qc/depth_per_site.ldepth.mean",
    # QC plots
    'qc_plot_sample_miss': f"{MOD_PROC}plots/sample_missingness_distribution.png",
    'qc_plot_het_miss':    f"{MOD_PROC}plots/het_vs_missingness.png",
    'qc_plot_maf':         f"{MOD_PROC}plots/maf_distribution.png",
    'qc_plot_snp_miss':    f"{MOD_PROC}plots/snp_missingness_distribution.png",
    'qc_plot_attrition':   f"{MOD_PROC}plots/filtering_attrition.png",
    'qc_plot_snp_density': f"{MOD_PROC}plots/snp_density_by_chr.png",
    'qc_plot_depth':       f"{MOD_PROC}plots/depth_distribution.png",
}

# K_BEST dependent paths (added dynamically when K_BEST is set)
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
    O['tajima'] = f"{MOD_STRUCTK}tables/pop_stats/tajima_d_by_pop.tsv"
    O['pi_div'] = f"{MOD_STRUCTK}tables/pop_stats/pi_diversity_by_pop.tsv"
    O['ibd_raw'] = f"{MOD_STRUCTK}tables/pop_stats/ibd_raw.tsv"
    O['ibd_pairs'] = f"{MOD_STRUCTK}tables/pop_stats/ibd_pairs.tsv"
    O['amova'] = f"{MOD_STRUCTK}tables/pop_stats/amova.tsv"

    # Plots - climate
    O['corr_heatmap'] = f"{MOD_CLIMATE}plots/correlation_heatmap.png"
    # Plots - structure_k
    O['mantel'] = f"{MOD_STRUCTK}plots/pop_stats/mantel_test.png"
    O['amova_plot'] = f"{MOD_STRUCTK}plots/pop_stats/amova.png"

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
    O['ld_decay_table'] = f"{MOD_STRUCTK}tables/ld_decay_half_distances.tsv"
    O['ld_decay_plot_gw'] = f"{MOD_STRUCTK}plots/ld_decay/ld_decay_genome_wide.png"
    O['ld_decay_plot_gw_svg'] = f"{MOD_STRUCTK}plots/ld_decay/ld_decay_genome_wide.svg"
    O['ld_decay_plot_chr'] = f"{MOD_STRUCTK}plots/ld_decay/ld_decay_per_chromosome.png"
    O['ld_decay_plot_chr_svg'] = f"{MOD_STRUCTK}plots/ld_decay/ld_decay_per_chromosome.svg"

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

    # EMMAX work directory and pre-computed TPED/kinship
    W['emmax_work'] = f"{WORK_FILT}emmax/"
    W['assoc_tped'] = f"{WORK_FILT}emmax/{VCF_BASE}.tped"
    W['assoc_tfam'] = f"{WORK_FILT}emmax/{VCF_BASE}.tfam"
    W['assoc_kinship'] = f"{WORK_FILT}emmax/{VCF_BASE}.aBN.kinf"

    # GAPIT numeric format (GD + GM)
    W['gapit_gd'] = f"{WORK_FILT}gapit/{VCF_BASE}_GD.tsv"
    W['gapit_gm'] = f"{WORK_FILT}gapit/{VCF_BASE}_GM.tsv"
    W['gapit_work'] = f"{INTER}gapit/association/"

    # Combined outputs - association
    O['selected_snps'] = f"{MOD_ASSOC}tables/selected_snps.tsv"
    O['regions_per_trait'] = f"{MOD_ASSOC}tables/regions_per_trait.tsv"
    O['regions_combined'] = f"{MOD_ASSOC}tables/regions_combined.tsv"
    O['genes_per_region'] = f"{MOD_ASSOC}tables/genes_per_region.tsv"
    O['genes_per_region_collapsed'] = f"{MOD_ASSOC}tables/genes_per_region_collapsed.tsv"
    O['genes_combined_regions'] = f"{MOD_ASSOC}tables/genes_combined.tsv"
    W['enrichment_done'] = f"{INTER}flags/enrichment_done.flag"

    # Combined Manhattan plots (all traits and methods)
    O['manhattan_combined'] = f"{MOD_ASSOC}plots/manhattan/combined/manhattan_combined_K{K_BEST}.png"
    O['manhattan_combined_regions'] = f"{MOD_ASSOC}plots/manhattan/combined/manhattan_combined_K{K_BEST}_regions.png"
    O['qq_combined'] = f"{MOD_ASSOC}plots/manhattan/combined/qq_combined_K{K_BEST}.png"

    # Regionplot outputs
    O['gff_topr'] = f"{INTER}annotation/topr_gene_annotation.tsv"

add_association_paths()

# Phenotype association paths
def add_pheno_association_paths():
    """Add phenotype association paths to W and O dictionaries."""
    if K_BEST is None or not PHENO_ASSOC_CONFIGS:
        return

    PHENO_WORK = f"{WORK_FILT}phenotypes/"

    # Working paths
    W['pheno_work'] = PHENO_WORK
    W['pheno_emmax_work'] = f"{PHENO_WORK}emmax/"
    W['pheno_gapit_work'] = f"{INTER}gapit/phenotype_association/"

    # Ensure GAPIT paths exist even if association mode isn't configured
    if PHENO_GAPIT_CONFIGS and 'gapit_gd' not in W:
        W['gapit_gd'] = f"{WORK_FILT}gapit/{VCF_BASE}_GD.tsv"
        W['gapit_gm'] = f"{WORK_FILT}gapit/{VCF_BASE}_GM.tsv"
        W['gapit_work'] = f"{INTER}gapit/association/"
    if PHENO_GAPIT_CONFIGS and 'assoc_kinship' not in W:
        W['emmax_work'] = f"{WORK_FILT}emmax/"
        W['assoc_tped'] = f"{WORK_FILT}emmax/{VCF_BASE}.tped"
        W['assoc_tfam'] = f"{WORK_FILT}emmax/{VCF_BASE}.tfam"
        W['assoc_kinship'] = f"{WORK_FILT}emmax/{VCF_BASE}.aBN.kinf"

    # Path A (MEAN/MEDIAN) — single set of TPED/kinship
    W['pheno_tped'] = f"{PHENO_WORK}emmax/{VCF_BASE}.tped"
    W['pheno_tfam'] = f"{PHENO_WORK}emmax/{VCF_BASE}.tfam"
    W['pheno_kinship'] = f"{PHENO_WORK}emmax/{VCF_BASE}.aBN.kinf"
    W['pheno_all_phenotypes'] = f"{PHENO_WORK}all_phenotypes.tsv"

    # Prep outputs
    O['pheno_missing_summary'] = f"{MOD_PHENO}tables/phenotype_missing_summary.tsv"

    # Combined analysis outputs
    O['pheno_selected_snps'] = f"{MOD_PHENO}tables/selected_snps.tsv"
    O['pheno_regions_per_trait'] = f"{MOD_PHENO}tables/regions_per_trait.tsv"
    O['pheno_regions_combined'] = f"{MOD_PHENO}tables/regions_combined.tsv"
    O['pheno_genes_per_region'] = f"{MOD_PHENO}tables/genes_per_region.tsv"
    O['pheno_genes_collapsed'] = f"{MOD_PHENO}tables/genes_per_region_collapsed.tsv"
    O['pheno_genes_combined'] = f"{MOD_PHENO}tables/genes_combined.tsv"
    W['pheno_enrichment_done'] = f"{INTER}flags/pheno_enrichment_done.flag"

    # Manhattan combined
    O['pheno_manhattan_combined'] = f"{MOD_PHENO}plots/manhattan/combined/manhattan_combined_K{K_BEST}.png"
    O['pheno_manhattan_combined_regions'] = f"{MOD_PHENO}plots/manhattan/combined/manhattan_combined_K{K_BEST}_regions.png"
    O['pheno_qq_combined'] = f"{MOD_PHENO}plots/manhattan/combined/qq_combined_K{K_BEST}.png"

add_pheno_association_paths()

# Overlap paths (GEA + GWAS combined analysis)
def add_overlap_paths():
    """Add overlap analysis paths to W and O dictionaries."""
    if K_BEST is None:
        return

    has_gea  = bool(ASSOC_CONFIGS)
    has_gwas = bool(PHENO_ASSOC_CONFIGS)

    # Pairwise trait overlap — available when at least one source exists
    if has_gea or has_gwas:
        O['pairwise_collapsed_snps'] = f"{MOD_OVERLAP}tables/pairwise_collapsed_snps.tsv"
        O['pairwise_overlap_table']  = f"{MOD_OVERLAP}tables/pairwise_overlap_table.tsv"

    # GEA × GWAS combined analysis — requires both sources
    if has_gea and has_gwas:
        O['overlap_selected_snps'] = f"{MOD_OVERLAP}tables/selected_snps_all.tsv"
        O['overlap_regions_per_trait'] = f"{MOD_OVERLAP}tables/regions_per_trait_all.tsv"
        O['overlap_regions_combined'] = f"{MOD_OVERLAP}tables/regions_combined.tsv"
        O['overlap_summary'] = f"{MOD_OVERLAP}tables/overlap_summary.tsv"
        O['overlap_genes_per_region'] = f"{MOD_OVERLAP}tables/genes_per_region.tsv"
        O['overlap_genes_collapsed'] = f"{MOD_OVERLAP}tables/genes_per_region_collapsed.tsv"
        O['overlap_genes_combined'] = f"{MOD_OVERLAP}tables/genes_combined.tsv"
        W['overlap_enrichment_done'] = f"{INTER}flags/overlap_enrichment_done.flag"
        # Miami plots (static GEA vs GWAS)
        O['overlap_miami'] = f"{MOD_OVERLAP}plots/miami_combined_K{K_BEST}.png"
        O['overlap_miami_svg'] = f"{MOD_OVERLAP}plots/miami_combined_K{K_BEST}.svg"
        O['overlap_miami_regions'] = f"{MOD_OVERLAP}plots/miami_combined_K{K_BEST}_regions.png"
        O['overlap_miami_regions_svg'] = f"{MOD_OVERLAP}plots/miami_combined_K{K_BEST}_regions.svg"

add_overlap_paths()

# Haplotype analysis paths
def add_haplotype_paths():
    """Add haplotype analysis paths to W and O dictionaries."""
    if K_BEST is None:
        return
    for tag in HAP_TAGS:
        W[f'hap_scan_done_{tag}'] = f"{INTER}flags/haplotype_{tag}_scan_done.flag"
        W[f'hap_viz_done_{tag}'] = f"{INTER}flags/haplotype_{tag}_viz_done.flag"
        O[f'hap_selected_regions_{tag}'] = f"{OUTDIR}haplotype_scan/{tag}/tables/selected_regions.tsv"

add_haplotype_paths()

# Maladaptation paths
def add_maladaptation_paths():
    """Add maladaptation-specific paths to W and O dictionaries."""
    if K_BEST is None or not MODELS_LIST:
        return

    SUFFIX = f"{GF_SUFFIX}_{PCNM_TAG}"

    # Per-model future climate rasters
    for model in MODELS_LIST:
        W[f'climate_future_{model}'] = f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_{model}.tif"

    # Merged future climate
    W['climate_future_raster'] = f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_rasterstack.tif"
    O['climate_future_all'] = f"{MOD_CLIMATE}tables/future/climate_future_year{YEAR}_ssp{SSP}_all.tsv"
    O['climate_future_site'] = f"{MOD_CLIMATE}tables/future/climate_future_year{YEAR}_ssp{SSP}_site.tsv"

    # Gradient Forest models
    W['gf_selected_snps'] = f"{INTER}gradient_forest/gf_selected_snps_{SUFFIX}.tsv"
    W['gf_adaptive'] = f"{INTER}gradient_forest/gradientforest_adaptive_{SUFFIX}.qs"
    W['gf_random'] = f"{INTER}gradient_forest/gradientforest_random_{SUFFIX}.qs"

    # Genetic offset outputs
    W['gf_offset_raster'] = f"{INTER}gradient_forest/genetic_offset_{SUFFIX}.tif"
    O['gf_offset_map_values'] = f"{MOD_MALAD}tables/{SUFFIX}/genetic_offset_map.tsv"
    O['gf_offset_site_values'] = f"{MOD_MALAD}tables/{SUFFIX}/genetic_offset_site.tsv"

    # Gradient Forest plots
    O['gf_offset_piemap'] = f"{MOD_MALAD}plots/{SUFFIX}/genetic_offset_piemap.png"
    O['gf_offset_piemap_tajima'] = f"{MOD_MALAD}plots/{SUFFIX}/genetic_offset_piemap_tajima_d.png"
    O['gf_offset_piemap_diversity'] = f"{MOD_MALAD}plots/{SUFFIX}/genetic_offset_piemap_pi_diversity.png"
    O['gf_cumimp'] = f"{MOD_MALAD}plots/{SUFFIX}/cumulative_importance.png"
    O['gf_importance'] = f"{MOD_MALAD}plots/{SUFFIX}/overall_importance.png"

    # Future climate density plot
    O['density_future'] = f"{MOD_CLIMATE}plots/density_plot_future_ssp{SSP}_{YEAR}.png"

    # Create suffix-based directories
    os.makedirs(f"{MOD_MALAD}plots/{SUFFIX}/", exist_ok=True)
    os.makedirs(f"{MOD_MALAD}tables/{SUFFIX}/", exist_ok=True)

add_maladaptation_paths()

# Templates for K-dependent outputs
def clusters_table(k): return f"{MOD_STRUCT}tables/K{k}/clusters_K{k}.tsv"
def structure_plot(k): return f"{MOD_STRUCT}plots/K{k}/structure_K{k}.png"
def pca_struct_plot(k): return f"{MOD_STRUCT}plots/K{k}/pca_structure_K{k}.png"
def pop_diff_plot(k): return f"{MOD_STRUCT}plots/K{k}/pop_diff_K{k}.png"

# Templates for climate/trait-dependent outputs
DENSITY_PLOT_COMBINED    = f"{MOD_CLIMATE}plots/density_plot_present.png"
DENSITY_PLOT_PHENOTYPES  = f"{MOD_CLIMATE}plots/density_plot_phenotypes.png"
def piemap_tajima(bio): return f"{MOD_STRUCTK}plots/piemap/piemap_{bio}_tajima_d.png"
def piemap_diversity(bio): return f"{MOD_STRUCTK}plots/piemap/piemap_{bio}_pi_diversity.png"
def piemap_notrait(bio): return f"{MOD_STRUCTK}plots/piemap/piemap_{bio}.png"

# Templates for association outputs
def assoc_pvalues(method): return f"{MOD_ASSOC}tables/methods/{method}/{method}_pvalues_K{K_BEST}.tsv"
def assoc_sigsnps(method, adjust): return f"{MOD_ASSOC}tables/methods/{method}/{method}_pvalues_K{K_BEST}_sig_snps_{adjust}.tsv"
def manhattan_plot(method, trait, adjust): return f"{MOD_ASSOC}plots/manhattan/{method}/manhattan_{trait}_K{K_BEST}_{adjust}.png"
def manhattan_plot_regions(method, trait, adjust): return f"{MOD_ASSOC}plots/manhattan/{method}/manhattan_{trait}_K{K_BEST}_{adjust}_regions.png"
def qq_plot(method, trait, adjust): return f"{MOD_ASSOC}plots/manhattan/{method}/qq_{trait}_K{K_BEST}_{adjust}.png"

# Templates for phenotype association outputs
def pheno_pvalues(method): return f"{MOD_PHENO}tables/methods/{method}/{method}_pvalues_K{K_BEST}.tsv"
def pheno_qvalues(method): return f"{MOD_PHENO}tables/methods/{method}/{method}_qvalues_K{K_BEST}.tsv"
def pheno_sigsnps(method, adjust): return f"{MOD_PHENO}tables/methods/{method}/{method}_pvalues_K{K_BEST}_sig_snps_{adjust}.tsv"
def pheno_manhattan(method, trait, adjust): return f"{MOD_PHENO}plots/manhattan/{method}/manhattan_{trait}_K{K_BEST}_{adjust}.png"
def pheno_manhattan_regions(method, trait, adjust): return f"{MOD_PHENO}plots/manhattan/{method}/manhattan_{trait}_K{K_BEST}_{adjust}_regions.png"
def pheno_qq(method, trait, adjust): return f"{MOD_PHENO}plots/manhattan/{method}/qq_{trait}_K{K_BEST}_{adjust}.png"

# Per-trait DROP mode paths
def pheno_trait_vcf(trait): return f"{WORK_FILT}phenotypes/{trait}/{VCF_BASE}.vcf"
def pheno_trait_tped(trait): return f"{WORK_FILT}phenotypes/{trait}/emmax/{VCF_BASE}.tped"
def pheno_trait_tfam(trait): return f"{WORK_FILT}phenotypes/{trait}/emmax/{VCF_BASE}.tfam"
def pheno_trait_kinship(trait): return f"{WORK_FILT}phenotypes/{trait}/emmax/{VCF_BASE}.aBN.kinf"
def pheno_trait_pvalues(method, trait): return f"{MOD_PHENO}tables/methods/{method}/{trait}_pvalues_K{K_BEST}.tsv"

#=============================================================================
# CREATE DIRECTORIES
#=============================================================================
dirs_to_create = [
    WORK, WORK_FILT, WORK_LD, INTER,
    f"{INTER}samples/", f"{INTER}annotation/", f"{INTER}flags/", f"{INTER}qc/",
    LOGDIR,
    # Module directories
    f"{MOD_PROC}tables/", f"{MOD_PROC}plots/",
    f"{MOD_STRUCT}plots/", f"{MOD_STRUCT}tables/",
    *[f"{MOD_STRUCT}plots/K{k}/" for k in k_range(K_START, K_END)],
    *[f"{MOD_STRUCT}tables/K{k}/" for k in k_range(K_START, K_END)],
    *([ f"{MOD_CLIMATE}plots/", f"{MOD_CLIMATE}tables/present/", f"{MOD_CLIMATE}rasters/present/" ] if CLIMATE_ENABLED else []),
    f"{MOD_STRUCTK}plots/piemap/", f"{MOD_STRUCTK}tables/",
    # Log subdirectories
    f"{LOGDIR}processing/", f"{LOGDIR}structure/", f"{LOGDIR}structure_k/",
    f"{LOGDIR}association/", f"{LOGDIR}maladaptation/",
]

# Add association directories for each method
for method in ASSOC_CONFIGS:
    dirs_to_create.append(f"{MOD_ASSOC}plots/manhattan/{method}/")
    dirs_to_create.append(f"{MOD_ASSOC}tables/methods/{method}/")

# Add enrichment directories
dirs_to_create.append(f"{MOD_ASSOC}tables/")
dirs_to_create.append(f"{MOD_ASSOC}tables/enrichment/")
dirs_to_create.append(f"{MOD_ASSOC}plots/")
dirs_to_create.append(f"{MOD_ASSOC}plots/enrichment/")
dirs_to_create.append(f"{MOD_ASSOC}plots/manhattan/combined/")
dirs_to_create.append(f"{INTER}enrichment/association/")

# Add GAPIT directories
if ASSOC_GAPIT_CONFIGS:
    dirs_to_create.append(f"{WORK_FILT}gapit/")
    dirs_to_create.append(f"{INTER}gapit/association/")
    dirs_to_create.append(f"{MOD_ASSOC}GAPIT_native_output/")
    for model in ASSOC_GAPIT_CONFIGS:
        dirs_to_create.append(f"{MOD_ASSOC}GAPIT_native_output/{model}/")

# Add maladaptation directories (require climate)
if CLIMATE_ENABLED:
    dirs_to_create.append(f"{INTER}gradient_forest/")
    dirs_to_create.append(f"{MOD_CLIMATE}rasters/future/")
    dirs_to_create.append(f"{MOD_CLIMATE}tables/future/")

# Add phenotype association directories
if PHENO_ASSOC_CONFIGS:
    dirs_to_create.append(f"{MOD_PHENO}tables/")
    dirs_to_create.append(f"{MOD_PHENO}tables/enrichment/")
    dirs_to_create.append(f"{MOD_PHENO}plots/")
    dirs_to_create.append(f"{MOD_PHENO}plots/enrichment/")
    dirs_to_create.append(f"{MOD_PHENO}plots/piemap/")
    dirs_to_create.append(f"{MOD_PHENO}plots/manhattan/combined/")
    dirs_to_create.append(f"{INTER}enrichment_phenotypes/")
    dirs_to_create.append(f"{LOGDIR}phenotype_association/")
    for method in PHENO_ASSOC_CONFIGS:
        dirs_to_create.append(f"{MOD_PHENO}plots/manhattan/{method}/")
        dirs_to_create.append(f"{MOD_PHENO}tables/methods/{method}/")
    if PHENO_MISSING == 'DROP':
        for trait in PHENO_TRAITS:
            dirs_to_create.append(f"{WORK_FILT}phenotypes/{trait}/emmax/")
    if PHENO_GAPIT_CONFIGS:
        dirs_to_create.append(f"{INTER}gapit/phenotype_association/")
        dirs_to_create.append(f"{MOD_PHENO}GAPIT_native_output/")
        for model in PHENO_GAPIT_CONFIGS:
            dirs_to_create.append(f"{MOD_PHENO}GAPIT_native_output/{model}/")
        if not ASSOC_GAPIT_CONFIGS:
            # GAPIT numeric dir needed even when only pheno uses GAPIT
            dirs_to_create.append(f"{WORK_FILT}gapit/")

# Add overlapping analysis directories
if ASSOC_CONFIGS and PHENO_ASSOC_CONFIGS:
    dirs_to_create.append(f"{MOD_OVERLAP}tables/")
    dirs_to_create.append(f"{MOD_OVERLAP}tables/enrichment/")
    dirs_to_create.append(f"{MOD_OVERLAP}plots/")
    dirs_to_create.append(f"{MOD_OVERLAP}plots/enrichment/")
    dirs_to_create.append(f"{INTER}enrichment_overlapping/")
    dirs_to_create.append(f"{LOGDIR}overlapping/")

# Add haplotype analysis directories
for tag in HAP_TAGS:
    dirs_to_create.append(f"{OUTDIR}haplotype_scan/{tag}/plots/clustree/")
    dirs_to_create.append(f"{OUTDIR}haplotype_scan/{tag}/tables/")
    dirs_to_create.append(f"{OUTDIR}haplotype/{tag}/plots/")
    dirs_to_create.append(f"{OUTDIR}haplotype/{tag}/tables/")
    dirs_to_create.append(f"{INTER}haplotype/{tag}/")
dirs_to_create.append(f"{LOGDIR}haplotype_scan/")
dirs_to_create.append(f"{LOGDIR}haplotype/")

# Add LD decay directories
if K_BEST is not None:
    dirs_to_create.append(W['ld_decay_work'])
    dirs_to_create.append(W['ld_decay_sample_lists'])
    dirs_to_create.append(W['ld_decay_chr_vcfs'])
    dirs_to_create.append(W['ld_decay_stat_gw'])
    dirs_to_create.append(W['ld_decay_stat_chr'])

# Add pop stats directories if enabled
if CALC_POP_STATS:
    dirs_to_create.append(f"{MOD_STRUCTK}plots/pop_stats/")
    dirs_to_create.append(f"{MOD_STRUCTK}tables/pop_stats/")

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
    
    elif mode == 'structure':
        ks = k_range(K_START, K_END)
        return (
            [clusters_table(k) for k in ks] +
            [structure_plot(k) for k in ks] +
            [pca_struct_plot(k) for k in ks] +
            [pop_diff_plot(k) for k in ks] +
            [O['pca'], O['tracy'], O['cross_entropy'], W['summary_done']]
        )
    
    elif mode == 'structure_K':
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

    elif mode == 'association':
        if not CLIMATE_ENABLED:
            raise ValueError("association mode requires climate.enabled: true (GEA uses climate predictors as traits)")
        check_numeric(K_BEST, 'K_BEST')
        predictors = get_predictors_list()
        if not predictors:
            raise ValueError("climate.predictors must be set for association mode")
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
                targets.append(qq_plot(method, trait, adjust))

        # Combined analysis targets
        targets.append(O['selected_snps'])
        targets.append(O['regions_per_trait'])
        targets.append(O['regions_combined'])
        targets.append(O['genes_per_region'])
        targets.append(O['genes_per_region_collapsed'])
        targets.append(O['genes_combined_regions'])

        # Combined Manhattan plots (all traits and methods)
        targets.append(O['manhattan_combined'])
        targets.append(O['manhattan_combined_regions'])
        targets.append(O['qq_combined'])

        # Enrichment (if GO_FIELD is specified)
        if GO_FIELD and GO_FIELD != 'NULL':
            targets.append(W['enrichment_done'])

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

        targets.append(W['summary_done'])
        return targets

    elif mode in ('association_phenotypes', 'phenotype_association'):
        check_numeric(K_BEST, 'K_BEST')
        if not PHENO_TRAITS:
            raise ValueError("No phenotype columns found in metadata (expected columns 5+ after site, sample, latitude, longitude)")
        if not PHENO_ASSOC_CONFIGS:
            raise ValueError("PHENO_ASSOC_CONFIGS must be set for association_phenotypes mode")
        if GFF:
            check_file_exists(INDIR, GFF, 'GFF')

        targets = [O['pheno_missing_summary']]

        # Per-method targets
        for method, adjust in PHENO_ASSOC_CONFIGS.items():
            targets.append(pheno_pvalues(method))
            targets.append(pheno_sigsnps(method, adjust))
            for trait in PHENO_TRAITS:
                targets.append(pheno_manhattan(method, trait, adjust))
                targets.append(pheno_manhattan_regions(method, trait, adjust))
                targets.append(pheno_qq(method, trait, adjust))

        targets.extend([
            O['pheno_selected_snps'],
            O['pheno_regions_per_trait'], O['pheno_regions_combined'],
            O['pheno_genes_per_region'], O['pheno_genes_collapsed'], O['pheno_genes_combined'],
            O['pheno_manhattan_combined'], O['pheno_manhattan_combined_regions'],
            O['pheno_qq_combined'],
        ])

        # Phenotype piemaps (require climate raster for background)
        if CLIMATE_ENABLED:
            for trait in PHENO_TRAITS:
                targets.append(f"{MOD_PHENO}plots/piemap/phenomap_{trait}.png")

        # Enrichment (if GO_FIELD is specified)
        if GO_FIELD and GO_FIELD != 'NULL':
            targets.append(W['pheno_enrichment_done'])

        targets.append(W['summary_done'])
        return targets

    elif mode == 'overlapping':
        check_numeric(K_BEST, 'K_BEST')
        has_gea  = bool(ASSOC_CONFIGS) and CLIMATE_ENABLED
        has_gwas = bool(PHENO_ASSOC_CONFIGS)
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

        # GEA × GWAS combined analysis (only when both sources available)
        if has_gea and has_gwas:
            targets.extend([
                O['overlap_selected_snps'],
                O['overlap_regions_per_trait'],
                O['overlap_regions_combined'],
                O['overlap_summary'],
                O['overlap_genes_per_region'],
                O['overlap_genes_collapsed'],
                O['overlap_genes_combined'],
            ])
            targets.extend([
                O['overlap_miami'], O['overlap_miami_svg'],
                O['overlap_miami_regions'], O['overlap_miami_regions_svg'],
            ])
            if GO_FIELD and GO_FIELD != 'NULL':
                targets.append(W['overlap_enrichment_done'])

        targets.append(W['summary_done'])
        return targets

    elif mode == 'haplotype_scan':
        check_numeric(K_BEST, 'K_BEST')
        targets = []
        for tag in HAP_TAGS:
            targets.append(O[f'hap_selected_regions_{tag}'])
            targets.append(W[f'hap_scan_done_{tag}'])
        targets.append(W['summary_done'])
        return targets

    elif mode == 'haplotype':
        check_numeric(K_BEST, 'K_BEST')
        if str(HAP_EPSILON_SELECTED) == 'NULL':
            raise ValueError("Set HAPLOTYPE_EPSILON_SELECTED after reviewing clustree plots from haplotype_scan mode")
        targets = []
        for tag in HAP_TAGS:
            targets.append(W[f'hap_viz_done_{tag}'])
        targets.append(W['summary_done'])
        return targets

    elif mode is None:
        raise ValueError("Specify mode: --config mode=processing or mode=structure or mode=structure_K or mode=association or mode=association_phenotypes or mode=phenotype_association or mode=overlapping or mode=haplotype_scan or mode=haplotype")
    else:
        raise ValueError(f"Unknown mode: {mode}")

# Helper: provide LD decay table as input when region_distance is "auto"
def ld_decay_input(auto_flag):
    """Return LD decay table path as input dependency when auto mode is enabled."""
    return O.get('ld_decay_table', []) if auto_flag else []

#=============================================================================
# MAIN RULE
