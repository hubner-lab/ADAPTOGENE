# ADAPTOGENE Pipeline - Refactored
# vim: filetype=python
import os
import sys
import subprocess
from pathlib import Path

# Method registries — import before config parsing so GAPIT_MODELS can be derived
sys.path.insert(0, os.path.join(workflow.basedir, "workflow"))
from methods.gea import GEA_METHODS, K_BEST_SENTINEL
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
PI_HAT = _cfg('Filter', 'relatedness', None); check_float(PI_HAT, 'filter.relatedness', allow_null=True)
# relatedness_action gates removal separately from the threshold: 'keep' (default) only
# colors/counts related pairs in the Processing histogram; 'remove' actually drops samples.
_relatedness_action_explicit = 'relatedness_action' in config.get('Filter', {})
RELATEDNESS_ACTION = str(_cfg('Filter', 'relatedness_action', 'keep')).lower()
check_in_list(RELATEDNESS_ACTION, ['keep', 'remove'], 'filter.relatedness_action')
if PI_HAT is not None and not _relatedness_action_explicit:
    print(
        f"WARNING: Filter.relatedness={PI_HAT} is set but Filter.relatedness_action is not — "
        "defaulting to 'keep' (related samples are visualized/counted but NOT removed). "
        "Set Filter.relatedness_action: remove to drop them.", file=sys.stderr)
RELATEDNESS_REMOVE = (PI_HAT is not None and RELATEDNESS_ACTION == 'remove')

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
CLIMATE_SOURCE = _cfg('Climate', 'source', 'worldclim')
check_in_list(CLIMATE_SOURCE, ['worldclim', 'custom'], 'Climate.source')
_clim_custom = config.get('Climate', {}).get('custom', {}) if CLIMATE_SOURCE == 'custom' else {}
CUSTOM_PRESENT_TABLE   = f"{INDIR}{_clim_custom.get('present_table', '')}"  if CLIMATE_SOURCE == 'custom' else ''
CUSTOM_FUTURE_TABLE    = f"{INDIR}{_clim_custom.get('future_table',  '')}"  if CLIMATE_SOURCE == 'custom' else ''
CUSTOM_CLIMATE_COLUMNS = _clim_custom.get('columns', PREDICTORS_SELECTED)
CUSTOM_GRID_RES        = _clim_custom.get('grid_resolution', 0.1)
CUSTOM_CLIMATE_KEY     = _clim_custom.get('key', 'site')

# Bioclimatic variable universe.
# WorldClim: all 19 variables (raster always ships bio_1..bio_19).
# Custom:    only the columns declared in Climate.custom.columns — narrows piemap/mantel fan-out
#            so no unbuildable targets are generated for layers the custom raster lacks.
if CLIMATE_SOURCE == 'custom':
    ALL_BIO = [c.strip() for c in str(CUSTOM_CLIMATE_COLUMNS).split(',') if c.strip()]
else:
    ALL_BIO = [f"bio_{i}" for i in range(1, 20)]
ALL_BIO_STR = ",".join(ALL_BIO)

# POP parameters
CALC_POP_STATS = _cfg_bool('Population', 'calc_stats', False)
METRICS_WINSIZE = _cfg('Population', 'window_size', 10000)
CUSTOM_TRAIT = _cfg('Population', 'custom_trait_file', 'NULL')

# PIEMAP parameters (palette fixed to viridis plasma)
# NOTE: use_points is NOT read here for the Structure/Maladaptation/GWAS geo
# maps — those always emit a pie file plus a "*_points" companion file, and
# the Shiny app toggles between them at runtime. Piemap.use_points still
# governs the on-demand haplotype-viz piemaps (see fct_regions.R).
PIEMAP_ALPHA = _cfg('Piemap', 'alpha', 0.6)
PIEMAP_SHOW_LABELS = 'T' if _cfg_bool('Piemap', 'show_labels', False) else 'F'
PIEMAP_LABEL_SIZE = _cfg('Piemap', 'label_size', 10)
PIEMAP_PIE_SCALE = _cfg('Piemap', 'pie_scale', 1.0)

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

_VALID_REGION_MODES = ('auto_per_chromosome', 'auto_genome_wide', 'auto')

def _migrate_assoc_config(group_dict):
    """Upgrade legacy region param keys to current names with deprecation warnings."""
    import warnings
    d = dict(group_dict)
    if 'region_distance' in d and 'snp_clumping_distance' not in d:
        warnings.warn(
            "config key 'region_distance' is deprecated — rename to 'snp_clumping_distance'",
            DeprecationWarning, stacklevel=2
        )
        d['snp_clumping_distance'] = d.pop('region_distance')
    if 'combine_gap' in d and 'snp_clumping_distance' not in d:
        warnings.warn(
            "config key 'combine_gap' is deprecated — absorbed into 'snp_clumping_distance'",
            DeprecationWarning, stacklevel=2
        )
        d.setdefault('snp_clumping_distance', d.pop('combine_gap'))
    elif 'combine_gap' in d:
        d.pop('combine_gap')  # silently drop; already superseded
    if 'region_r2_threshold' in d and 'clumping_r2_threshold' not in d:
        warnings.warn(
            "config key 'region_r2_threshold' is deprecated — rename to 'clumping_r2_threshold'",
            DeprecationWarning, stacklevel=2
        )
        d['clumping_r2_threshold'] = d.pop('region_r2_threshold')
    if 'region_ld_decay_group' in d and 'ld_decay_group' not in d:
        warnings.warn(
            "config key 'region_ld_decay_group' is deprecated — rename to 'ld_decay_group'",
            DeprecationWarning, stacklevel=2
        )
        d['ld_decay_group'] = d.pop('region_ld_decay_group')
    return d

def resolve_region_params(group_dict, defaults=None):
    """Parse snp_clumping_distance, clumping_r2_threshold, ld_decay_group, promoter_length."""
    defaults = defaults or {}
    group_dict = _migrate_assoc_config(group_dict)
    defaults   = _migrate_assoc_config(defaults)
    raw_rdist = group_dict.get('snp_clumping_distance', defaults.get('snp_clumping_distance', 100000))
    raw_str   = str(raw_rdist).lower().strip()
    if raw_str in _VALID_REGION_MODES:
        rdist_mode = 'auto_genome_wide' if raw_str == 'auto' else raw_str
        rdist_val  = raw_str  # pass string spec to R script
        if raw_str == 'auto':
            import warnings
            warnings.warn("snp_clumping_distance: 'auto' is deprecated; use 'auto_genome_wide' or 'auto_per_chromosome'")
    else:
        rdist_mode = 'fixed'
        rdist_val  = int(raw_rdist)
    return {
        'clumping_distance':      rdist_val,
        'clumping_distance_mode': rdist_mode,
        'clumping_r2_threshold':  float(group_dict.get('clumping_r2_threshold', defaults.get('clumping_r2_threshold', 0.2))),
        'ld_decay_group':         str(group_dict.get('ld_decay_group', defaults.get('ld_decay_group', 'All'))),
        'promoter_length':        int(group_dict.get('promoter_length',  defaults.get('promoter_length', 10000))),
    }

_gea_rdp = resolve_region_params(_assoc)
PROMOTER_LENGTH       = _gea_rdp['promoter_length']
CLUMPING_DISTANCE     = _gea_rdp['clumping_distance']
CLUMPING_DISTANCE_MODE = _gea_rdp['clumping_distance_mode']
CLUMPING_R2_THRESHOLD = _gea_rdp['clumping_r2_threshold']
CLUMPING_LD_DECAY_GROUP = _gea_rdp['ld_decay_group']
# Legacy aliases — keep so existing references in summary.smk etc. don't break
SIGSNPS_GAP          = CLUMPING_DISTANCE
SNP_DISTANCE         = CLUMPING_DISTANCE
REGION_DISTANCE      = CLUMPING_DISTANCE
REGION_DISTANCE_MODE = CLUMPING_DISTANCE_MODE
REGION_R2_THRESHOLD  = CLUMPING_R2_THRESHOLD
REGION_LD_DECAY_GROUP = CLUMPING_LD_DECAY_GROUP
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
_mala_cfg = config.get('Maladaptation', {})
_gf = _mala_cfg.get('methods', {}).get('gradient_forest', {})
NTREE = _gf.get('ntree', '500')
COR_THRESHOLD = _gf.get('cor_threshold', '0.5')
SPATIAL_CORRECTION = _gf.get('spatial_correction', 'with')
GF_RANDOM_MODEL = _gf.get('random_model', True)
# Legacy keys — ignored now (SNP sets are named by the user in the Shiny GEA tab)
_GF_LEGACY = {k: _gf.get(k) for k in ('run_label', 'combine_method', 'combine_gap') if k in _gf}
if _GF_LEGACY:
    import sys; print(
        f"WARNING: Maladaptation.methods.gradient_forest keys {list(_GF_LEGACY)} are no longer "
        "used. SNP sets are now saved from the GEA tab in the Shiny app. "
        "See Maladaptation.snp_sets in config.", file=sys.stderr)

if SPATIAL_CORRECTION not in ('with', 'without', 'both'):
    raise ValueError(
        f"Maladaptation.methods.gradient_forest.spatial_correction must be "
        f"'with', 'without', or 'both'; got {SPATIAL_CORRECTION!r}"
    )
if SPATIAL_CORRECTION == 'both':
    ACTIVE_SPATIAL_TAGS = ['spatial', 'nospatial']
elif SPATIAL_CORRECTION == 'without':
    ACTIVE_SPATIAL_TAGS = ['nospatial']
else:  # 'with'
    ACTIVE_SPATIAL_TAGS = ['spatial']

# GEOMETRIC OFFSET parameters
_go = _mala_cfg.get('methods', {}).get('geometric_offset', {})
GO_CANDIDATE_ONLY = _go.get('candidate_only', True)
_go_k = _go.get('k', '')
GO_K = str(_go_k) if _go_k != '' else ''   # '' → use K_BEST at rule time
GO_SCALE = str(_go.get('scale', True)).upper()  # 'TRUE'/'FALSE' for R

# Raw SNP-sets config (cheap parse — actual glob/validate happens only inside get_targets)
SNP_SETS_CFG = _mala_cfg.get('snp_sets', 'all')

# Validate and derive active maladaptation methods from config
_mala_methods_cfg = config.get('Maladaptation', {}).get('methods', {})
_unknown_mala = set(_mala_methods_cfg.keys()) - set(MALADAPTATION_METHODS.keys())
if _unknown_mala:
    raise ValueError(f"Unknown maladaptation methods in config: {_unknown_mala}. "
                     f"Registered: {list(MALADAPTATION_METHODS.keys())}")
ACTIVE_MALA_METHODS = [m for m in _mala_methods_cfg if m in MALADAPTATION_METHODS] or ['gradient_forest']

def mala_spatial_tags(method):
    """Return the list of spatial tags this method supports.

    geometric_offset is nospatial-only; gradient_forest follows ACTIVE_SPATIAL_TAGS
    (which reflects the Maladaptation.methods.gradient_forest.spatial_correction setting).
    """
    if MALADAPTATION_METHODS.get(method, {}).get('supports_spatial', True):
        return ACTIVE_SPATIAL_TAGS
    return ['nospatial']

def _resolve_param_default(default):
    """K_BEST_SENTINEL -> K_BEST; everything else passes through unchanged.
    This is the sole backward-compat mechanism for per-method hyperparameters:
    every existing config (no 'params:' key) resolves LFMM K / EMMAX #PCs /
    GAPIT PCA.total / RDA condition_pcs to today's K_BEST, byte-identical."""
    return K_BEST if default == K_BEST_SENTINEL else default

def _coerce_param(method, name, spec, raw):
    """Cast a user-supplied param value to spec['type'] and enforce min/max/choices.
    Raises with method+param named so a YAML typo is a parse-time error."""
    ptype = spec["type"]
    try:
        if ptype == "int":
            val = int(raw)
        elif ptype == "float":
            val = float(raw)
        elif ptype == "bool":
            val = raw if isinstance(raw, bool) else str(raw).strip().lower() not in ('false', 'f', '0', 'no', '')
        else:  # "str" | "enum"
            val = str(raw)
    except (TypeError, ValueError) as e:
        raise ValueError(f"GEA method '{method}' param '{name}': cannot cast {raw!r} to {ptype}: {e}")
    if ptype == "enum" and val not in spec.get("choices", []):
        raise ValueError(f"GEA method '{method}' param '{name}' must be one of "
                         f"{spec.get('choices')}, got: {val!r}")
    if ptype in ("int", "float"):
        if "min" in spec and val < spec["min"]:
            raise ValueError(f"GEA method '{method}' param '{name}'={val} is below min={spec['min']}")
        if "max" in spec and val > spec["max"]:
            raise ValueError(f"GEA method '{method}' param '{name}'={val} is above max={spec['max']}")
    return val

def resolve_method_params(method, registry, user_params):
    """Registry defaults (sentinel-resolved) <- user YAML overrides <- validate.
    Special-cases 'auto'/'@k_best'-style string sentinels for str-typed params
    (e.g. RDA's axes/predictor_set/fit_mode) by passing them through unchanged —
    only int/float/bool/enum params get type coercion; str params are free text
    that downstream scripts interpret themselves."""
    spec_dict = registry.get(method, {}).get("params", {})
    resolved = {name: _resolve_param_default(spec["default"]) for name, spec in spec_dict.items()}
    user_params = user_params or {}
    unknown = set(user_params) - set(spec_dict)
    if unknown:
        raise ValueError(f"GEA method '{method}': unknown param(s) {sorted(unknown)}. "
                         f"Valid params: {sorted(spec_dict.keys())}")
    for name, raw in user_params.items():
        spec = spec_dict[name]
        resolved[name] = raw if spec["type"] == "str" else _coerce_param(method, name, spec, raw)
    return resolved

def parse_association_configs(configs_list, registry=None):
    """Parse association configs list into method -> adjust_threshold dict.
    When `registry` is given (the GEA_METHODS/GWAS_METHODS dict), also returns
    a sibling method -> resolved-params dict as a second tuple element."""
    configs = {}
    params = {}
    for cfg in (configs_list or []):
        method = cfg['method']
        adjust = f"{cfg['adjust']}_{cfg['threshold']}"
        if method in configs:
            raise ValueError(f"Method '{method}' appears multiple times in configs")
        configs[method] = adjust
        if registry is not None:
            params[method] = resolve_method_params(method, registry, cfg.get('params', {}) or {})
    return configs, params

GEA_CONFIGS, GEA_PARAMS = parse_association_configs(_assoc.get('configs', []), GEA_METHODS)

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
GWAS_CONFIGS, GWAS_PARAMS = parse_association_configs(_pheno.get('configs', []), GWAS_METHODS)
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

def _validate_rda_semantics(configs, params):
    """RDA-specific config validation, run at parse time so a bad setup fails
    loudly before any expensive rule runs — not mid-run inside rda.R."""
    for m, adj in configs.items():
        if GEA_METHODS[m]["engine"] != "rda":
            continue
        p = params[m]
        rule_name = adj.split("_", 1)[0]
        if rule_name not in ("bonf", "qval"):
            print(
                f"WARNING: GEA method '{m}': significance rule '{rule_name}' is outside "
                "the RDA literature. Capblancq & Forester (2021) use Bonferroni 0.01/m; "
                "Capblancq et al. (2018) use q<0.1 (docs/rda_research.md A.4). "
                f"'{rule_name}' is a valid pipeline rule but has no published RDA "
                "calibration — the candidates table and diagnostics record it as "
                "candidate_rule so the deviation stays visible downstream.",
                file=sys.stderr)
        if not CLIMATE_ENABLED:
            raise ValueError(f"GEA method '{m}': RDA requires Climate.enabled: true")
        preds_raw = p["predictor_set"] if p["predictor_set"] != "auto" else PREDICTORS_SELECTED
        preds = [x.strip() for x in str(preds_raw).split(',') if x.strip()]
        if len(preds) < 2:
            raise ValueError(f"GEA method '{m}': RDA needs >=2 climate predictors, got {preds}")
        if p["axes"] != "auto":
            try:
                axes_int = int(p["axes"])
            except ValueError:
                raise ValueError(f"GEA method '{m}': axes must be 'auto' or an integer, got {p['axes']!r}")
            if axes_int < 2:
                raise ValueError(f"GEA method '{m}': axes must be >= 2 (covRob needs >= 2 "
                                 f"loading columns — verified empirically: K=1 crashes). Got {axes_int}.")
        if p["fit_mode"] == "pruned":
            print(
                f"WARNING: GEA method '{m}': fit_mode='pruned' is an ENGINEERING FALLBACK, "
                "not a scientific alternative. Capblancq & Forester (2021) explicitly chose "
                "NOT to pre-prune (docs/rda_research.md A6) — results deviate from the "
                "canonical method.", file=sys.stderr)

if GEA_CONFIGS:
    _validate_rda_semantics(GEA_CONFIGS, GEA_PARAMS)

# Inherit from GEA.* with optional override
PHENO_COMBINE_METHOD = 'All'  # pipeline always uses All; combine strategy moved to gradient_forest config
_gwas_rdp = resolve_region_params(_pheno, defaults=_assoc)
PHENO_CLUMPING_DISTANCE      = _gwas_rdp['clumping_distance']
PHENO_CLUMPING_DISTANCE_MODE = _gwas_rdp['clumping_distance_mode']
PHENO_CLUMPING_R2_THRESHOLD  = _gwas_rdp['clumping_r2_threshold']
PHENO_CLUMPING_LD_DECAY_GROUP = _gwas_rdp['ld_decay_group']
PHENO_PROMOTER_LENGTH        = _gwas_rdp['promoter_length']
# Legacy aliases
PHENO_COMBINE_GAP           = PHENO_CLUMPING_DISTANCE
PHENO_SNP_DISTANCE          = PHENO_CLUMPING_DISTANCE
PHENO_REGION_DISTANCE       = PHENO_CLUMPING_DISTANCE
PHENO_REGION_DISTANCE_MODE  = PHENO_CLUMPING_DISTANCE_MODE
PHENO_REGION_R2_THRESHOLD   = PHENO_CLUMPING_R2_THRESHOLD
PHENO_REGION_LD_DECAY_GROUP = PHENO_CLUMPING_LD_DECAY_GROUP

# OVERLAP parameters (GEA + GWAS combined analysis)
_overlap = config.get('GEAxGWAS', {})
# support both new name and legacy 'region_distance'
_overlap_raw = _migrate_assoc_config(_overlap)
_overlap_rdist_raw = _overlap_raw.get('snp_clumping_distance', None)
_gea_is_auto  = CLUMPING_DISTANCE_MODE != 'fixed'
_gwas_is_auto = PHENO_CLUMPING_DISTANCE_MODE != 'fixed'
if _overlap_rdist_raw is None:
    # Default: inherit auto if either source is auto, otherwise use max of both fixed values
    if _gea_is_auto or _gwas_is_auto:
        OVERLAP_CLUMPING_DISTANCE_MODE = CLUMPING_DISTANCE_MODE  # prefer GEA mode
        OVERLAP_CLUMPING_DISTANCE = CLUMPING_DISTANCE
    else:
        OVERLAP_CLUMPING_DISTANCE_MODE = 'fixed'
        OVERLAP_CLUMPING_DISTANCE = max(int(CLUMPING_DISTANCE), int(PHENO_CLUMPING_DISTANCE))
else:
    _raw_str = str(_overlap_rdist_raw).lower().strip()
    if _raw_str in _VALID_REGION_MODES:
        OVERLAP_CLUMPING_DISTANCE_MODE = 'auto_genome_wide' if _raw_str == 'auto' else _raw_str
        OVERLAP_CLUMPING_DISTANCE = _overlap_rdist_raw
    else:
        OVERLAP_CLUMPING_DISTANCE_MODE = 'fixed'
        OVERLAP_CLUMPING_DISTANCE = int(_overlap_rdist_raw)
OVERLAP_CLUMPING_R2_THRESHOLD   = float(_overlap_raw.get('clumping_r2_threshold', CLUMPING_R2_THRESHOLD))
OVERLAP_CLUMPING_LD_DECAY_GROUP = str(_overlap_raw.get('ld_decay_group', CLUMPING_LD_DECAY_GROUP))
# Legacy aliases
OVERLAP_REGION_DISTANCE      = OVERLAP_CLUMPING_DISTANCE
OVERLAP_REGION_DISTANCE_MODE = OVERLAP_CLUMPING_DISTANCE_MODE
OVERLAP_REGION_R2_THRESHOLD  = OVERLAP_CLUMPING_R2_THRESHOLD
OVERLAP_REGION_LD_DECAY_GROUP = OVERLAP_CLUMPING_LD_DECAY_GROUP
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
# PREGEA parameters (hyperparameter exploration mode — LD-pruned only)
#=============================================================================
_pregea = config.get('PreGEA', {})
PREGEA_SEED = _pregea.get('seed', 42)

_pg_lfmm  = _pregea.get('LFMM', {})
PREGEA_LFMM_ENABLED = _pg_lfmm.get('enabled', True)
PREGEA_LFMM_GC      = _pg_lfmm.get('genomic_control', False)   # C.0: FALSE so lambda is informative
PREGEA_LFMM_FDR     = _pg_lfmm.get('fdr', 0.1)
PREGEA_LFMM_BONF    = _pg_lfmm.get('bonf_alpha', 0.05)
PREGEA_LFMM_KS = ([] if K_BEST is None else
    sorted({k for k in range(K_BEST + int(_pg_lfmm.get('k_offset_low', -1)),
                             K_BEST + int(_pg_lfmm.get('k_offset_high', 2)) + 1)
            if k >= int(_pg_lfmm.get('k_min', 1))}))
PREGEA_LFMM_K_STR = ",".join(map(str, PREGEA_LFMM_KS))

_pg_emmax = _pregea.get('EMMAX', {})
PREGEA_EMMAX_ENABLED = _pg_emmax.get('enabled', True)
PREGEA_EMMAX_KINSHIP = _pg_emmax.get('kinship', 'BN')
check_in_list(PREGEA_EMMAX_KINSHIP, ['BN', 'IBS'], 'PreGEA.EMMAX.kinship')
PREGEA_EMMAX_NPC_MIN = int(_pg_emmax.get('n_pcs_min', 0))
PREGEA_EMMAX_NPC_MAX = int(_pg_emmax.get('n_pcs_max', 10))
PREGEA_EMMAX_NPCS    = list(range(PREGEA_EMMAX_NPC_MIN, PREGEA_EMMAX_NPC_MAX + 1))
PREGEA_EMMAX_NPC_STR = ",".join(map(str, PREGEA_EMMAX_NPCS))
PREGEA_EMMAX_FDR  = _pg_emmax.get('fdr', 0.1)
PREGEA_EMMAX_BONF = _pg_emmax.get('bonf_alpha', 0.05)
PREGEA_EMMAX_LAMBDA_TOL      = _pg_emmax.get('lambda_tol', 0.15)
PREGEA_EMMAX_DEFLATION_FLOOR = _pg_emmax.get('deflation_floor', 0.90)

# 'auto' -> Climate.predictors (shared resolver for both ladders)
def _pregea_predictors(block):
    raw = block.get('predictors', 'auto')
    src = PREDICTORS_SELECTED if raw in (None, 'auto') else raw
    return [p.strip() for p in str(src).split(',') if p.strip()]
PREGEA_LFMM_PREDICTORS  = _pregea_predictors(_pg_lfmm)
PREGEA_EMMAX_PREDICTORS = _pregea_predictors(_pg_emmax)

_pg_rda = _pregea.get('RDA', {});      PREGEA_RDA_ENABLED = _pg_rda.get('enabled', True)
PREGEA_RDA_COND_MIN   = int(_pg_rda.get('condition_pcs_min', 0))
PREGEA_RDA_COND_MAX   = int(_pg_rda.get('condition_pcs_max', 6))
PREGEA_RDA_COLLIN_R   = float(_pg_rda.get('collinearity_r', 0.7))
PREGEA_RDA_VIF_MAX    = float(_pg_rda.get('vif_max', 10.0))
PREGEA_RDA_MIN_PREDS  = int(_pg_rda.get('min_predictors', 2))
PREGEA_RDA_AXIS_ALPHA = float(_pg_rda.get('axis_alpha', 0.05))
PREGEA_RDA_PERMS      = int(_pg_rda.get('permutations', 999))

_pg_vp = _pregea.get('Varpart', {});   PREGEA_VP_ENABLED = _pg_vp.get('enabled', True)
PREGEA_VP_RESPONSE = _pg_vp.get('response', 'pcs')
check_in_list(PREGEA_VP_RESPONSE, ['pcs', 'snps'], 'PreGEA.Varpart.response')
PREGEA_VP_STRUCT = _pg_vp.get('structure_table', 'qmatrix')
check_in_list(PREGEA_VP_STRUCT, ['qmatrix', 'none'], 'PreGEA.Varpart.structure_table')
PREGEA_VP_LEVEL = _pg_vp.get('spatial_level', 'auto')
check_in_list(PREGEA_VP_LEVEL, ['auto', 'site', 'sample'], 'PreGEA.Varpart.spatial_level')
PREGEA_VP_MIN_SITES  = int(_pg_vp.get('min_sites', 6))
PREGEA_VP_VAR_CUTOFF = float(_pg_vp.get('response_var_cutoff', 0.8))
PREGEA_VP_MAX_PCS    = int(_pg_vp.get('response_max_pcs', 20))
PREGEA_VP_MIN_PCS    = int(_pg_vp.get('response_min_pcs', 2))
PREGEA_VP_PIN        = float(_pg_vp.get('ordir2step_pin', 0.01))
PREGEA_VP_SELPERMS   = int(_pg_vp.get('ordir2step_permutations', 999))
PREGEA_VP_R2PERMS    = int(_pg_vp.get('r2_permutations', 1000))
PREGEA_VP_PERMS      = int(_pg_vp.get('varpart_permutations', 999))
PREGEA_VP_CONFOUND   = _pg_vp.get('confounding_flag', True)

_pg_tg = _pregea.get('TransferGuard', {})
PREGEA_TRANSFER_GUARD = _pg_tg.get('enabled', False)
PREGEA_TG_LFMM_K  = _pg_tg.get('lfmm_k', 'auto')
PREGEA_TG_EMMAX_NPCS = _pg_tg.get('emmax_n_pcs', 'auto')

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
_rel_tag = f"_rel{PI_HAT}" if RELATEDNESS_REMOVE else ""
FILT_TAG = f"maf{MAF}_miss{MISS}_smiss{SAMPLE_MISS}{_depth_tag}{_het_tag}{_rel_tag}"
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
MOD_PREGEA     = f"{OUTDIR}PreGEA/"
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
    # Coord-valid subset (samples with non-NA latitude/longitude) - feeds climate/GEA/piemap paths
    'coord_valid_samples': f"{INTER}samples/coord_valid_samples.list",
    'metadata_climate': f"{INTER}samples/metadata_climate.tsv",
    # Climate-valid subset: coord-valid samples further narrowed by filter_climate_valid_samples
    # to always exclude samples whose raster extraction returned NA (e.g. ocean/NoData pixel;
    # see download_climate_present.R). Feeds every climate-VALUE-dependent rule (GEA/
    # gradient_forest/geometric_offset/Mantel); coordinate-only rules (IBD, piemaps) stay on the
    # wider coord_valid_samples/metadata_climate. Only applies to CLIMATE_SOURCE 'worldclim'
    # (custom climate can't produce ocean-NA by construction) — for 'custom', these alias
    # directly to the coord-valid paths so consumers don't need source-specific branching.
    'climate_valid_samples': (f"{INTER}samples/climate_valid_samples.list" if CLIMATE_SOURCE == 'worldclim'
                               else f"{INTER}samples/coord_valid_samples.list"),
    'metadata_climate_valid': (f"{INTER}samples/metadata_climate_valid.tsv" if CLIMATE_SOURCE == 'worldclim'
                                else f"{INTER}samples/metadata_climate.tsv"),
    # Row-subset LFMM-format matrices matching the coord-valid sample set (see subset_lfmm_climate
    # in processing.smk). EMMAX/LFMM/gradient_forest/geometric_offset all bind climate values to
    # genotype-matrix rows positionally, so once climate drops coord-missing samples these must match.
    'lfmm_imp_climate': f"{INTER}climate_subset/lfmm_imp_climate.lfmm",
    'lfmm_imp_full_climate': f"{INTER}climate_subset/lfmm_imp_full_climate.lfmm",
    'lfmm_full_climate': f"{INTER}climate_subset/lfmm_full_climate.lfmm",
    # Optional: samples list after het outlier removal (used by filter_vcf when HET_OUTLIER_SD set)
    'samples_het_filtered': f"{INTER}samples/samples_het_filtered.list",
    # Optional: samples list after relatedness (IBD pi-hat) removal (used by filter_vcf when PI_HAT set)
    'samples_rel_filtered': f"{INTER}samples/samples_rel_filtered.list",
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
    # Relatedness (IBD) intermediate — always computed (feeds QC histogram; removal is optional)
    'relatedness_genome': f"{INTER}qc/relatedness.genome",
}

# Output paths (organized by module)
O = {
    'summary': f"{OUTDIR}pipeline_summary.tsv",
    'metadata': f"{MOD_PROCESSING}tables/metadata.tsv",
    'coord_missing_summary': f"{MOD_PROCESSING}tables/coord_missing_summary.tsv",
    'sample_missing_stats': f"{MOD_PROCESSING}tables/sample_missing_stats.tsv",
    'pca': f"{MOD_PRESTRUCT}plots/pca.png",
    'pca_svg': f"{MOD_PRESTRUCT}plots/pca.svg",
    'tracy': f"{MOD_PRESTRUCT}plots/tracy_widom.png",
    'cross_entropy': f"{MOD_PRESTRUCT}plots/cross_entropy_K{K_START}-{K_END}.png",
    # QC intermediate stats
    'qc_raw_summary':      f"{MOD_PROCESSING}tables/vcf_raw_summary.tsv",
    'qc_filtering_summary':f"{MOD_PROCESSING}tables/filtering_summary.tsv",
    'qc_sample_het':       f"{MOD_PROCESSING}tables/sample_heterozygosity.tsv",
    'qc_relatedness_pairs':   f"{MOD_PROCESSING}tables/relatedness_pairs.tsv",
    'qc_relatedness_removed': f"{MOD_PROCESSING}tables/relatedness_removed.tsv",
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
    'qc_plot_relatedness': f"{MOD_PROCESSING}plots/relatedness_distribution.png",
    'qc_plot_relatedness_mds': f"{MOD_PROCESSING}plots/relatedness_mds.png",
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
O['climate_na_excluded']  = _ph('climate_na_excluded')
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
W['vcf_filt_climate']      = _ph('vcf_filt_climate')
W['assoc_tped_climate']    = _ph('assoc_tped_climate')
W['assoc_tfam_climate']    = _ph('assoc_tfam_climate')
W['assoc_kinship_climate'] = _ph('assoc_kinship_climate')
W['assoc_kinship_climate_ibs'] = _ph('assoc_kinship_climate_ibs')  # EMMAX kinship='IBS' flavor
W['gapit_gd']             = _ph('gapit_gd')
W['gapit_gm']             = _ph('gapit_gm')
W['gapit_work']           = _ph('gapit_work')
O['rda_candidates']       = _ph('rda_candidates')
O['rda_diagnostics']      = _ph('rda_diagnostics')
O['rda_anova']            = _ph('rda_anova')
O['rda_screeplot']        = _ph('rda_screeplot')
O['rda_screeplot_svg']    = _ph('rda_screeplot_svg')
O['rda_pval_hist']        = _ph('rda_pval_hist')
O['rda_pval_hist_svg']    = _ph('rda_pval_hist_svg')
O['rda_biplot']           = _ph('rda_biplot')
O['rda_biplot_svg']       = _ph('rda_biplot_svg')
O['selected_snps']        = _ph('selected_snps')
O['regions_per_trait']    = _ph('regions_per_trait')
O['regions_combined']     = _ph('regions_combined')
O['genes_per_region']     = _ph('genes_per_region')
O['genes_per_region_collapsed'] = _ph('genes_per_region_collapsed')
O['genes_combined_regions']     = _ph('genes_combined_regions')
O['manhattan_combined']   = _ph('manhattan_combined')
O['qq_combined']          = _ph('qq_combined')
# --- from add_pheno_association_paths() ---
W['pheno_work']              = _ph('pheno_work')
W['pheno_emmax_work']        = _ph('pheno_emmax_work')
W['pheno_gapit_work']        = _ph('pheno_gapit_work')
W['pheno_tped']              = _ph('pheno_tped')
W['pheno_tfam']              = _ph('pheno_tfam')
W['pheno_kinship']           = _ph('pheno_kinship')
W['pheno_all_phenotypes']    = _ph('pheno_all_phenotypes')
O['pheno_missing_summary']   = _ph('pheno_missing_summary')
O['pheno_selected_snps']     = _ph('pheno_selected_snps')
O['pheno_regions_per_trait'] = _ph('pheno_regions_per_trait')
O['pheno_regions_combined']  = _ph('pheno_regions_combined')
O['pheno_genes_per_region']  = _ph('pheno_genes_per_region')
O['pheno_genes_collapsed']   = _ph('pheno_genes_collapsed')
O['pheno_genes_combined']    = _ph('pheno_genes_combined')
O['pheno_manhattan_combined']= _ph('pheno_manhattan_combined')
O['pheno_qq_combined']       = _ph('pheno_qq_combined')
# --- from add_overlap_paths() ---
O['pairwise_collapsed_snps'] = _ph('pairwise_collapsed_snps')
O['pairwise_overlap_table']  = _ph('pairwise_overlap_table')
O['overlap_miami']           = _ph('overlap_miami')
O['overlap_miami_svg']       = _ph('overlap_miami_svg')
# --- from add_maladaptation_paths() ---
W['climate_future_raster']= _ph('climate_future_raster')
O['climate_future_all']   = _ph('climate_future_all')
O['climate_future_site']  = _ph('climate_future_site')
O['climate_future_na_excluded'] = _ph('climate_future_na_excluded')
O['density_future']       = _ph('density_future')
# --- from add_pregea_paths() ---
for _k in ('pregea_screeplot', 'pregea_screeplot_svg', 'pregea_screeplot_tsv',
           'pregea_lfmm_ladder', 'pregea_lfmm_hist', 'pregea_lfmm_qq',
           'pregea_lfmm_lambda', 'pregea_lfmm_hits',
           'pregea_emmax_ladder', 'pregea_emmax_hist', 'pregea_emmax_qq',
           'pregea_emmax_lambda', 'pregea_emmax_hits',
           'pregea_rda_collin', 'pregea_rda_ladder', 'pregea_rda_axis', 'pregea_rda_fwd',
           'pregea_rda_collin_png', 'pregea_rda_scree_png', 'pregea_rda_ladder_png',
           'pregea_rda_fwd_png', 'pregea_rda_biplot_png',
           'pregea_dbmem_vectors', 'pregea_dbmem_diag', 'pregea_dbmem_png', 'pregea_dbmem_svg',
           'pregea_vp_selection', 'pregea_vp_selected', 'pregea_vp_fractions',
           'pregea_vp_anova', 'pregea_vp_px', 'pregea_vp_path_png', 'pregea_vp_venn_png',
           'pregea_vp_frac_png', 'pregea_vp_px_png',
           'pregea_recommendations', 'pregea_transfer_guard', 'pregea_transfer_png'):
    O[_k] = _ph(_k)

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
    O['climate_na_excluded'] = f"{MOD_CLIMATE}tables/present/climate_na_excluded.tsv"
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

def add_pregea_paths():
    """preGEA-specific paths. K_BEST-dependent because the LFMM K sweep is
    anchored on it and the varpart structure table is the Q-matrix at K_BEST."""
    if K_BEST is None:
        return
    _pl, _tb = f"{MOD_PREGEA}plots/", f"{MOD_PREGEA}tables/"
    # LD-pruned + climate-valid EMMAX chain (does not exist elsewhere — preGEA
    # is the first mode that needs EMMAX-format inputs on the LD-pruned set;
    # WORK_LD prefix avoids collision with gea.smk's WORK_FILT-prefixed rules,
    # so both can coexist in one DAG (needed by the transfer guard).
    W['vcf_ld_climate']     = f"{WORK_LD}climate/{VCF_BASE}.vcf"
    W['pregea_tped']        = f"{WORK_LD}climate/emmax/{VCF_BASE}.tped"
    W['pregea_tfam']        = f"{WORK_LD}climate/emmax/{VCF_BASE}.tfam"
    W['pregea_kinship']     = f"{WORK_LD}climate/emmax/{VCF_BASE}.aBN.kinf"
    W['pregea_kinship_ibs'] = f"{WORK_LD}climate/emmax/{VCF_BASE}.aIBS.kinf"
    O['pregea_screeplot']     = f"{_pl}structure/pruned_pca_screeplot.png"
    O['pregea_screeplot_svg'] = f"{_pl}structure/pruned_pca_screeplot.svg"
    O['pregea_screeplot_tsv'] = f"{_tb}structure/pruned_pca_screeplot.tsv"
    O['pregea_lfmm_ladder']   = f"{_tb}lfmm/lfmm_ladder.tsv"
    O['pregea_lfmm_hist']     = f"{_pl}lfmm/lfmm_pvalue_histogram_grid.png"
    O['pregea_lfmm_qq']       = f"{_pl}lfmm/lfmm_qq_grid.png"
    O['pregea_lfmm_lambda']   = f"{_pl}lfmm/lfmm_lambda_vs_K.png"
    O['pregea_lfmm_hits']     = f"{_pl}lfmm/lfmm_hits_vs_K.png"
    O['pregea_emmax_ladder']  = f"{_tb}emmax/emmax_ladder.tsv"
    O['pregea_emmax_hist']    = f"{_pl}emmax/emmax_pvalue_histogram_grid.png"
    O['pregea_emmax_qq']      = f"{_pl}emmax/emmax_qq_grid.png"
    O['pregea_emmax_lambda']  = f"{_pl}emmax/emmax_lambda_vs_n_pcs.png"
    O['pregea_emmax_hits']    = f"{_pl}emmax/emmax_hits_vs_n_pcs.png"
    O['pregea_rda_collin']    = f"{_tb}rda/rda_predictor_collinearity.tsv"
    O['pregea_rda_ladder']    = f"{_tb}rda/rda_condition_ladder.tsv"
    O['pregea_rda_axis']      = f"{_tb}rda/rda_axis_anova.tsv"
    O['pregea_rda_fwd']       = f"{_tb}rda/rda_ordir2step_path.tsv"
    O['pregea_rda_collin_png']= f"{_pl}rda/rda_predictor_collinearity.png"
    O['pregea_rda_scree_png'] = f"{_pl}rda/rda_axis_screeplot.png"
    O['pregea_rda_ladder_png']= f"{_pl}rda/rda_condition_ladder.png"
    O['pregea_rda_fwd_png']   = f"{_pl}rda/rda_ordir2step_path.png"
    O['pregea_rda_biplot_png']= f"{_pl}rda/rda_biplot_best.png"
    O['pregea_dbmem_vectors'] = f"{_tb}spatial/dbmem_vectors.tsv"
    O['pregea_dbmem_diag']    = f"{_tb}spatial/dbmem_diagnostics.tsv"
    O['pregea_dbmem_png']     = f"{_pl}spatial/dbmem_screeplot.png"
    O['pregea_dbmem_svg']     = f"{_pl}spatial/dbmem_screeplot.svg"
    O['pregea_vp_selection']  = f"{_tb}varpart/dbmem_selection_path.tsv"
    O['pregea_vp_selected']   = f"{_tb}varpart/dbmem_selected.tsv"
    O['pregea_vp_fractions']  = f"{_tb}varpart/varpart_fractions.tsv"
    O['pregea_vp_anova']      = f"{_tb}varpart/varpart_anova.tsv"
    O['pregea_vp_px']         = f"{_tb}varpart/px_per_variable.tsv"
    O['pregea_vp_path_png']   = f"{_pl}varpart/dbmem_selection_path.png"
    O['pregea_vp_venn_png']   = f"{_pl}varpart/varpart_venn.png"
    O['pregea_vp_frac_png']   = f"{_pl}varpart/varpart_fractions_bar.png"
    O['pregea_vp_px_png']     = f"{_pl}varpart/px_barplot.png"
    O['pregea_recommendations']= f"{_tb}pregea_recommendations.tsv"
    O['pregea_transfer_guard'] = f"{_tb}pregea_transfer_guard.tsv"
    O['pregea_transfer_png']   = f"{_pl}transfer/transfer_guard_lambda.png"

add_pregea_paths()

def pregea_kinship_path():
    """BN by default, IBS when configured. Mirrors emmax_kinship_climate_path()."""
    return W['pregea_kinship_ibs'] if PREGEA_EMMAX_KINSHIP == 'IBS' else W['pregea_kinship']

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
        # Coord-valid subset (climate coordinate NA handling, GEA only) — samples with
        # missing lat/lon are dropped for GEA; TPED/kinship rebuilt to match the reduced
        # climate table, since EMMAX binds climate/PCA to genotypes positionally.
        W['vcf_filt_climate']      = f"{WORK_FILT}climate/{VCF_BASE}.vcf"
        W['assoc_tped_climate']    = f"{WORK_FILT}climate/emmax/{VCF_BASE}.tped"
        W['assoc_tfam_climate']    = f"{WORK_FILT}climate/emmax/{VCF_BASE}.tfam"
        W['assoc_kinship_climate'] = f"{WORK_FILT}climate/emmax/{VCF_BASE}.aBN.kinf"
        # GAPIT always consumes the BN kinship above, regardless of EMMAX.kinship —
        # only EMMAX's own kinship becomes flavor-keyed. IBS path is a SEPARATE
        # file (never overwrites the BN one), built by a distinct rule in gea.smk.
        # Default 'BN' -> assoc_kinship_climate is reused verbatim (byte-identical
        # to today; no re-run of existing *_results/).
        if GEA_PARAMS.get("EMMAX", {}).get("kinship") == "IBS":
            W['assoc_kinship_climate_ibs'] = f"{WORK_FILT}climate/emmax/{VCF_BASE}.aIBS.kinf"
    if "gapit" in _assoc_engines:
        W['gapit_gd']   = f"{WORK_FILT}gapit/{VCF_BASE}_GD.tsv"
        W['gapit_gm']   = f"{WORK_FILT}gapit/{VCF_BASE}_GM.tsv"
        W['gapit_work'] = f"{INTER}gapit/gea/"
    # lfmm uses W['lfmm_imp'] / W['lfmm_imp_full'] from add_kbest_paths()
    if "rda" in _assoc_engines:
        _rda_tbl = f"{MOD_GEA}tables/methods/RDA/"
        _rda_plt = f"{MOD_GEA}plots/rda/"
        O['rda_candidates']    = f"{_rda_tbl}RDA_candidates_K{K_BEST}.tsv"
        O['rda_diagnostics']   = f"{_rda_tbl}RDA_diagnostics_K{K_BEST}.tsv"
        O['rda_anova']         = f"{_rda_tbl}RDA_anova_K{K_BEST}.tsv"
        O['rda_screeplot']     = f"{_rda_plt}rda_screeplot_K{K_BEST}.png"
        O['rda_screeplot_svg'] = f"{_rda_plt}rda_screeplot_K{K_BEST}.svg"
        O['rda_pval_hist']     = f"{_rda_plt}rda_pvalue_histogram_K{K_BEST}.png"
        O['rda_pval_hist_svg'] = f"{_rda_plt}rda_pvalue_histogram_K{K_BEST}.svg"
        O['rda_biplot']        = f"{_rda_plt}rda_loadings_biplot_K{K_BEST}.png"
        O['rda_biplot_svg']    = f"{_rda_plt}rda_loadings_biplot_K{K_BEST}.svg"

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

def snp_set_dir(set_name):
    """Directory for a curated SNP set (produced by Shiny, ancestor-less source)."""
    return f"{INTER}snp_sets/{set_name}/"

def snp_set_file(set_name):
    """Path to the curated SNP set TSV (SNPID, chr, pos, min_pvalue)."""
    return f"{INTER}snp_sets/{set_name}/selected_snps.tsv"

def resolve_active_snp_sets():
    """Return the list of active SNP-set names for mode=maladaptation.

    Called ONLY from get_targets('maladaptation') — never at module top level,
    because this function runs for every mode and would crash unrelated modes.
    Raises a friendly ValueError if no sets are found / named sets are missing.
    """
    import glob as _glob
    store = f"{INTER}snp_sets/"
    if SNP_SETS_CFG == 'all' or SNP_SETS_CFG is None:
        found = sorted(
            os.path.basename(os.path.dirname(p))
            for p in _glob.glob(f"{store}*/selected_snps.tsv")
        )
        if not found:
            raise ValueError(
                f"No curated SNP sets found under {store}. "
                "Open the GEA tab in the ADAPTOGENE Shiny app, curate SNPs with your "
                "desired threshold/strategy/regime, and click 'Save SNP set for "
                "maladaptation' before running mode=maladaptation. "
                "(Set Maladaptation.snp_sets to a list of names to select specific sets.)"
            )
        return found
    # explicit list
    names = list(SNP_SETS_CFG)
    if not names:
        raise ValueError(
            "Maladaptation.snp_sets is an empty list. Set it to \"all\" or provide "
            ">=1 saved set name from the GEA tab."
        )
    missing = [n for n in names if not os.path.exists(f"{store}{n}/selected_snps.tsv")]
    if missing:
        raise ValueError(
            f"Maladaptation.snp_sets names not found under {store}: {missing}. "
            "Save them in the GEA tab or fix the names in config."
        )
    return names

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

def mala_offset_piemap_points(method, run_label, spatial_tag):
    """Points ('clear map') companion — trait-independent, one per method/run_label/spatial_tag."""
    return f"{mala_plot_dir(method, run_label, spatial_tag)}genetic_offset_piemap_points.png"

# Maladaptation paths
def add_maladaptation_paths():
    """Add maladaptation-specific paths to W and O dictionaries."""
    # Custom-climate path needs no MODELS_LIST — skip only when worldclim and models absent
    if K_BEST is None or (CLIMATE_SOURCE == 'worldclim' and not MODELS_LIST):
        return

    # Per-model future climate rasters (worldclim only — custom stages a single future table)
    if CLIMATE_SOURCE == 'worldclim':
        for model in MODELS_LIST:
            W[f'climate_future_{model}'] = f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_{model}.tif"

    # Merged future climate (both sources — identical output paths)
    W['climate_future_raster'] = f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_rasterstack.tif"
    O['climate_future_all'] = f"{MOD_CLIMATE}tables/future/climate_future_year{YEAR}_ssp{SSP}_all.tsv"
    O['climate_future_site'] = f"{MOD_CLIMATE}tables/future/climate_future_year{YEAR}_ssp{SSP}_site.tsv"
    O['climate_future_na_excluded'] = f"{MOD_CLIMATE}tables/future/climate_future_year{YEAR}_ssp{SSP}_na_excluded.tsv"

    # Future climate density plot (method-agnostic)
    O['density_future'] = f"{MOD_CLIMATE}plots/density_plot_future_ssp{SSP}_{YEAR}.png"

    # Snakemake auto-creates declared-output parent dirs; no makedirs needed here.

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
def piemap_notrait_points(bio): return f"{MOD_STRUCT}plots/piemap/piemap_{bio}_points.png"

# Templates for association outputs
def assoc_pvalues(method): return f"{MOD_GEA}tables/methods/{method}/{method}_pvalues_K{K_BEST}.tsv"
def assoc_sigsnps(method, adjust): return f"{MOD_GEA}tables/methods/{method}/{method}_pvalues_K{K_BEST}_sig_snps_{adjust}.tsv"

def emmax_kinship_climate_path():
    """The kinship file the GEA EMMAX rule actually consumes — BN (default,
    W['assoc_kinship_climate'], byte-identical to pre-hyperparameter-system
    behavior) or IBS (W['assoc_kinship_climate_ibs'], a separate file, never
    overwriting the BN one) per EMMAX.kinship."""
    if GEA_PARAMS.get("EMMAX", {}).get("kinship") == "IBS":
        return W['assoc_kinship_climate_ibs']
    return W['assoc_kinship_climate']
def manhattan_plot(method, trait, adjust): return f"{MOD_GEA}plots/manhattan/{method}/manhattan_{trait}_K{K_BEST}_{adjust}.png"
def qq_plot(method, trait, adjust): return f"{MOD_GEA}plots/manhattan/{method}/qq_{trait}_K{K_BEST}_{adjust}.png"

# Per-factor intermediate cache (input-driven: K in name ensures cascade on K change)
def assoc_pvalues_trait(method, trait):
    """Per-trait intermediate p-value file under _intermediate/gea_per_trait/.
    Named with K so that changing K changes this path, invalidating the per-trait cache."""
    return f"{INTER}gea_per_trait/{method}/{trait}_pvalues_K{K_BEST}.tsv"

def pheno_pvalues_trait(method, trait):
    """Per-trait intermediate p-value file for GWAS Path A (MEAN/MEDIAN mode).
    Named with K so that changing K changes this path, invalidating the per-trait cache."""
    return f"{INTER}gwas_per_trait/{method}/{trait}_pvalues_K{K_BEST}.tsv"

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

# --- Registry-derived helpers for multivariate ("pseudo-trait") methods like RDA ---
RDA_MODELS = {n for n, c in GEA_METHODS.items() if c["engine"] == "rda"}

def method_pseudo_trait(method):
    """The single p-value column name a multivariate method emits, or None."""
    return GEA_METHODS.get(method, {}).get("pseudo_trait")

def method_traits(source, method):
    """Trait/column list a given (source, method) actually emits in its wide
    p-value table. Multivariate methods (RDA) emit exactly one pseudo-trait
    column instead of one column per configured predictor — using the plain
    predictor list for them would request e.g. manhattan_bio_1_* targets that
    plot_manhattan.R hard-stops on ('Trait not found')."""
    pt = method_pseudo_trait(method)
    if pt:
        return [pt]
    return get_predictors_list() if source == "GEA" else PHENO_TRAITS

def wza_methods(source):
    """Methods (of a source) whose p-values should be fed through compute_wza.R."""
    configs = GEA_CONFIGS if source == "GEA" else GWAS_CONFIGS
    return [m for m in configs if GEA_METHODS.get(m, {}).get("supports_wza", True)]

def combine_predictors(source):
    """PREDICTORS_SELECTED (or PHENO_PREDICTORS) + the pseudo-traits of every
    configured multivariate method for this source, comma-joined. Used ONLY as
    the predictor-name argument threaded into combine_selected_snps.R and the
    Manhattan plot scripts — without this, a multivariate method's pseudo-trait
    (e.g. RDA's climate_multivariate) is invisible to trait-name filters that
    otherwise assume every column name is a literal configured predictor."""
    base_str = PREDICTORS_SELECTED if source == "GEA" else PHENO_PREDICTORS
    base = [p.strip() for p in (base_str or '').split(',') if p.strip()]
    configs = GEA_CONFIGS if source == "GEA" else GWAS_CONFIGS
    extra = [pt for m in configs if (pt := method_pseudo_trait(m))]
    return ",".join(base + extra)

# WZA config parsing
def _parse_wza_config(group_dict, defaults=None):
    """Parse WZA sub-config (window_size, fallback_window_bp)."""
    defaults = defaults or {}
    wza = group_dict.get('wza', {})
    wza_def = defaults.get('wza', {})
    raw_ws = wza.get('window_size', wza_def.get('window_size', 'auto_genome_wide'))
    raw_str = str(raw_ws).lower().strip()
    if raw_str in _VALID_REGION_MODES or raw_str == 'auto_genome_wide':
        ws_mode = 'auto_genome_wide' if raw_str in ('auto', 'auto_genome_wide') else raw_str
        ws_val  = raw_str
    else:
        ws_mode = 'fixed'
        ws_val  = int(raw_ws)
    return {
        'window_size':       ws_val,
        'window_size_mode':  ws_mode,
        'fallback_window_bp': int(wza.get('fallback_window_bp', wza_def.get('fallback_window_bp', 10000))),
    }

_gea_wza = _parse_wza_config(_assoc)

if K_BEST is not None and GEA_CONFIGS:
    ASSOC_SOURCES["GEA"] = {
        "mod":             MOD_GEA,
        "logdir":          f"{LOGDIR}GEA/",
        "configs":         GEA_CONFIGS,
        "method_regex":    GEA_METHOD_REGEX,
        "trait_regex":     r"bio_\d+",
        "predictors":      PREDICTORS_SELECTED,
        "combine_predictors": combine_predictors("GEA"),
        "params":          GEA_PARAMS,
        "clumping_distance":      CLUMPING_DISTANCE,
        "clumping_distance_mode": CLUMPING_DISTANCE_MODE,
        "clumping_r2_threshold":  CLUMPING_R2_THRESHOLD,
        "ld_decay_group":         CLUMPING_LD_DECAY_GROUP,
        "promoter_length":        PROMOTER_LENGTH,
        "combine_method":  SIGSNPS_METHOD,
        "wza":             _gea_wza,
        "pvalues_fn":      assoc_pvalues,
        "sigsnps_fn":      assoc_sigsnps,
    }

_gwas_wza = _parse_wza_config(_pheno, defaults=_assoc)

if K_BEST is not None and GWAS_CONFIGS:
    ASSOC_SOURCES["GWAS"] = {
        "mod":             MOD_GWAS,
        "logdir":          f"{LOGDIR}GWAS/",
        "configs":         GWAS_CONFIGS,
        "method_regex":    GWAS_METHOD_REGEX,
        "trait_regex":     r"[a-zA-Z]\w*",
        "predictors":      PHENO_PREDICTORS,
        "combine_predictors": combine_predictors("GWAS"),
        "params":          GWAS_PARAMS,
        "clumping_distance":      PHENO_CLUMPING_DISTANCE,
        "clumping_distance_mode": PHENO_CLUMPING_DISTANCE_MODE,
        "clumping_r2_threshold":  PHENO_CLUMPING_R2_THRESHOLD,
        "ld_decay_group":         PHENO_CLUMPING_LD_DECAY_GROUP,
        "promoter_length":        PHENO_PROMOTER_LENGTH,
        "combine_method":  PHENO_COMBINE_METHOD,
        "wza":             _gwas_wza,
        "pvalues_fn":      pheno_pvalues,
        "sigsnps_fn":      pheno_sigsnps,
    }

def _src(source, key):
    """Look up a source-indexed config value in unified downstream rule bodies."""
    return ASSOC_SOURCES[source][key]

# Wildcard regex matching all configured sources (used in _assoc_downstream.smk)
SOURCE_REGEX = "|".join(ASSOC_SOURCES.keys()) if ASSOC_SOURCES else "gea"
# Trait wildcard: union of GEA (bio_\d+) and GWAS (identifier) patterns.
# {source} in the output path disambiguates which branch a match belongs to.
TRAIT_REGEX_ANY = r"(?!wza_)(bio_\d+|[a-zA-Z]\w*)"


def assoc_out(source, key):
    """Return the per-source output path for a logical downstream key."""
    _templates = {
        "selected_snps":              "tables/selected_snps.tsv",
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
        # WZA combined outputs
        "wza_combined_png":           f"plots/manhattan/combined/manhattan_wza_combined_K{K_BEST}.png",
        "wza_combined_svg":           f"plots/manhattan/combined/manhattan_wza_combined_K{K_BEST}.svg",
        "wza_qq_combined_png":        f"plots/manhattan/combined/qq_wza_combined_K{K_BEST}.png",
        "wza_qq_combined_svg":        f"plots/manhattan/combined/qq_wza_combined_K{K_BEST}.svg",
        "wza_combined_bg":            f"plots/manhattan/combined/manhattan_wza_combined_K{K_BEST}_background.png",
        "wza_combined_coords":        f"plots/manhattan/combined/manhattan_wza_combined_K{K_BEST}_coords.json",
    }
    return f"{ASSOC_SOURCES[source]['mod']}{_templates[key]}"


def wza_pvalues(source, method):
    """WZA p-values TSV path for a (source, method)."""
    return f"{ASSOC_SOURCES[source]['mod']}tables/methods/{method}/{method}_wza_K{K_BEST}.tsv"

def wza_sig_windows(source, method, adjust):
    """WZA significant windows TSV path."""
    return f"{ASSOC_SOURCES[source]['mod']}tables/methods/{method}/{method}_wza_K{K_BEST}_sig_windows_{adjust}.tsv"

def wza_manhattan(source, method, trait, adjust):
    """WZA per-method Manhattan background PNG path."""
    return f"{ASSOC_SOURCES[source]['mod']}plots/manhattan/{method}/manhattan_wza_{trait}_K{K_BEST}_{adjust}_background.png"


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

# preGEA directories (LD-pruned hyperparameter-exploration mode)
dirs_to_create += [
    f"{MOD_PREGEA}plots/structure/", f"{MOD_PREGEA}plots/lfmm/", f"{MOD_PREGEA}plots/emmax/",
    f"{MOD_PREGEA}plots/rda/", f"{MOD_PREGEA}plots/spatial/", f"{MOD_PREGEA}plots/varpart/",
    f"{MOD_PREGEA}plots/transfer/",
    f"{MOD_PREGEA}tables/structure/", f"{MOD_PREGEA}tables/lfmm/", f"{MOD_PREGEA}tables/emmax/",
    f"{MOD_PREGEA}tables/rda/", f"{MOD_PREGEA}tables/spatial/", f"{MOD_PREGEA}tables/varpart/",
    f"{INTER}pregea/lfmm/", f"{INTER}pregea/emmax/",
    f"{LOGDIR}pregea/",
]
if K_BEST is not None:
    dirs_to_create += [f"{INTER}pregea/lfmm/K{k}/" for k in PREGEA_LFMM_KS]
    dirs_to_create += [f"{INTER}pregea/emmax/npc{n}/" for n in PREGEA_EMMAX_NPCS]

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

# Per-factor intermediate cache dirs (created for each configured GEA method)
if GEA_CONFIGS and K_BEST is not None:
    for _m in GEA_CONFIGS:
        dirs_to_create.append(f"{INTER}gea_per_trait/{_m}/")

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
    # Per-factor intermediate cache dirs for Path A GWAS (MEAN/MEDIAN)
    if K_BEST is not None and PHENO_MISSING != 'DROP':
        for _m in GWAS_CONFIGS:
            dirs_to_create.append(f"{INTER}gwas_per_trait/{_m}/")
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
    targets = []
    for method, adjust in src["configs"].items():
        # Per-method trait/column list — multivariate methods (RDA) emit a single
        # pseudo-trait column, not one column per configured predictor.
        m_traits = method_traits(source, method)
        # Per-SNP outputs
        targets.append(src["pvalues_fn"](method))
        targets.append(src["sigsnps_fn"](method, adjust))
        for trait in m_traits:
            targets.append(f"{src['mod']}plots/manhattan/{method}/manhattan_{trait}_K{K_BEST}_{adjust}.png")
            targets.append(f"{src['mod']}plots/manhattan/{method}/qq_{trait}_K{K_BEST}_{adjust}.png")
        # WZA outputs — opt-out via the method's supports_wza registry flag
        if GEA_METHODS.get(method, {}).get("supports_wza", True):
            targets.append(f"{OUTDIR}{source}/tables/methods/{method}/{method}_wza_K{K_BEST}.tsv")
            targets.append(f"{OUTDIR}{source}/tables/methods/{method}/{method}_wza_K{K_BEST}_sig_windows_{adjust}.tsv")
            for trait in m_traits:
                targets.append(f"{src['mod']}plots/manhattan/{method}/manhattan_wza_{trait}_K{K_BEST}_{adjust}.png")
                targets.append(f"{src['mod']}plots/manhattan/{method}/qq_wza_{trait}_K{K_BEST}_{adjust}.png")
        # RDA-specific side outputs (diagnostics/candidates/anova/plots) — see
        # add_association_paths()'s "rda" engine block for these O[] keys.
        if method in RDA_MODELS:
            targets += [O['rda_candidates'], O['rda_diagnostics'], O['rda_anova'],
                        O['rda_screeplot'], O['rda_pval_hist'], O['rda_biplot']]
    for key in (
        "selected_snps", "regions_per_trait", "regions_combined",
        "genes_per_region", "genes_per_region_collapsed", "genes_combined",
        "manhattan_combined_png", "qq_combined_png",
    ):
        targets.append(assoc_out(source, key))
    # Combined WZA targets only if at least one configured method supports WZA
    if wza_methods(source):
        targets.append(assoc_out(source, "wza_combined_png"))
        targets.append(assoc_out(source, "wza_qq_combined_png"))
    return targets


def get_targets(mode):
    if mode == 'processing':
        targets = [
            W['samples_missing_stats'], W['samples_removed'],  # Sample missingness outputs
            W['vcf_filt'], W['vcf_ld'], W['geno'], W['lfmm'], O['metadata'],
            O['coord_missing_summary'],
            # QC outputs
            O['qc_raw_summary'], O['qc_filtering_summary'], O['qc_sample_het'],
            O['qc_plot_sample_miss'], O['qc_plot_het_miss'],
            O['qc_plot_maf'], O['qc_plot_snp_miss'],
            O['qc_plot_attrition'], O['qc_plot_snp_density'],
            O['qc_plot_relatedness'], O['qc_plot_relatedness_mds'], O['qc_relatedness_pairs'],
            W['summary_done']
        ]
        if GFF:
            targets.append(W['gff_normalized'])
        if HAS_FORMAT_DP:
            targets.append(O['qc_plot_depth'])
        if PI_HAT is not None:
            # Always compute the would-remove list (even in 'keep' mode) — the Shiny
            # hover note needs it to preview the removal count before the user opts in.
            targets.append(O['qc_relatedness_removed'])
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

    elif mode == 'pregea':
        check_numeric(K_BEST, 'K_BEST')
        if not CLIMATE_ENABLED:
            raise ValueError("pregea mode requires Climate.enabled: true "
                             "(every ladder uses climate predictors)")
        if not get_predictors_list():
            raise ValueError("Climate.predictors must be set for pregea mode")
        if not (PREGEA_LFMM_ENABLED or PREGEA_EMMAX_ENABLED or
                PREGEA_RDA_ENABLED or PREGEA_VP_ENABLED):
            raise ValueError("pregea mode requires at least one PreGEA block enabled")
        targets = [O['pregea_screeplot'], O['pregea_screeplot_tsv']]
        if PREGEA_LFMM_ENABLED:
            if not PREGEA_LFMM_KS:
                raise ValueError("PreGEA.LFMM: empty K range (check k_offset_low/high/k_min)")
            targets += [O['pregea_lfmm_ladder'], O['pregea_lfmm_hist'],
                        O['pregea_lfmm_qq'], O['pregea_lfmm_lambda'], O['pregea_lfmm_hits']]
        if PREGEA_EMMAX_ENABLED:
            targets += [O['pregea_emmax_ladder'], O['pregea_emmax_hist'],
                        O['pregea_emmax_qq'], O['pregea_emmax_lambda'], O['pregea_emmax_hits']]
        if PREGEA_RDA_ENABLED:
            targets += [O['pregea_rda_collin'], O['pregea_rda_ladder'],
                        O['pregea_rda_axis'], O['pregea_rda_fwd'],
                        O['pregea_rda_collin_png'], O['pregea_rda_scree_png'],
                        O['pregea_rda_ladder_png'], O['pregea_rda_fwd_png'],
                        O['pregea_rda_biplot_png']]
        if PREGEA_VP_ENABLED:
            targets += [O['pregea_dbmem_vectors'], O['pregea_dbmem_diag'],
                        O['pregea_dbmem_png'],
                        O['pregea_vp_selection'], O['pregea_vp_selected'],
                        O['pregea_vp_fractions'], O['pregea_vp_anova'], O['pregea_vp_px'],
                        O['pregea_vp_path_png'], O['pregea_vp_venn_png'],
                        O['pregea_vp_frac_png'], O['pregea_vp_px_png']]
        targets.append(O['pregea_recommendations'])
        if PREGEA_TRANSFER_GUARD:
            targets += [O['pregea_transfer_guard'], O['pregea_transfer_png']]
        targets.append(W['summary_done'])
        return targets

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
            targets += [piemap_notrait_points(bio) for bio in ALL_BIO]
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
        if CLIMATE_SOURCE == 'worldclim' and not MODELS_LIST:
            raise ValueError("MODELS must be set for maladaptation mode (worldclim source requires GCM models)")
        # GF-specific checks: only validate when GF is active
        if 'gradient_forest' in ACTIVE_MALA_METHODS:
            check_numeric(NTREE, 'NTREE')
            check_float(COR_THRESHOLD, 'COR_THRESHOLD')

        # Resolve saved SNP sets (glob/validate happens here, not at module top level)
        ACTIVE_SNP_SETS = resolve_active_snp_sets()

        targets = [
            # Future climate (method-agnostic)
            O['climate_future_site'],
            O['climate_future_all'],
            W['climate_future_raster'],
            O['density_future'],
        ]

        for set_name in ACTIVE_SNP_SETS:
            for method in ACTIVE_MALA_METHODS:
                _mflags = MALADAPTATION_METHODS[method]
                for spatial_tag in mala_spatial_tags(method):
                    # Separate model artifact (GF only — geometric_offset is single-call)
                    if _mflags['builds_model']:
                        targets.append(mala_model(method, set_name, spatial_tag, 'adaptive'))
                    # Core offset outputs (all methods)
                    targets += [
                        mala_offset_map_values(method, set_name, spatial_tag),
                        mala_offset_site_values(method, set_name, spatial_tag),
                        mala_importance(method, set_name, spatial_tag),
                        mala_offset_piemap(method, set_name, spatial_tag, 'notrait'),
                        mala_offset_piemap_points(method, set_name, spatial_tag),
                    ]
                    # Cumulative importance (GF only)
                    if _mflags['supports_cumulative_importance']:
                        targets.append(mala_cumimp(method, set_name, spatial_tag))
                    # Random/neutral model (GF + config flag)
                    if _mflags['supports_random_model'] and GF_RANDOM_MODEL:
                        targets.append(mala_model(method, set_name, spatial_tag, 'random'))
                    # Pop-stats scaled piemaps (all methods, optional)
                    if CALC_POP_STATS:
                        targets += [
                            mala_offset_piemap(method, set_name, spatial_tag, 'tajima_d'),
                            mala_offset_piemap(method, set_name, spatial_tag, 'pi_diversity'),
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
            # Points ("clear map") companion — trait-independent, one per project
            targets.append(f"{MOD_GWAS}plots/piemap/phenomap_points.png")

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
        raise ValueError("Specify mode: --config mode=processing or mode=prestructure or "
                         "mode=pregea or mode=structure or mode=gea or mode=gwas or "
                         "mode=gea_x_gwas or mode=maladaptation")
    else:
        raise ValueError(f"Unknown mode: {mode}")

# Helper: provide LD decay table as input when distance mode is auto_*
def ld_decay_input(mode):
    """Return LD decay table path as input dependency when mode is auto_* or auto."""
    needs_table = isinstance(mode, str) and mode.lower() in _VALID_REGION_MODES
    return O.get('ld_decay_table', []) if needs_table else []

def wza_ld_decay_input(wza_cfg):
    """Return LD decay table as input for WZA when window_size is auto_*."""
    needs = isinstance(wza_cfg.get('window_size', ''), str) and \
            wza_cfg.get('window_size_mode', 'fixed') != 'fixed'
    return O.get('ld_decay_table', []) if needs else []

#=============================================================================
# MAIN RULE
