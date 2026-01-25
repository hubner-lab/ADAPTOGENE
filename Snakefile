# ADAPTOGENE Pipeline - Refactored
# vim: filetype=python
import os
from pathlib import Path

#=============================================================================
# CONFIGURATION
#=============================================================================
configfile: "/pipeline/config.yaml"

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
CUSTOM_TRAIT = config.get('CUSTOM_TRAIT', 'NULL')

# PieMap plot parameters
PIEMAP_PALETTE = config.get('PieMap_PALLETE_MAP', 1022614)
PIEMAP_PALETTE_REV = config.get('PieMap_PALLETE_MAP_reverse', False)
PIEMAP_PIE_SIZE = config.get('PieMap_PIE_SIZE', 0.1)
PIEMAP_PIE_ALPHA = config.get('PieMap_PIE_ALPHA', 0.6)
PIEMAP_IBD_ALPHA = config.get('PieMap_IBD_ALPHA', 0.6)
PIEMAP_IBD_COLOR = config.get('PieMap_IBD_COLOR', 'black')
PIEMAP_POP_LABEL = config.get('PieMap_POP_LABEL', 'F')
PIEMAP_PIE_RESCALE = config.get('PieMap_PIE_RESCALE', 1)
PIEMAP_POP_LABEL_SIZE = config.get('PieMap_POP_LABEL_SIZE', 10)

#=============================================================================
# PATH DEFINITIONS
#=============================================================================
# Directory tags (easy to modify if adding new parameters)
FILT_TAG = f"maf{MAF}_miss{MISS}"
LD_TAG = f"ld{LD_R2}"

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
    'samples_order': f"{INTER}samples_order.list",
    # Filtered VCF (in FILT_TAG directory)
    'vcf_filt': f"{WORK_FILT}{VCF_BASE}.vcf",
    # LD-pruned files (in LD_TAG subdirectory)
    'vcf_ld': f"{WORK_LD}{VCF_BASE}.vcf",
    'eigenvec': f"{WORK_LD}{VCF_BASE}.eigenvec",
    'eigenval': f"{WORK_LD}{VCF_BASE}.eigenval",
    'prune_in': f"{WORK_LD}{VCF_BASE}.prune.in",
    'geno': f"{WORK_LD}{VCF_BASE}.geno",
    'lfmm': f"{WORK_LD}{VCF_BASE}.lfmm",
    'vcfsnp': f"{WORK_LD}{VCF_BASE}.vcfsnp",
    'removed': f"{WORK_LD}{VCF_BASE}.removed",
    'snmf': f"{WORK_LD}{VCF_BASE}.snmfProject",
}

# Output paths (organized results)
O = {
    'metadata': f"{TABLES}metadata.tsv",
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
    W['climate_raster'] = f"{INTER}climate_present.grd"
    
    # Tables
    O['climate_site'] = f"{TABLES}climate_present_site.tsv"
    O['climate_site_scaled'] = f"{TABLES}climate_present_site_scaled.tsv"
    O['climate_all'] = f"{TABLES}climate_present_all.tsv"
    O['tajima'] = f"{TABLES}TajimaD_byPop.tsv"
    O['pi_div'] = f"{TABLES}Pi_diversity_byPop.tsv"
    O['ibd_raw'] = f"{TABLES}IBD_raw.tsv"
    O['ibd_pairs'] = f"{TABLES}IBD_notIsolated.tsv"
    O['amova'] = f"{TABLES}AMOVA.tsv"
    
    # Plots
    O['corr_heatmap'] = f"{PLOTS}climate/CorrelationHeatmap.png"
    O['mantel'] = f"{PLOTS}structure/MantelTest.png"
    O['amova_plot'] = f"{PLOTS}structure/AMOVA.png"

add_kbest_paths()

# Templates for K-dependent outputs
def clusters_table(k): return f"{TABLES}clusters_K{k}.tsv"
def structure_plot(k): return f"{PLOTS}structure/structure_K{k}.png"
def pca_struct_plot(k): return f"{PLOTS}pca/pca_structure_K{k}.png"
def pop_diff_plot(k): return f"{PLOTS}structure/pop_diff_K{k}.png"

# Templates for climate/trait-dependent outputs
def density_plot(bio): return f"{PLOTS}climate/DensityPlot_{bio}.png"
def piemap_tajima(bio): return f"{PLOTS}piemap/PieMap_{bio}_TajimaD.png"
def piemap_diversity(bio): return f"{PLOTS}piemap/PieMap_{bio}_PiDiversity.png"
def piemap_notrait(bio): return f"{PLOTS}piemap/PieMap_{bio}.png"

#=============================================================================
# CREATE DIRECTORIES
#=============================================================================
for d in [WORK, WORK_FILT, WORK_LD, PLOTS, 
          f"{PLOTS}pca/", f"{PLOTS}structure/", f"{PLOTS}climate/", f"{PLOTS}piemap/",
          TABLES, INTER, LOGDIR]:
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
        return [W['vcf_filt'], W['vcf_ld'], W['geno'], W['lfmm'], W['eigenvec'], O['metadata']]
    
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
        
        return (
            # Imputed data
            [W['lfmm_imp'], W['vcf_imp']] +
            # Climate data
            [O['climate_site'], O['climate_site_scaled'], O['climate_all'], W['climate_raster']] +
            # Density plots for each predictor
            [density_plot(bio) for bio in predictors] +
            # Population metrics
            [O['tajima'], O['pi_div'], O['ibd_raw'], O['ibd_pairs']] +
            # Summary plots
            [O['corr_heatmap'], O['mantel'], O['amova'], O['amova_plot']] +
            # PieMaps for each predictor
            [piemap_tajima(bio) for bio in predictors] +
            [piemap_diversity(bio) for bio in predictors]
        )
    
    elif mode is None:
        raise ValueError("Specify mode: --config mode=processing or mode=structure or mode=structure_K")
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

rule filter_vcf:
    """Filter VCF by MAF, missingness, and sample list."""
    input:  vcf = f"{INDIR}{VCF_RAW}", samples = W['samples_list']
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
    """LD prune and compute PCA."""
    input:  vcf = W['vcf_filt']
    output: vcf = W['vcf_ld'], eigenvec = W['eigenvec'], eigenval = W['eigenval'], prune = W['prune_in']
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
            --make-bed --pca --recode vcf \
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
    """Generate PCA and Tracy-Widom plots."""
    input:  lfmm = W['lfmm'], meta = O['metadata']
    output: pca = O['pca'], pca_svg = O['pca_svg'], tracy = O['tracy']
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
    """Generate PCA with structure pie charts for specific K."""
    input:
        lfmm = W['lfmm'],
        clusters = clusters_table("{k}"),
        eigenvec = W['eigenvec'],
        eigenval = W['eigenval']
    output: pca_struct_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k, plot_dir = f"{PLOTS}pca/", inter_dir = INTER
    log:    f"{LOGDIR}pca_structure_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pca_structure.R \
            {input.lfmm} {input.clusters} {input.eigenvec} {input.eigenval} \
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
        tables_dir = TABLES
    log: f"{LOGDIR}download_climate_present.log"
    shell:
        """
        Rscript /pipeline/scripts/download_climate_present.R \
            {input.meta} {params.crop} {params.gap} {params.data_dir} {params.resolution} \
            {params.inter_dir} {params.tables_dir} > {log} 2>&1
        """

rule density_plot:
    """Generate density plot for a climate predictor."""
    input:  climate = O['climate_site']
    output: density_plot("{bio}")
    wildcard_constraints: bio = r"bio_\d+"
    params: bio = lambda wc: wc.bio, plot_dir = f"{PLOTS}climate/", inter_dir = INTER
    log:    f"{LOGDIR}density_plot_{{bio}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_density.R \
            {input.climate} {params.bio} {params.plot_dir} {params.inter_dir} > {log} 2>&1
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
    params: tables_dir = TABLES
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
    params: plot_dir = f"{PLOTS}climate/", inter_dir = INTER
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
    params: plot_dir = f"{PLOTS}structure/", tables_dir = TABLES, inter_dir = INTER
    log:    f"{LOGDIR}amova.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/amova.R \
            {input.vcf} {input.meta} {threads} \
            {params.plot_dir} {params.tables_dir} {params.inter_dir} > {log} 2>&1
        """

rule piemap_plot:
    """Generate PieMap plots for a climate predictor."""
    input:
        meta = O['metadata'],
        clusters = clusters_table(K_BEST),
        raster = W['climate_raster'],
        ibd = O['ibd_pairs'],
        tajima = O['tajima'],
        diversity = O['pi_div']
    output:
        tajima_plot = piemap_tajima("{bio}"),
        diversity_plot = piemap_diversity("{bio}")
    wildcard_constraints: bio = r"bio_\d+"
    params:
        bio = lambda wc: wc.bio,
        custom_trait = f"{INDIR}{CUSTOM_TRAIT}",
        palette = PIEMAP_PALETTE,
        palette_rev = PIEMAP_PALETTE_REV,
        pie_size = PIEMAP_PIE_SIZE,
        pie_alpha = PIEMAP_PIE_ALPHA,
        ibd_alpha = PIEMAP_IBD_ALPHA,
        ibd_color = PIEMAP_IBD_COLOR,
        pop_label = PIEMAP_POP_LABEL,
        pie_rescale = PIEMAP_PIE_RESCALE,
        pop_label_size = PIEMAP_POP_LABEL_SIZE,
        plot_dir = f"{PLOTS}piemap/",
        inter_dir = INTER
    log: f"{LOGDIR}piemap_{{bio}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.meta} {input.clusters} {input.raster} {input.ibd} \
            {input.tajima} {input.diversity} {params.custom_trait} \
            {params.bio} {params.palette} {params.palette_rev} \
            {params.pie_size} {params.pie_alpha} {params.ibd_alpha} {params.ibd_color} \
            {params.pop_label} {params.pie_rescale} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """
