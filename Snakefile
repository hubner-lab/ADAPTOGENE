#TODO train_N and SNPID_N count the same - SNPs number
# vim: filetype=python
import os
from collections import defaultdict
from itertools import product
import re
from pathlib import Path


# Create dirs
########################### Functions
def check_numeric(value, arg_name, allowed_null = False):
	if allowed_null:
		return True	
	try:
		int(value)
	except:
		raise ValueError(f'The {arg_name} is not numeric. Please provide correct number')

def check_in_list(value, allowed_values_lst, arg_name):
	if value in allowed_values_lst:
		return True
	else:
		allowed_values_str = '/'.join(allowed_values_lst)
		raise ValueError(f"The {arg_name} has not proper value. Should be on of the {allowed_values_str}.")

def check_float(value, arg_name, allowed_null = False):
	if allowed_null:
		return True
	try:
		float(value)
	except:
		raise ValueError(f'The {arg_name} is not float. Please provide correct number')

def check_comma_separated(value, arg_name, numeric = False, allowed_single_values = None):
	if allowed_single_values is not None and value in allowed_single_values:
		return True
	if " " in value:
		raise ValueError(f"Invalid input: {value} for argument {arg_name}. Spaces are not allowed.")
	if ',' not in value:
		raise ValueError(f"Invalid input: {value} for argument {arg_name}. Input must be a comma-separated string with at least 2 arguments.")
	if numeric:
		try:
			[check_float(x, arg_name) for x in value.split(',')]
		except:
			raise ValueError(f"Invalid input: {value} for argument {arg_name}. Input must contain only numbers separated by commas")

def check_vcf_name(dir, value, arg_name):
	if not value.endswith('.vcf') and not value.endswith('vcf.gz'):
		raise ValueError(f"Please provide file for {arg_name} in .vcf or .vcf.gz format")

def check_file_exist(dir_path, value, arg_name):
	
	file_path = os.path.join(dir_path, value)
	
	if not os.path.exists(file_path):
		raise ValueError(f"File {value} doesn't exist in {dir_path}, please provide correct file name for argument {arg_name} with correct INDIR argument")
 
# Parse association config list
# Parse association config list
def parse_association_configs(config, varname = 'ASSOCIATION_CONFIGS'):
    """Parse association configurations from config file"""
    assoc_configs = config.get(varname, [])
    
    # If no configs specified, use legacy parameter
    configs = defaultdict(str)
    method_seen = set()
    
    for cfg in assoc_configs:
        method = cfg["METHOD"]
        adjust = cfg["ADJUST"] + '_' + cfg["THRESHOLD"]
        
        # Check if method already exists
        if method in method_seen:
            raise ValueError(f"Method '{method}' appears multiple times in ASSOCIATION_CONFIGS. Only one ADJUST+THRESHOLD combination per method is allowed.")
        
        method_seen.add(method)
        configs[method] = adjust
    
    return configs



# Generate targets for all association configurations
def get_association_targets(assoc_configs, config):  # Accept pre-parsed configs
    targets = []
    # Remove the internal parsing line
    # assoc_configs = parse_association_configs(config)  # DELETE THIS LINE
    
    for method in assoc_configs:
        if method == "EMMAX":
            base = f"tables/EMMAX/EMMAX_pvalues_K{K_BEST}"
        elif method == "LFMM":
            base = f"tables/LFMM/LFMM_pvalues_K{K_BEST}"
        
        adjust = assoc_configs[method]
        
        targets.append(f"{base}_sigSNPs_{adjust}.tsv")
        targets.append(f"{base}_{adjust}_genesAround{GENE_DISTANCE}.tsv")
        targets.append(f"{base}_{adjust}_genesAround{GENE_DISTANCE}_collapsed.tsv")
        
        targets += [f"plots/{method}/Rectangular-Manhattan.{trait}_K{K_BEST}_{adjust}.pdf" 
                   for trait in PREDICTORS_SELECTED.split(',')]
    
    return targets
#def parse_association_configs(config, varname = 'ASSOCIATION_CONFIGS'):
#	"""Parse association configurations from config file"""
#	assoc_configs = config.get(varname, [])
#
#	# If no configs specified, use legacy parameter
#	configs = defaultdict(list)
#	for cfg in assoc_configs:
#		method = cfg["METHOD"]
#		adjust = cfg["ADJUST"] + '_' + cfg["THRESHOLD"]
#		configs[method].append(adjust)
#
#	return configs
#
## Generate targets for all association configurations
#def get_association_targets(config):
#	
#	targets = []
#	assoc_configs = parse_association_configs(config)
#	
#	for method in assoc_configs:
#		if method == "EMMAX":
#			base = f"tables/EMMAX/EMMAX_pvalues_K{K_BEST}"
#		elif method == "LFMM":
#			base = f"tables/LFMM/LFMM_pvalues_K{K_BEST}"
#		#TODO add another method        	
#
#		# Add significant SNPs target
#		sigsnps = [ f"{base}_sigSNPs_{adjust}.tsv" for adjust in assoc_configs[method] ]
#		targets += sigsnps
#		
#		# Add genes around target
#		genes = [f"{base}_{adjust}_genesAround{GENE_DISTANCE}.tsv" for adjust in assoc_configs[method] ]
#		targets += genes
#
#		# Add genes around target
#		genes = [f"{base}_{adjust}_genesAround{GENE_DISTANCE}_collapsed.tsv" for adjust in assoc_configs[method] ]
#		targets += genes
#
#		# Add manhatan plots
#		targets += ["plots/" + method + "/" + f"Rectangular-Manhattan.{trait}_K{K_BEST}_{adjust}.pdf" for trait, adjust in product(PREDICTORS_SELECTED.split(','), assoc_configs[method]) ]
#	return targets

#def sig_snps_files(wildcards):
#    pattern = re.compile(
#        r"tables/(?P<method>[A-Za-z]+)/(?P=method)_pvalues_(?:PC|K)\d+_[^_]+_genesAround\.tsv"
#    )
#    files = []
#    for f in Path("tables").rglob("*_genesAround{GENE_DISTANCE}.tsv"):
#        if pattern.fullmatch(str(f)):
#            files.append(str(f))
#    if not files:
#        raise ValueError("No matching sigSNPs files found!")
#    return files
# [f"plots/Rectangular-Manhattan.{trait}_K{K_BEST}_{adjust}.pdf" for trait, adjust in product(PREDICTORS_SELECTED.split(','), ADJUST['EMMAX']) ]

######################### Basic arguments for all modes
CPU = config['CPU'] ; check_numeric(CPU, 'CPU')
INDIR = '/pipeline/' + config["INDIR"]
PROJECTNAME = config["PROJECTNAME"]
VCF_RAW = config['VCF_RAW'] ; check_vcf_name(INDIR, VCF_RAW, 'VCF_RAW') ; check_file_exist(INDIR, VCF_RAW, 'VCF_RAW')
FILENAME = VCF_RAW[:-4] if VCF_RAW.endswith('.vcf') else VCF_RAW[:-7] # extract basename
SAMPLES = config['SAMPLES'] ; check_file_exist(INDIR, SAMPLES, 'SAMPLES') # filename of table with columns site,sample,latitude,longitude

OUTDIR = '/pipeline/' + PROJECTNAME + '_results/'
LOGDIR = '/pipeline/' + PROJECTNAME + '_logs/'
# Create dirs
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(OUTDIR + 'tables', exist_ok=True)
os.makedirs(OUTDIR + 'plots', exist_ok=True)
os.makedirs(OUTDIR + 'intermediate', exist_ok=True)
workdir: OUTDIR

# Process VCF arguments
MAF = config['MAF'] ; check_float(MAF, 'MAF')
MISS = config['MISS'] ; check_float(MISS, 'MISS')
LDwin = config['LDwin'] ; check_numeric(LDwin, 'LDwin')
LDstep = config['LDstep'] ; check_numeric(LDstep, 'LDstep')
LDr2 = config['LDr2'] ; check_float(LDr2, 'LDr2')

# SNMF arguments
K_START = config['K_START'] ; check_numeric(K_START, 'K_START')
K_END = config['K_END'] ; check_numeric(K_END, 'K_END')
PLOIDY = config['PLOIDY']; check_numeric(PLOIDY, 'PLOIDY')
REPEAT = config['REPEAT']; check_numeric(REPEAT, 'REPEAT')
PROJECT = config['PROJECT'] # could be new or continue

# structure_K
K_BEST = config['K_BEST'] ; check_numeric(K_BEST, 'K_BEST', allowed_null = True) #TODO relising the range of K applied, in parallel imputing and GEA/GWAS, for comparing association results (it's very obvious in comparison of different K how well population structure was accounted (possibly sometimes required different K for different traits, for example precipitation in case of Israel highly linked with GEO, so there is some masking effect of bigger K)
#TODO also possibly different Ks for Impute,EMMAX,LFMM (LFMM works strangenly with low Ks, while for EMMAX from my experience lower K often better)
CUSTOM_TRAIT = INDIR + config['CUSTOM_TRAIT'] #TODO could be No or filename (check that file exists in INDIR directory
CROP_REGION = config['CROP_REGION'] ; check_comma_separated(CROP_REGION, 'CROP_REGION', numeric = True, allowed_single_values = ['auto', 'world'])
RESOLUTION = config['RESOLUTION'] ; check_in_list(RESOLUTION, [0.5, 2.5, 5, 10 ], 'RESOLUTION')
GAP = config['GAP'] ; check_float(GAP, 'GAP')# GAP for automatically identified map region
METRICS_WINSIZE = config['METRICS_WINSIZE']; check_numeric(METRICS_WINSIZE, 'METRICS_WINSIZE')
PieMap_PALLETE_MAP = config['PieMap_PALLETE_MAP']
PieMap_PALLETE_MAP_reverse = config['PieMap_PALLETE_MAP_reverse']
PieMap_PIE_SIZE = config['PieMap_PIE_SIZE']
PieMap_PIE_ALPHA = config['PieMap_PIE_ALPHA']
PieMap_IBD_ALPHA = config['PieMap_IBD_ALPHA']
PieMap_IBD_COLOR = config['PieMap_IBD_COLOR']
PieMap_POP_LABEL = config['PieMap_POP_LABEL'] ; check_in_list(PieMap_POP_LABEL, ['T', 'F'], 'PieMap_POP_LABEL')
PieMap_PIE_RESCALE = config['PieMap_PIE_RESCALE']; check_numeric(PieMap_PIE_RESCALE, 'PieMap_PIE_RESCALE')
PieMap_POP_LABEL_SIZE = config['PieMap_POP_LABEL_SIZE'] ; check_numeric(PieMap_POP_LABEL_SIZE, 'PieMap_POP_LABEL_SIZE')
# Association
#ADJUST = config['ADJUST'] #TODO add checking values bonf or qval - in all_input part
ADJUST = parse_association_configs(config) # dict where we can call EMMAX, LFMM or ANOTHER assoc method

print("DEBUG ADJUST:", ADJUST)
for method in ADJUST:
    print(f"DEBUG {method}:", ADJUST[method])

PREDICTORS_SELECTED = config['PREDICTORS_SELECTED']
GFF = config['GFF'] #TODO check function file exists
GFF_FEATURE = config['GFF_FEATURE'] #TODO add check function which check that value is one from the list of values(2 args)
GENE_DISTANCE = config['GENE_DISTANCE']
SNP_DISTANCE = config['SNP_DISTANCE']
PROMOTER_LENGTH = config['PROMOTER_LENGTH'] ; check_numeric(PROMOTER_LENGTH, 'PROMOTER_LENGTPROMOTER_LENGTHH')
sigSNPs_METHOD = config['sigSNPs_METHOD'] ; check_in_list(sigSNPs_METHOD, ['EMMAX', 'LFMM', 'Sum', 'Overlap', 'PairOverlap'], 'sigSNPs_METHOD')
sigSNPs_GAP = config['sigSNPs_GAP'] ; check_numeric(sigSNPs_GAP, 'sigSNPs_GAP')

# regionPlot
#REGIONPLOT_REQUEST = config['REGIONPLOT_REQUEST'] #TODO
GFF_GENENAME = config['GFF_GENENAME'] #TODO
GFF_BIOTYPE = config['GFF_BIOTYPE'] #TODO
REGIONPLOT_TRAITS = config['REGIONPLOT_TRAITS'] #TODO add cheking if there is that trait in metadata table
REGIONPLOT_TRAITS_SAFE = REGIONPLOT_TRAITS.replace(',', '_')

REGIONPLOT_ASSOCMETHOD = config['REGIONPLOT_ASSOCMETHOD'] ; check_in_list(REGIONPLOT_ASSOCMETHOD, ['EMMAX', 'LFMM'], REGIONPLOT_ASSOCMETHOD)
REGIONPLOT_REGION = config['REGIONPLOT_REGION'] # check is in script
REGIONPLOT_REGION_SAFE = config['REGIONPLOT_REGION'].replace(':', '_').replace('-', '_') # for saving in files
GENES_TO_HIGHLIGHT = config['GENES_TO_HIGHLIGHT']
#REGIONPLOT_REQUEST = config['REGIONPLOT_REQUEST'] 

# Maladaptation
SSP = config['SSP']
YEAR = config['YEAR']
MODELS = config['MODELS']
NTREE = config['NTREE']  # number of tree parameter
COR_THRESHOLD = config['COR_THRESHOLD'] # Correlation threshold parameter
PCNM = config['PCNM'] #TODO check values with or without
GF_SUFFIX = config['GF_SUFFIX']


# Rule to process all files
####################################
MODE = config.get('mode', 'Absence')
if MODE == 'structure':
	all_inputs = [f"tables/SNMF_clusters_K{K}.tsv" for K in range(K_START, K_END + 1, 1)] + ["plots/PCA.png"] # Add another plots
elif MODE == 'structure_K':	
	
	check_numeric(K_BEST, 'K_BEST') # null not allowed
	# Define required output file for mode
	all_inputs = [f"plots/PieMap_{trait}_TajimaD.png" for trait in PREDICTORS_SELECTED.split(',') ] + [f"plots/PieMap_{trait}_PiDiversity.png" for trait in PREDICTORS_SELECTED.split(',') ] + [f"plots/DensityPlot_{bio}.png" for bio in PREDICTORS_SELECTED.split(',') ] + ['plots/CorrelationHeatmap.png'] + ["plots/MantelTest.png"] + [f"plots/PieMap_{trait}_PiDiversity.svg" for trait in PREDICTORS_SELECTED.split(',')] + ["plots/AMOVA.png", "tables/AMOVA.tsv"]  
	# further should be PieMap plot
elif MODE == 'association':
	# Load specific arguments
	#TODO ADJUST add checking values bonf or qval
	check_numeric(K_BEST, 'K_BEST') # null not allowed
	check_comma_separated(PREDICTORS_SELECTED, 'PREDICTORS_SELECTED')
	#TODO GFF check function file exists
	#TODO GFF_FEATURE add check function which check that value is one from the list of values(2 args)
	check_numeric(GENE_DISTANCE, 'GENE_DISTANCE')
	check_numeric(SNP_DISTANCE, 'SNP_DISTANCE')

	# Create output directories for manhatan plots
	for method in ADJUST: # iterate through keys
		os.makedirs(OUTDIR + 'plots/' + method, exist_ok=True)
		os.makedirs(OUTDIR + 'tables/' + method, exist_ok=True)
		

	# Define required output file for mode
	all_inputs = get_association_targets(ADJUST, config) + [f"tables/Selected_SNPs_redundant_{sigSNPs_METHOD}_Gap{sigSNPs_GAP}.tsv", f"tables/Selected_SNPs_unique_{sigSNPs_METHOD}_Gap{sigSNPs_GAP}.list" ] #, f"tables/Selected_SNPs_redundant_{sigSNPs_METHOD}_Gap{sigSNPs_GAP}_genesAround{GENE_DISTANCE}.tsv"] #TODO removed for testing on trifolium
	print(all_inputs) #TODO TEMP

elif MODE == 'regionplot':
	all_inputs = ['topr_gene_annotation.tsv',
f"plots/regionPlot_{REGIONPLOT_TRAITS_SAFE}_{REGIONPLOT_REGION_SAFE}_{ADJUST[REGIONPLOT_ASSOCMETHOD]}.png"	
]

elif MODE == 'maladaptation':
	check_numeric(SSP, 'SSP') #TODO add checking that it in list
	check_numeric(NTREE, 'NTREE')
	check_float(COR_THRESHOLD, 'COR_THRESHOLD')
	#TODO add checking YEAR in list
	#TODO add checking MODELS that all elements in list
	
	all_inputs = [
			#TODO adjust final files, make + not new list if another mode selected
		      #get_association_targets(ADJUST, config),
		      f"tables/Climate_future_year{YEAR}_ssp{SSP}_site.tsv",
		      f"tables/gradientForest_GeneticOffesetValues_{GF_SUFFIX}_{PCNM}PCNM.tsv",
		      f"plots/gradientForest_GeneticOffesetPieMap_{GF_SUFFIX}_{PCNM}PCNM.png"] #TODO maybe add MODELS with some formating in the name of the file?
	
	
elif MODE == 'Absence':
	raise ValueError(f"No mode selected. Please provide any with --config mode=<MODE> argument")
else:
	raise ValueError(f"Unknown mode: {mode}. Please provide: --config mode=structure or mode=association")


##################################### structure
rule all:
	input:
		all_inputs

rule extract_samples:
	input:
		samples=f"{INDIR}{SAMPLES}"
	output:
		samples_list = "samples_metadata.list"
	threads: 1
	log:
		f"{LOGDIR}extract_samples.log"
	shell:
		"(tail -n +2 {input.samples} | awk '{{print 0,$2}}'  > {output.samples_list}) > {log} 2>&1"
	
rule process_vcf:
	input: 
		vcf=lambda wildcards: f"{INDIR}{FILENAME}.vcf.gz" if VCF_RAW.endswith('.vcf.gz') else f"{INDIR}{FILENAME}.vcf",
		samples="samples_metadata.list"
	output:
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.vcf",
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcf",
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.eigenvec"
	threads: CPU
	log:
		f"{LOGDIR}process_vcf.log"
	shell:
		"/pipeline/scripts/process_vcf.sh {input.vcf} {input.samples} {MAF} {MISS} {LDwin} {LDstep} {LDr2} > {log} 2>&1"

rule extract_samples_order_vcf:
	input:
		vcf = f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcf"
	output:
		"samples_vcf.list"
	threads: 1
	log:
		f"{LOGDIR}extract_samples_order_vcf.log"
	shell:
		"grep -m1 CHROM {input.vcf} | cut -f10- | tr '\t' '\n' > samples_vcf.list "

rule filter_arrange_metadata:
	input:
		metadata = f"{INDIR}{SAMPLES}",
		samples_vcf = "samples_vcf.list"
	output:
		"Metadata_filtered_arranged.tsv"
	threads: 1
	log:
		f"{LOGDIR}filter_arrange_metadata.log"
	shell:
		"Rscript /pipeline/scripts/filter_arrange_metadata.R {input.metadata} {input.samples_vcf} {output} > {log} 2>&1"	

# Convert only LD for PopStructure inference
rule vcf2lfmm_LD:
	input:
		vcf=f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.vcf"
	output:
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.geno",
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.lfmm",
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.vcfsnp", #TODO maybe it's already removed so simplify the rule in gradForest
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.removed"
	threads: 1
	log:
		f"{LOGDIR}vcf2lfmm_LD.log"
	shell:
		"Rscript /pipeline/scripts/vcf2lfmm.R {input.vcf} > {log} 2>&1"

# Run SNMF structure pipeline
rule SNMF_LD:
  input:
          geno=f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.geno"
  output:
          f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.snmfProject"
  threads: CPU
  log:
          f"{LOGDIR}SNMF_LD.log"
  shell:
          "Rscript /pipeline/scripts/SNMF.R {input.geno} {K_START} {K_END} {PLOIDY} {REPEAT} {threads} {PROJECT} > {log} 2>&1"

rule Impute_LD: 
	input:
		snmf = f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.snmfProject",
		lfmm = f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.lfmm"
	output:
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}_K{K_BEST}imputed.lfmm"
	log:
		f"{LOGDIR}Impute_LD.log"
	run:
		if config.get("run_Impute_LD", True): #TEMP
			shell(f"Rscript /pipeline/scripts/Impute.R {input.snmf} {input.lfmm} {K_BEST} > {log} 2>&1")
		else:
			print("Skipping Impute_LD rule as per config")

rule lfmm2vcf_LD:
	input:
		vcf = f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.vcf",
		lfmm_imp = f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}_K{K_BEST}imputed.lfmm"
	output:
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}_K{K_BEST}imputed.vcf"
	params:
		blocksize = 20000 # Blocksize for transposeBigData, set as deafult but could be changed if face some resource problems
	log:
		f"{LOGDIR}lfmm2vcf_LD.log"
	shell:
		"Rscript /pipeline/scripts/lfmm2vcf.R {input.lfmm_imp} {input.vcf} 20000 > {log} 2>&1" 
rule SNMF_LD_plots:
	input:
		snmfProject=f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.snmfProject",
		lfmm=f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.lfmm",
		samples="Metadata_filtered_arranged.tsv"
	output:
		[f"tables/SNMF_clusters_K{K}.tsv" for K in range(K_START, K_END + 1, 1)],
		'plots/PCA.png'
	threads: 1
	log:
		f"{LOGDIR}SNMF_LD_plots.log"
	shell:
		"Rscript /pipeline/scripts/SNMF_plots.R {input.snmfProject} {input.lfmm} {K_START} {K_END} {PLOIDY} {input.samples} > {log} 2>&1"

rule AMOVA: # for now it subset 10k SNPs from LD prunned file so it fast return approximate AMOVA results (too long with sufficient number of SNPs)
	input:
		vcf = f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}_K{K_BEST}imputed.vcf", # use imputed data for more stable performance (failed with miss = 0.2)
		sample = "Metadata_filtered_arranged.tsv" 
	output:
		"plots/AMOVA.png",
		"tables/AMOVA.tsv"
	threads: 1
	log:
		f"{LOGDIR}AMOVA.log"
	shell:
		"Rscript /pipeline/scripts/AMOVA.R {input.vcf} {input.sample} {threads} > {log} 2>&1"

##################### structure_K required Kbest

# Convert whole vcf to lfmm
rule vcf2lfmm:
	input:
		vcf=f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcf"
	output:
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.geno",
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.lfmm",
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcfsnp",
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.removed"
	log:
		f"{LOGDIR}vcf2lfmm.log"
	shell:
		"Rscript /pipeline/scripts/vcf2lfmm.R {input.vcf} > {log} 2>&1"

# run SNMF for whole dataset for imputation missing values
rule SNMF:
	input:
		geno=f"{FILENAME}_MAF{MAF}_MISS{MISS}.geno" # after correctly define structure (if not structure flag with trun)
	output: 
		f"{FILENAME}_MAF{MAF}_MISS{MISS}.snmfProject"
	log:
		f"{LOGDIR}SNMF.log"
	shell:
		"Rscript /pipeline/scripts/SNMF.R {input.geno} {K_START} {K_END} {PLOIDY} {REPEAT} {CPU} {PROJECT} > {log} 2>&1"
# 

# After Imputation we can calculate neutrality test and IBD (maybe it's better to calculate it on missing values data)
rule TajimaD:	#TODO add exception for only 1 individual in population
	input: 
		vcf=f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcf",
		samples= "Metadata_filtered_arranged.tsv"
	params:
		winsize = METRICS_WINSIZE		
	output:
		"tables/TajimaD_TotalByPop.tsv"
	log:
		f"{LOGDIR}TajimaD.log"
	threads: CPU
	shell:
		"Rscript /pipeline/scripts/TajimaD.R {input.vcf} {input.samples} {params.winsize} {CPU} > {log} 2>&1"

rule Pi_diversity: #TODO add exception for only 1 individual in population
	input:
		vcf=f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcf",
		samples= "Metadata_filtered_arranged.tsv"
	params:
		winsize = METRICS_WINSIZE
	output:
		"tables/Pi_diversity_TotalByPop.tsv"
	log:
		f"{LOGDIR}Pi_diversity.log"
	threads: CPU
	shell:
		"Rscript /pipeline/scripts/Pi_diversity.R {input.vcf} {input.samples} {params.winsize} {CPU} > {log} 2>&1"
		
rule IBD: #TODO not very reprodicible, CHECK IT!
	input:
		cluster=f"tables/SNMF_clusters_K{K_BEST}.tsv",
		samples="Metadata_filtered_arranged.tsv"
	output:
		"tables/IBD_raw.tsv",
		"tables/IBD_notIsolated.tsv" # for plots
	log:
		f"{LOGDIR}IBD.log"
	shell:
		"Rscript /pipeline/scripts/IBD.R {input.cluster} {input.samples} {CPU} > {log} 2>&1"

rule download_climate_present:
	input:
		lfmm=f"{FILENAME}_MAF{MAF}_MISS{MISS}.lfmm",
		samples="Metadata_filtered_arranged.tsv"
	params:
		resolution=RESOLUTION
	output:
		"tables/Climate_present_region.tsv",
		"tables/Climate_present_region_scaled.tsv",
		"tables/Climate_present_all.tsv",
		"intermediate/Climate_present_RasterStack.grd",
		[f"plots/DensityPlot_{bio}.png" for bio in PREDICTORS_SELECTED.split(',') ]
	log:
		f"{LOGDIR}download_climate_present.log"
	shell:
		"Rscript /pipeline/scripts/Download_PresentClimate.R {input.samples} {CROP_REGION} {GAP} /pipeline/data/ {params.resolution} > {log} 2>&1" #TODO manual climate dir

rule pie_map_plot:
	input:
		samples="Metadata_filtered_arranged.tsv",
		clusters=f"tables/SNMF_clusters_K{K_BEST}.tsv",
		raster=f"intermediate/Climate_present_RasterStack.grd",
		ibd="tables/IBD_notIsolated.tsv",
		tajima="tables/TajimaD_TotalByPop.tsv",
		diversity="tables/Pi_diversity_TotalByPop.tsv"
	params:
		custom_trait=CUSTOM_TRAIT,
		palette_map=PieMap_PALLETE_MAP,
		palette_map_reverse=PieMap_PALLETE_MAP_reverse,
		pie_size=PieMap_PIE_SIZE,
		pie_alpha=PieMap_PIE_ALPHA,
		ibd_alpha=PieMap_IBD_ALPHA,
		ibd_color=PieMap_IBD_COLOR,
		pop_label=PieMap_POP_LABEL,
		pop_label_size=PieMap_POP_LABEL_SIZE,
		pie_rescale=PieMap_PIE_RESCALE
	output:
		[f"plots/PieMap_{trait}_TajimaD.png" for trait in PREDICTORS_SELECTED.split(',')],
		[f"plots/PieMap_{trait}_TajimaD.svg" for trait in PREDICTORS_SELECTED.split(',')],
		[f"plots/PieMap_{trait}_PiDiversity.png" for trait in PREDICTORS_SELECTED.split(',')],
		[f"plots/PieMap_{trait}_PiDiversity.svg" for trait in PREDICTORS_SELECTED.split(',')]
	log:
		f"{LOGDIR}pie_map_plot.log"
	shell:
		"Rscript /pipeline/scripts/PieMap_plot.R {input.samples} {input.clusters} {input.raster} {input.ibd} {input.tajima} {input.diversity} {params.custom_trait} {params.palette_map} {params.palette_map_reverse} {params.pie_size} {params.pie_alpha} {params.ibd_alpha} {params.ibd_color} {params.pop_label} {params.pie_rescale} {params.pop_label_size} > {log} 2>&1"

rule Plot_CorrelationHM:
	input:
		samples = "Metadata_filtered_arranged.tsv",
		climate = "tables/Climate_present_region.tsv"
	output:
		'plots/CorrelationHeatmap.png'
	log:
		f"{LOGDIR}Plot_CorrelationHM.log"
	shell:
		"Rscript /pipeline/scripts/Plot_CorrelationHeatmap.R {input.climate} {input.samples} > {log} 2>&1"


#TODO one-tail behavior of vegan::mantel - change on another one mantel.testr()
#TODO handle scaling?
rule MantelTest:
	input:
		samples ="Metadata_filtered_arranged.tsv", 
		climate = "tables/Climate_present_region.tsv",
		clusters= f"tables/SNMF_clusters_K{K_BEST}.tsv"
	output:
		"plots/MantelTest.png"
	log:
		f"{LOGDIR}MantelTest.log"
	shell:
		"Rscript /pipeline/scripts/MantelTest.R {input.samples} {input.clusters} {input.climate} {PREDICTORS_SELECTED} > {log} 2>&1"
################################ association mode ################################
# Impute also LD dataset for LFMM model training 
if 'LFMM' in ADJUST:

	rule Impute:
		input:
			snmf = f"{FILENAME}_MAF{MAF}_MISS{MISS}.snmfProject",
			lfmm = f"{FILENAME}_MAF{MAF}_MISS{MISS}.lfmm"
		output:
			f"{FILENAME}_MAF{MAF}_MISS{MISS}_K{K_BEST}imputed.lfmm"
		log:
			f"{LOGDIR}Impute.log"
		run:
			if config.get('run_Impute', True): #TEMP
				shell(f"Rscript /pipeline/scripts/Impute.R {input.snmf} {input.lfmm} {K_BEST} > {log} 2>&1")
			else:
				print("Skipping Impute rule as per config")
	#		"Rscript /pipeline/scripts/Impute.R {input.snmf} {input.lfmm} {K_BEST} > {log} 2>&1"
	rule lfmm2vcf:
		input:
			vcf = f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcf",
			lfmm_imp = f"{FILENAME}_MAF{MAF}_MISS{MISS}_K{K_BEST}imputed.lfmm"
		output:
			f"{FILENAME}_MAF{MAF}_MISS{MISS}_K{K_BEST}imputed.vcf"
		log:
			f"{LOGDIR}lfmm2vcf.log"
		shell:
			"Rscript /pipeline/scripts/lfmm2vcf.R {input.lfmm_imp} {input.vcf} 20000 > {log} 2>&1" # Blocksize for transposeBigData, set as deafult but could be changed if face some resource problems
	
	rule LFMM:
		input: 
			lfmm = f"{FILENAME}_MAF{MAF}_MISS{MISS}_K{K_BEST}imputed.lfmm",
			lfmm_LD = f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}_K{K_BEST}imputed.lfmm",
			climate="tables/Climate_present_region_scaled.tsv",
			vcfsnp=f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcfsnp"
		output:
			f"tables/LFMM/LFMM_pvalues_K{K_BEST}.tsv"
		log:
			f"{LOGDIR}LFMM.log"
		shell:
			"Rscript /pipeline/scripts/LFMM.R {input.lfmm_LD} {input.lfmm} {input.climate} {K_BEST} {PREDICTORS_SELECTED} {input.vcfsnp} > {log} 2>&1"

	rule LFMM_Manhattan_plot:
		input:
			assoc = f"tables/LFMM/LFMM_pvalues_K{K_BEST}.tsv"
		output:
			[f"plots/LFMM/Rectangular-Manhattan.{trait}_K{K_BEST}_{ADJUST['LFMM']}.pdf" for trait in PREDICTORS_SELECTED.split(',') ]
		params:
			K = K_BEST,
			adjust = ADJUST['LFMM'],
			method = 'LFMM' # 'plots/LFMM'
		log:
			f"{LOGDIR}LFMM_manhattan_plots.log"
		shell:
			"Rscript /pipeline/scripts/Plot_manhat_CMplot.R {input.assoc} {params.adjust} {params.K} {params.method} > {log} 2>&1"

	rule LFMM_sigSNP:
		input:
			assoc = f"tables/LFMM/LFMM_pvalues_K{K_BEST}.tsv"
		output:
			f"tables/LFMM/LFMM_pvalues_K{K_BEST}_sigSNPs_{ADJUST['LFMM']}.tsv"
		params:
			adjust = ADJUST['LFMM'],
			snp_dist = SNP_DISTANCE,
			cpu = CPU,
			method = 'LFMM'
		log:
			f"{LOGDIR}LFMM_sigSNP.log"
		shell:
			"Rscript /pipeline/scripts/find_sigSNPs.R {input.assoc} {params.adjust} {params.snp_dist} {params.method} {params.cpu} > {log} 2>&1"

	rule LFMM_sigSNP_genesAround:
		input:
			sigsnps = f"tables/LFMM/LFMM_pvalues_K{K_BEST}_sigSNPs_{ADJUST['LFMM']}.tsv",
			gff = f"{INDIR}{GFF}",
			vcfsnp = f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcfsnp"
		output:
			f"tables/LFMM/LFMM_pvalues_K{K_BEST}_{ADJUST['LFMM']}_genesAround{GENE_DISTANCE}.tsv",
			f"tables/LFMM/LFMM_pvalues_K{K_BEST}_{ADJUST['LFMM']}_genesAround{GENE_DISTANCE}_collapsed.tsv" 
		params:
			feature = GFF_FEATURE,
			distance = GENE_DISTANCE,
			promoter_len = PROMOTER_LENGTH
		log:
			f"{LOGDIR}EMMAX_sigSNP_genesAround.log" #TODO maybe add DISTANCE to name of 
		shell:
			"Rscript /pipeline/scripts/find_genes_around_sigSNPs.R {input.gff} {input.sigsnps} {params.feature} {params.distance} {params.promoter_len} {input.vcfsnp} {CPU} > {log} 2>&1"

if 'EMMAX' in ADJUST:
	rule EMMAX:
		input:
			vcf = f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcf",
			traits = "tables/Climate_present_region_scaled.tsv",
			covariates = f"{FILENAME}_MAF{MAF}_MISS{MISS}.LD{LDr2}.eigenvec",
			metadata ="Metadata_filtered_arranged.tsv" #TODO for now it's unscalled phen traits
		output:
			f"tables/EMMAX/EMMAX_pvalues_K{K_BEST}.tsv"
		log:
			f"{LOGDIR}EMMAX.log"
		shell:
			"Rscript /pipeline/scripts/EMMAX.R {input.vcf} {K_BEST} {input.traits} {input.covariates} {PREDICTORS_SELECTED} intermediate/ {input.metadata} > {log} 2>&1"

	rule EMMAX_Manhattan_plot:
		input:
			assoc = f"tables/EMMAX/EMMAX_pvalues_K{K_BEST}.tsv"
		output:
			[f"plots/EMMAX/Rectangular-Manhattan.{trait}_K{K_BEST}_{ADJUST['EMMAX']}.pdf" for trait in PREDICTORS_SELECTED.split(',')]
		params:
			K = K_BEST,
			adjust = ADJUST['EMMAX'],
			method = 'EMMAX'

		log:
			f"{LOGDIR}EMMAX_manhattan_plots.log"
		shell:
			"Rscript /pipeline/scripts/Plot_manhat_CMplot.R {input.assoc} {params.adjust} {params.K} {params.method} > {log} 2>&1"

	#TODO make it universal not only for EMMAX (add ASSOC_METHOD variable)
	rule EMMAX_sigSNP:
		input:
			assoc = f"tables/EMMAX/EMMAX_pvalues_K{K_BEST}.tsv"
		output:
			f"tables/EMMAX/EMMAX_pvalues_K{K_BEST}_sigSNPs_{ADJUST['EMMAX']}.tsv"
		params:
			adjust = ADJUST['EMMAX'],
			snp_dist = SNP_DISTANCE,
			cpu = CPU,
			method = 'EMMAX'
		log:
			f"{LOGDIR}EMMAX_sigSNP.log" 
		shell:
			"Rscript /pipeline/scripts/find_sigSNPs.R {input.assoc} {params.adjust} {params.snp_dist} {params.method} {params.cpu} > {log} 2>&1"

	rule EMMAX_sigSNP_genesAround:
		input:
			sigsnps = f"tables/EMMAX/EMMAX_pvalues_K{K_BEST}_sigSNPs_{ADJUST['EMMAX']}.tsv",
			gff = f"{INDIR}{GFF}",
			vcfsnp = f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcfsnp"
		output:
			f"tables/EMMAX/EMMAX_pvalues_K{K_BEST}_{ADJUST['EMMAX']}_genesAround{GENE_DISTANCE}.tsv",
			f"tables/EMMAX/EMMAX_pvalues_K{K_BEST}_{ADJUST['EMMAX']}_genesAround{GENE_DISTANCE}_collapsed.tsv"
		params:
			feature = GFF_FEATURE,
			distance = GENE_DISTANCE,
			promoter_len = PROMOTER_LENGTH
		log:
			f"{LOGDIR}EMMAX_sigSNP_genesAround.log" #TODO maybe add DISTANCE to name of file
		shell:
			"Rscript /pipeline/scripts/find_genes_around_sigSNPs.R {input.gff} {input.sigsnps} {params.feature} {params.distance} {params.promoter_len} {input.vcfsnp} {CPU} > {log} 2>&1"

# Summary of association analysis
#TODO make some summing\comparing plots 
#TODO circular plot with only significant to see overlaps (think about tool for vizualization (probably interactive)
rule find_sigSNP_overlap: # between EMMAX and LFMM and #TODO for another method (RDA?)
	input: # will return error if some of them doesn't exist
	  sigsnps_lst = [f"tables/{method}/{method}_pvalues_K{K_BEST}_sigSNPs_{ADJUST[method]}.tsv" for method in ADJUST]
	params:
	  sigsnps_str = lambda wildcards, input: ','.join(input.sigsnps_lst),
	  method = sigSNPs_METHOD,
	  gap = sigSNPs_GAP,
	  predictors_selected = PREDICTORS_SELECTED
	output:
		f"tables/Selected_SNPs_redundant_{sigSNPs_METHOD}_Gap{sigSNPs_GAP}.tsv",
		f"tables/Selected_SNPs_unique_{sigSNPs_METHOD}_Gap{sigSNPs_GAP}.list"
	log:
		f"{LOGDIR}find_sigSNP_overlap.log"
	shell:
		"Rscript /pipeline/scripts/find_sigSNPs_overlap.R {params.sigsnps_str} {params.method} {params.gap} {params.predictors_selected} > {log} 2>&1"

rule Overlap_sigSNP_genesAround:
	input:
		sigsnps = f"tables/Selected_SNPs_redundant_{sigSNPs_METHOD}_Gap{sigSNPs_GAP}.tsv",
		gff = f"{INDIR}{GFF}",
		vcfsnp = f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcfsnp"
	params:
		feature = GFF_FEATURE,
		distance = GENE_DISTANCE,
		promoter_len = PROMOTER_LENGTH
	output:
		f"tables/Selected_SNPs_redundant_{sigSNPs_METHOD}_Gap{sigSNPs_GAP}_genesAround{GENE_DISTANCE}.tsv",
		f"tables/Selected_SNPs_redundant_{sigSNPs_METHOD}_Gap{sigSNPs_GAP}_genesAround{GENE_DISTANCE}_collapsed.tsv"
	log:
		f"{LOGDIR}Overlap_sigSNP_genesAround.log" #TODO maybe add DISTANCE to name of file
	shell:
		"Rscript /pipeline/scripts/find_genes_around_sigSNPs.R {input.gff} {input.sigsnps} {params.feature} {params.distance} {params.promoter_len} {input.vcfsnp} {CPU} > {log} 2>&1"
		
		 
##################### regionalPlot #########################
#TODO make some variable to choose what to print, by default some top 10 regions? make available to drow subset of SNPs provided 

rule gff2topr:
	input:
		gff = f"{INDIR}{GFF}"
	params:
		feature=GFF_FEATURE,
		genename=GFF_GENENAME,
		biotype=GFF_BIOTYPE
	output:
		outfile='topr_gene_annotation.tsv'
	log:
		f"{LOGDIR}gff2topr.log"
	shell:
		"python3 /pipeline/scripts/gff2topr.py {input.gff} {params.feature} {params.genename} {params.biotype} {output.outfile} > {log} 2>&1" 

#TODO add option for climate factors EMMAX/LFMM to build regional plot, probably joined 
#TODO for now only EMMAX and only 1 GWAS, rework it later
rule regionPlot:
	input:
		assoc = f"tables/{REGIONPLOT_ASSOCMETHOD}/{REGIONPLOT_ASSOCMETHOD}_pvalues_K{K_BEST}.tsv",
		gff_topr = 'topr_gene_annotation.tsv'
	params:
		region = REGIONPLOT_REGION,
		adjust = ADJUST[REGIONPLOT_ASSOCMETHOD],
		traits = REGIONPLOT_TRAITS,
		genes_to_highlight = GENES_TO_HIGHLIGHT
	output:
		f"plots/regionPlot_{REGIONPLOT_TRAITS_SAFE}_{REGIONPLOT_REGION_SAFE}_{ADJUST[REGIONPLOT_ASSOCMETHOD]}.png"	
	log:
		f"{LOGDIR}regionPlot.log"
	shell:
		"Rscript /pipeline/scripts/Plot_regionManh.R {input.assoc} {input.gff_topr} {params.region} {params.adjust} {params.traits} {params.genes_to_highlight} > {log} 2>&1"
##################### Maladaptation  ####################### part4
rule download_climate_future:
	input:
		lfmm=f"{FILENAME}_MAF{MAF}_MISS{MISS}.lfmm", # when starts downloading, can keep it like this
		samples="Metadata_filtered_arranged.tsv"
	params:
		resolution = RESOLUTION
	output:
		f"tables/Climate_future_year{YEAR}_ssp{SSP}_site.tsv",
		f"tables/Climate_future_year{YEAR}_ssp{SSP}_all.tsv",
		f"intermediate/Climate_future_year{YEAR}_ssp{SSP}_RasterStack.grd"
	log:
		f"{LOGDIR}download_climate_future_year{YEAR}_ssp{SSP}.log"
	shell:
		"Rscript /pipeline/scripts/Download_FutureClimate.R {input.samples} {CROP_REGION} {GAP} {INDIR} {SSP} {YEAR} {MODELS} {CPU} {params.resolution} > {log} 2>&1"

#TODO change random subsetting add while tryCatch or smth like this and advice to increase number of SNPs
rule gradientForest:
	input:
		lfmm = f"{FILENAME}_MAF{MAF}_MISS{MISS}.lfmm", #TODO add choice to make GF on imputed (if confidence) or notImp data
		sigsnps = f"tables/Selected_SNPs_unique_{sigSNPs_METHOD}_Gap{sigSNPs_GAP}.list",
		vcfsnp = f"{FILENAME}_MAF{MAF}_MISS{MISS}.vcfsnp",
		removed = f"{FILENAME}_MAF{MAF}_MISS{MISS}.removed",
		samples ="Metadata_filtered_arranged.tsv", 
		climate = "tables/Climate_present_region.tsv" # Here we need not scaled climate values for better interpretability
	output:
		f'intermediate/gradientForest_random_{GF_SUFFIX}_{PCNM}PCNM.qs', #PC3_sigSNPs_bonf_0.05
		f'intermediate/gradientForest_adaptive_{GF_SUFFIX}_{PCNM}PCNM.qs'
	log:
		f"{LOGDIR}gradientForest.log"
	shell:
		"Rscript /pipeline/scripts/gradientForest.R {input.lfmm} {input.sigsnps} {input.vcfsnp} {input.removed} {input.samples} {input.climate} {PREDICTORS_SELECTED} {NTREE} {COR_THRESHOLD} {PCNM} {GF_SUFFIX} > {log} 2>&1"
#Selected_SNPs_unique_EMMAX_Gap10000.list

rule gradientForest_plot:
	input:
		gf=f'intermediate/gradientForest_adaptive_{GF_SUFFIX}_{PCNM}PCNM.qs',
		gf_random=f'intermediate/gradientForest_random_{GF_SUFFIX}_{PCNM}PCNM.qs',
		future_climate_all=f"tables/Climate_future_year{YEAR}_ssp{SSP}_all.tsv",
		present_climate_all=f"tables/Climate_present_all.tsv",
		raster_present="intermediate/Climate_present_RasterStack.grd",
		clusters=f"tables/SNMF_clusters_K{K_BEST}.tsv",
		ibd="tables/IBD_notIsolated.tsv",
		tajima="tables/TajimaD_TotalByPop.tsv",
		samples="Metadata_filtered_arranged.tsv"
	output:
		f"plots/gradientForest_GeneticOffesetPieMap_{GF_SUFFIX}_{PCNM}PCNM.png",
		f"tables/gradientForest_GeneticOffesetValues_{GF_SUFFIX}_{PCNM}PCNM.tsv"
	log:
		f"{LOGDIR}gradientForest_plot.log"
	shell:
		"Rscript /pipeline/scripts/gradientForest_plots.R {input.gf} {input.gf_random} {PREDICTORS_SELECTED} {input.future_climate_all} {input.present_climate_all} {input.raster_present} {input.samples} {input.clusters} {input.ibd} {input.tajima} {PieMap_PALLETE_MAP} {PieMap_PIE_SIZE} {PieMap_PIE_ALPHA} {PieMap_IBD_ALPHA} {PieMap_IBD_COLOR} {PCNM} {GF_SUFFIX} > {log} 2>&1"	

