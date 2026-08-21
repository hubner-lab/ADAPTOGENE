
library(LEA)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
######################
SNMF = args[1]      # path to snmfProject
LFMM = args[2]      # path to LFMM with missing values
Kbest = args[3] %>% as.numeric
OUTPUT = args[4]    # output path for imputed lfmm file
NMISSING = args[5]  # path to the has-missing flag written by vcf2lfmm.R ("0" or "1")

######################
has_missing <- as.logical(as.integer(readLines(NMISSING, n = 1)))

if (!has_missing) {
    # LEA::impute() requires at least one missing (9-coded) genotype and errors
    # otherwise -- a genuinely complete dataset (e.g. SLiM simulation output) has
    # nothing to impute, so just pass the LFMM through unchanged.
    message('INFO: No missing genotypes (see ', NMISSING, ') -- skipping LEA::impute(), copying LFMM as-is')
    file.copy(from = LFMM, to = OUTPUT, overwrite = TRUE)
} else {
    project = load.snmfProject(SNMF)
    # Select the run with the lowest cross-entropy value
    Entropy.best = which.min(cross.entropy(project, K = Kbest))
    # Impute the missing genotypes
    impute(project, LFMM, method = 'mode', K = Kbest, run = Entropy.best)
    # Rename to expected output path
    file.rename(from = paste0(LFMM, '_imputed.lfmm'), to = OUTPUT)
}

