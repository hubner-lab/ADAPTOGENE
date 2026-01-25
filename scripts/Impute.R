
library(LEA)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
######################
SNMF = args[1]   # path to snmfProject
LFMM = args[2]   # path to LFMM with missing values
Kbest = args[3] %>% as.numeric
OUTPUT = args[4] # output path for imputed lfmm file

######################
project = load.snmfProject(SNMF)
# Select the run with the lowest cross-entropy value
Entropy.best = which.min(cross.entropy(project, K = Kbest))
# Impute the missing genotypes
impute(project, LFMM, method = 'mode', K = Kbest, run = Entropy.best)
# Rename to expected output path
file.rename(from = paste0(LFMM, '_imputed.lfmm'), to = OUTPUT)

