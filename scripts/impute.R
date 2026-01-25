library(LEA)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#################################
SNMF = args[1]    # path to snmfProject
LFMM = args[2]    # path to LFMM with missing values
K = args[3] %>% as.numeric
OUTPUT = args[4]  # explicit output path
#################################

# Load project
project <- load.snmfProject(SNMF)

# Select the run with the lowest cross-entropy value
best <- which.min(cross.entropy(project, K = K))

# Impute the missing genotypes
# LEA creates output as {LFMM}_imputed.lfmm in same directory
impute(project, LFMM, method = 'mode', K = K, run = best)

# LEA outputs to {LFMM}_imputed.lfmm, rename to desired output
lea_output <- paste0(LFMM, '_imputed.lfmm')
file.rename(from = lea_output, to = OUTPUT)
