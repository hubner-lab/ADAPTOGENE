library(LEA)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
set.seed(42)

####################
GENO=args[1] ; LFMM=sub('\\.geno$', '.lfmm', GENO)
K_START=args[2] %>% as.numeric
K_END=args[3] %>% as.numeric
PLOIDY=args[4] %>% as.numeric
REPEAT=args[5] %>% as.numeric
CPU=args[6] %>% as.numeric
PROJECT=args[7] # new or continue
####################

# FUN SNMF analysis
FUN_snmf <- function(Ks,
                     Ke,
                     GENO,
		     ploidy = 2,
                     repetions = 10,
		     entropy = T,
                     project = 'new', # with force run TRUE
                     I = 10000,
                     CPU = 24){

  SNMF = paste0(gsub('.geno', '', GENO), '.snmfProject')

  

   project = snmf(GENO,
                    CPU = CPU,
                    K = Ks:Ke,
                    entropy = TRUE,
                    repetitions = repetions,
                    ploidy = ploidy,
                    I = I,
                    project = project)                    

  return(project)
}


################### MAIN #########################

# Define the number of SNPs in GENO
nSNP = fread(GENO) %>% nrow
# Run SNMF
FUN_snmf(K_START,
	 K_END,
	 GENO,
	 ploidy = PLOIDY,
	 repetions = REPEAT,
	 project = PROJECT,
	 I = min(10000, nSNP),
	 CPU = CPU)
	
