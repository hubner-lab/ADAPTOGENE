library(dplyr)
library(data.table)
library(parallel)
args = commandArgs(trailingOnly=TRUE)
#################################
VCF = args[1]
SAMPLES = args[2]
WINSIZE = args[3]
CPU = args[4]
OUTPUT = args[5]
INTER_DIR = args[6]
INTERMEDIATE = paste0(INTER_DIR, '/')
VCFTOOLS = '/app/vcftools-0.1.16/bin/vcftools'
#################################

################### Functions ###################
FUN_pi <- function(LFMM_vcf, samples, INTERMEDIATE, VCFTOOLS, WINSIZE, CPU){

  mclapply(samples$site %>% unique, function(site){
	COMMAND = paste(VCFTOOLS, '--vcf', LFMM_vcf,
                   '--window-pi', WINSIZE, 
                   paste('--indv', samples$sample[samples$site == site], collapse = ' '),
                   '--out', paste0(INTERMEDIATE, 'pop', site)
                   )
							 	message(COMMAND)
	system(COMMAND)
    }, mc.cores = CPU)

  FILENAMES = paste0(INTERMEDIATE, 'pop', samples$site %>% unique, '.windowed.pi')
  pi_diversity <-
    lapply(FILENAMES, function(f) {
      popname = gsub(paste0(INTERMEDIATE, 'pop'), '', f) %>% 
	      	gsub('.windowed.pi', '', .) # extract popname
      fread(f, header = T) %>%
        dplyr::mutate(site = as.factor(popname)) %>%
	dplyr::group_by(site) %>%
        dplyr::summarise(PI = mean(PI, na.rm = T))
      }) %>% do.call(rbind, .)

  return(pi_diversity)
}

################## Main ##########################
samples <- fread(SAMPLES, 
		 colClasses = c("site" = "character", 'sample' = 'character'))
							message(samples %>% str)
							message("INFO: RUN DIVERSITY TEST FOR WHOLE GENOME")
FUN_pi(VCF, samples, INTERMEDIATE, VCFTOOLS, WINSIZE, CPU) %>%
	fwrite(OUTPUT, sep = '\t')
