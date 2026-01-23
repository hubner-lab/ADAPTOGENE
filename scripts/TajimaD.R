library(dplyr)
library(data.table)
library(parallel)
args = commandArgs(trailingOnly=TRUE)
#################################
VCF = args[1]
SAMPLES = args[2]
WINSIZE = args[3]
CPU = args[4]
INTERMEDIATE = 'intermediate/'
VCFTOOLS = '/app/vcftools-0.1.16/bin/vcftools'
#VCFTOOLS = 'vcftools'
#################################

################### Functions ###################
FUN_tajima <- function(LFMM_vcf, samples, INTERMEDIATE, VCFTOOLS, WINSIZE, CPU){
  # Don't work for some reason (Calculated not from docker)

  mclapply(samples$site %>% unique, function(site){
	COMMAND = paste(VCFTOOLS, '--vcf', LFMM_vcf,
                   '--TajimaD', WINSIZE, # in each position 
                   paste('--indv', samples$sample[samples$site == site], collapse = ' '),
                   '--out', paste0(INTERMEDIATE, 'pop', site)
                   )
							 	message(COMMAND)
	system(COMMAND)
    }, mc.cores = CPU)

  FILENAMES = paste0(INTERMEDIATE, 'pop', samples$site %>% unique, '.Tajima.D')
  tajima <-
    lapply(FILENAMES, function(f) {
      popname = gsub(paste0(INTERMEDIATE, 'pop'), '', f) %>% 
	      	gsub('.Tajima.D', '', .) # extract popname
      fread(f, header = T) %>%
        dplyr::mutate(site = as.factor(popname)) %>%
	dplyr::group_by(site) %>%
        dplyr::summarise(TajimaD = mean(TajimaD, na.rm = T))
      }) %>% do.call(rbind, .)

  return(tajima)
}

################## Main ##########################
samples <- fread(SAMPLES, 
		 colClasses = c("site" = "character", 'sample' = 'character'))
							message(samples %>% str)
							message("INFO: RUN NEUTRALITY TEST FOR WHOLE GENOME")
FUN_tajima(VCF, samples, INTERMEDIATE, VCFTOOLS, WINSIZE, CPU) %>%
	fwrite('tables/TajimaD_TotalByPop.tsv', sep = '\t')
