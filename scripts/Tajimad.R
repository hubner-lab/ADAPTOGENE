library(dplyr)
library(data.table)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
#################################
VCF = args[1]
SAMPLES = args[2]
WINSIZE = args[3]
CPU = args[4] %>% as.numeric
OUTPUT = args[5]       # output TSV path
INTER_DIR = args[6]    # directory for intermediate vcftools files
#################################

VCFTOOLS = '/app/vcftools-0.1.16/bin/vcftools'

################### Functions ###################
FUN_tajima <- function(vcf_file, samples, inter_dir, vcftools_path, winsize, cpu) {

  mclapply(samples$site %>% unique, function(site) {
    COMMAND = paste(vcftools_path, '--vcf', vcf_file,
                    '--TajimaD', winsize,
                    paste('--indv', samples$sample[samples$site == site], collapse = ' '),
                    '--out', paste0(inter_dir, 'pop', site))
    message(COMMAND)
    system(COMMAND)
  }, mc.cores = cpu)

  FILENAMES = paste0(inter_dir, 'pop', samples$site %>% unique, '.Tajima.D')
  
  tajima <- lapply(FILENAMES, function(f) {
    popname = gsub(paste0(inter_dir, 'pop'), '', f) %>%
      gsub('.Tajima.D', '', .)
    fread(f, header = TRUE) %>%
      dplyr::mutate(site = as.factor(popname)) %>%
      dplyr::group_by(site) %>%
      dplyr::summarise(TajimaD = mean(TajimaD, na.rm = TRUE))
  }) %>% do.call(rbind, .)

  return(tajima)
}

################## Main ##########################
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character", 'sample' = 'character'))
message(samples %>% str)
message("INFO: RUN NEUTRALITY TEST FOR WHOLE GENOME")

FUN_tajima(VCF, samples, INTER_DIR, VCFTOOLS, WINSIZE, CPU) %>%
  fwrite(OUTPUT, sep = '\t')
