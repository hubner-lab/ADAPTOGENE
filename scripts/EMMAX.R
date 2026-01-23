library(dplyr)
library(data.table)
library(magrittr)
library(purrr)
library(stringr)
library(qvalue)
########################
args = commandArgs(trailingOnly=TRUE)
VCF=args[1]
Kbest=args[2] %>% as.numeric
trait=fread(args[3], sep = '\t', header = T) # tsv file with set of traits
covariates=fread(args[4], sep = ' ', header = F)[,3:(2 + Kbest)] %>% setNames(paste0('PC', 1:Kbest))# MAKE SURE THAT SAMPLES IN THE SAME ORDER AS IN VCF
PREDICTORS_SELECTED=args[5] %>% str_split(',') %>% unlist
OUT=args[6]
SAMPLES = args[7]
########################
FUN_emmax <- function(VCF,  #TODO everytime make convertion of input file, should create separate function
                      trait, # Trait (1 col df)
                      covariates, # NULL if not use correction
		      OUT,
                      kinship = T, # use kinship correction
                      force = F){

  f = paste0(gsub('.vcf', '', basename(VCF)))
  tfam = paste0(f, '.tfam')
  kinship.f = paste0(f, '.aIBS.kinf')
  											message("INFO: Read VCF and extract SNPID")
  snpid = fread(cmd = paste('grep -v \"##\"', VCF, '| cut -f1-2')) %>% 
	  setNames(c('CHROM', 'POS')) %>% 
	  dplyr::mutate(SNPID = paste0(CHROM, ':', POS)) %>%
	  .$SNPID
 											message("INFO: Read VCF and extract SampleNames") 
  SampleName = fread(cmd = paste("head -1000", VCF, "| grep \'#CHROM\' | cut -f10- "), header = F) %>%
    as.character()


  if(!file.exists(kinship.f) | force){
                                                                                                                      message(paste0('EMMAX: Create ',f,' file'))

    # Prepare TPED/TFAM files
    system(paste0('plink --vcf ', VCF,
                  '  --allow-extra-chr --recode12 transpose --output-missing-genotype 0 --out ', f))
    # Edit family and individual names in .tfam
#    system(paste("cat", tfam, " | awk \'{ print $1 \"_\" $2, $1 \"_\" $2,$3,$4,$5,$6}\' > tmp.tfam ; mv tmp.tfam ", tfam))

    system(paste("cat", tfam, " | awk \'{split($1, a, \"_\"); split($2, b, \"_\"); if (a[1] == b[1]) {$1 = a[1]; $2 = a[1]} print}\' > tmp.trfam && mv tmp.tfam ", tfam, " && rm tmp.tfam"))	

    message(paste0('EMMAX: Calculate kinship matrix'))
    # Run Kinship matrix calculation
    system(paste('../scripts/emmax-kin-intel64 -v -s -d 10 -x', f))
  }

  # FILES produced
  traitname = colnames(trait)[1]
  PHEN = paste0(OUT,'EMMAX_phenotype_', '_', traitname, '.tsv')
  COVAR = paste0(OUT,'EMMAX_covariates_', '_', traitname, '.tsv')
  EMMAXOUT = paste0(OUT,'EMMAX_OUT_', traitname)
	
  # Prepare phenotype file
  emmax.phen <-
    trait %>%
    dplyr::mutate(FAMID = SampleName,
                  INDID = SampleName) %>%
    dplyr::select(FAMID, INDID, everything()) %T>%
    write.table(PHEN, sep = '\t', col.names = F, row.names = F, quote = F)

  if(!is.null(covariates) & !is.null(kinship)){
                                                                                      message('RUN emmax with PCA and kinship corrections')

    # Prepare CV file
    emmax.phen %>%
      dplyr::select(FAMID, INDID) %>%
      dplyr::mutate(smth = 1) %>% # intercept
      cbind(covariates) %>%
      write.table(COVAR, sep = '\t', col.names = F, row.names = F, quote = F)

    # Run EMMAX
    system(paste('../scripts/emmax-intel64 -v -d 10 -t', f, 
		 ' -p', PHEN,
		' -k', kinship.f, 
		'-c', COVAR,
                 '-o', EMMAXOUT)
           )
  }

  if(is.null(covariates) & !is.null(kinship)){

                                                                                              message('RUN emmax without PCA correction')
      # Run EMMAX
    system(paste('../scripts/emmax-intel64 -v -d 10 -t', f, 
		 ' -p', PHEN,
		 '-k', kinship.f,
                 '-o', EMMAXOUT)
           )

  }

  if(is.null(covariates) & is.null(kinship)){
                                                                                  message('RUN emmax without PCA and kinship corrections')

      # Run EMMAX
    system(paste('../scripts/emmax-intel64 -v -d 10 -t', f, 
		 ' -p', PHEN,
                 '-o', EMMAXOUT)
           )

  }

  EMMAXOUT.ps = paste0(EMMAXOUT, '.ps')

  # load df
  df <-
    fread(EMMAXOUT.ps) %>% 
    	setNames(c('SNPID', 'beta', 'SE.beta', traitname)) %>%
    	dplyr::mutate(SNPID = !!snpid,
			  chr = SNPID %>% str_extract("(.*):", group = 1),
			  pos = SNPID %>% str_extract(":(.*)", group = 1) %>% as.numeric ) %>%
    	dplyr::select(SNPID, chr, pos, everything() ) %>%
	dplyr::select(-beta, -SE.beta)

  return(df)
}

# Run each trait one by one
trait <-
	trait %>% dplyr::select(!!PREDICTORS_SELECTED)

# Add pheno traits if we have it
samples = fread(SAMPLES)

if(ncol(samples) > 4){
	trait <-
		cbind(trait, 
		      samples %>% dplyr::select(-site, -sample, -latitude, -longitude)
		      )
}

pval_dt <-
	lapply(1:ncol(trait), function(i) {
		FUN_emmax(VCF, trait[, ..i], covariates, OUT)
		   }) %>%
	reduce(function(x, y) {
	    left_join(x, y, by = c("SNPID", "chr", "pos"))
	  })
								message(pval_dt %>% str)
qval_dt <-
	lapply(pval_dt %>% dplyr::select(-SNPID, -chr, -pos), function(biovec){
		qvalue(biovec)$qvalues
	  }) %>%
	do.call(cbind, .) %>%
	cbind(pval_dt %>% dplyr::select(SNPID, chr, pos), 
	      .)

pval_dt %>% 
	fwrite(paste0("tables/EMMAX/EMMAX_pvalues_K", Kbest, ".tsv"), sep = '\t') #TODO fix to return Kbest as PC number of covariate used
qval_dt %>% 
	fwrite(paste0("tables/EMMAX/EMMAX_qvalues_K", Kbest, ".tsv"), sep = '\t')
#TODO also add choice of kinship matrix algorithm

