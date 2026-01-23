library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(parallel)
args = commandArgs(trailingOnly=TRUE)
################
GFF = args[1]
sigSNPs_list = args[2] %>% str_split(',') %>% unlist # list of sigSNPs files with different ADJUST
GFF_FEATURE = args[3]
DISTANCE = args[4] %>% as.numeric; DISTANCE_str = args[4]
PROMOTER_LENGTH = args[5] %>% as.numeric
allSNPs = args[6] #<SNPID chr pos trait1 trait2>
CPU = args[7] %>% as.numeric

################
message(sigSNPs_list %>% str)
######################## Function

# Count exon\promoter SNPs

# Extract gene ID from Parent or fallback to ID field (cleaned)
extract_gene_id <- function(attr) {
	id <- str_extract(attr, "(?<=Parent=)[^;]+")                # 1. Extract Parent if exists
	id <- ifelse(is.na(id), str_extract(attr, "(?<=ID=)[^;]+"), id)  # 2. If no Parent, take ID
	id <- str_remove(id, "\\.[0-9]+(_exon_[0-9]+)?$")            # 3. Remove .1 and/or _exon_#
	return(id)
}

FUN_find_overlapped_genes <- function(GFF, # path to gff file
				      gff_feature = 'mRNA', # WITHOUT HEADER and comment chars!
				      sigSNPs,  # or FDR or another adjusted pval => Should contain chr and pos col
				      DISTANCE = 1e6, # DISTANCE / 2 from start and from end of the genes, should be based on LD calculations
				      CPU = 1
				      ){

	snps_sig = fread(sigSNPs)

	traits = snps_sig$trait %>% unique


	genes_gff <- 
		fread(cmd = paste("grep -v \'#\'", 
				  GFF), 
		      header = F) %>%
		dplyr::filter(V3 == !!gff_feature) %>%
		dplyr::select(V1, V4, V5, V9) %>%
		setNames(c('chr', 'start', 'end', 'description')) %>%
		dplyr::mutate(id = description %>% extract_gene_id)

	# extract description fields
	# --- extract keys (field names) robustly from first description string ---
	fields <- genes_gff$description[1] %>%
	  str_extract_all("(?<=^|;)[^=;]+(?==)", simplify = TRUE) %>%
	  as.character() %>%
	  unique() %>%
	  .[. != ""]   # drop any empty

	# split description by ';' into those keys (values will still be like "ID=xxx")
	genes_gff <- genes_gff %>%
	  tidyr::separate(col = "description", into = fields, sep = ";", fill = "right", extra = "drop") %>%
	  dplyr::mutate(across(all_of(fields), ~ na_if(.x, "NA"))) %>%
	  dplyr::mutate(start = as.integer(start), end = as.integer(end))

	#message(genes_gff %>% str)	

	genes_gr <- 
		genes_gff %>%
		GRanges()
	
	genes_extended_gr <- 
		GRanges(
			seqnames = seqnames(genes_gr),
			ranges = IRanges(
					 start = pmax(1, start(genes_gr) - DISTANCE / 2),
					 end = end(genes_gr) + DISTANCE / 2)
		)
	
	###############
	
	# In GRanges
	snps_sig_gr <-
		snps_sig %>%
		dplyr::mutate(start = pos,
			      end = pos) %>%
		dplyr::select(chr, start, end, trait, method) %>%
		dplyr::mutate(SNPID = paste0(chr, ':', start)) %>%
		GRanges()
	
	#message(snps_sig_gr %>% str)
	#message(genes_extended_gr %>% str)
	#message(genes_gff %>% str)	
	#message(snps_sig_gr$trait)
	#message(snps_sig %>% str)
	# Find overlaps
	genes_overlapped <-
		mclapply(1:length(snps_sig_gr), function(i){
	
				 overlaps = findOverlaps(snps_sig_gr[i], genes_extended_gr)

				 genes_gff[subjectHits(overlaps),] %>%
					 dplyr::mutate(SNPID = snps_sig_gr$SNPID[i],
						       trait = snps_sig_gr$trait[i],
					 	       method  = snps_sig_gr$method[i])
	
			      }, mc.cores = CPU) %>% 
			do.call(rbind, .) %>%
			dplyr::select(trait, SNPID, method, everything())

  return(genes_overlapped)
}


# Function for counting SNPs in EXONs of the genes and in promoter region
# Return data.table ID,exonSNP_N, exonSNP_IDs, promoterSNP_N, promoterSNP_IDs
count_exon_promoter_snps <- function(genesID, # vec of genesID
				     snps, # chr,pos
				     GFF,
				     promoter_len = PROMOTER_LENGTH, # define length of promoter region
				     gff_feature	
				     ) {

	gff_like_dt <- 
		fread(cmd = paste("grep -v \'#\' ", GFF)) %>% 
		dplyr::select(V1, V3, V4, V5, V9) %>%
		setNames(c('chr', 'feature', 'start', 'end', 'description')) %>%
		dplyr::filter(feature %in% c(!!gff_feature, 'exon')) %>%
		dplyr::mutate(id = description %>% extract_gene_id) %>%
		dplyr::select(-description)

	exons_gr <-
		gff_like_dt %>%
		dplyr::filter(id %in% !!genesID & feature == 'exon') %>%
		dplyr::select(-feature) %>%
		GRanges()

	promoters_gr <-
		gff_like_dt %>%
		dplyr::filter(id %in% !!genesID & feature == !!gff_feature) %>%
		dplyr::mutate(end = start,
			      start = start - promoter_len
			      ) %>%
		dplyr::select(-feature) %>%
		GRanges()

	snps_gr <-
		snps %>%
		dplyr::mutate(start = pos,
			      end = pos) %>%
		dplyr::select(chr, start, end) %>%
		GRanges()

	res_exon <-
		findOverlapPairs(snps_gr, exons_gr) %>% 
		as.data.table %>%
		dplyr::select(first.seqnames, first.start, second.id) %>%
		setNames(c('chr', 'pos', 'id')) %>%
		dplyr::mutate(SNPID = paste0(chr, ':', pos)) %>%
		group_by(id) %>%
		dplyr::summarise(exonSNP_N = n(),
				 exonSNP_IDs = SNPID %>% paste(collapse = ','))
		res_promoter <-
			findOverlapPairs(snps_gr, promoters_gr) %>%
			as.data.table %>%
			dplyr::select(first.seqnames, first.start, second.id) %>%
			setNames(c('chr', 'pos', 'id')) %>%
			dplyr::mutate(SNPID = paste0(chr, ':', pos)) %>%
			group_by(id) %>%
			dplyr::summarise(promoterSNP_N = n(),
					 promoterSNP_IDs = SNPID %>% paste(collapse = ','))

			return(full_join(res_promoter, res_exon, by = 'id'))
}
#################### MAIN

message('start')
message(allSNPs)
message(args)

allsnps <- 
	fread(allSNPs) %>%
	dplyr::select(V1,V2) %>%
	setNames(c('chr', 'pos'))

message(allsnps %>% str)
message('end')

lapply(sigSNPs_list, function(sigSNPs){
	       message('FUN_overlap')
	       genes_overlapped_dt <-
		       FUN_find_overlapped_genes(GFF = GFF,
						 gff_feature = GFF_FEATURE,
						 sigSNPs = sigSNPs,
						 DISTANCE = DISTANCE,
						 CPU = CPU)
	       message(genes_overlapped_dt %>% str)
		message('FUN_exon')
	       exon_promoter_snp_counts <- 
		       count_exon_promoter_snps(genes_overlapped_dt$id, # use ID extracted from ID or Parent field of description (update extract ID function for more flexability)
						allsnps,
						GFF,
						PROMOTER_LENGTH,
		       				GFF_FEATURE) 

	       message(exon_promoter_snp_counts %>% str)
		message('JOINING')
	       # Join
	       res_dt <- 
		       genes_overlapped_dt %>%
		       left_join(exon_promoter_snp_counts, by = 'id') %>% # Parent column was created based on ID or Parent of gff records 
		       dplyr::arrange(desc(exonSNP_N), desc(promoterSNP_N)) %>%
		       as.data.table
		
	       # Collapse by gene (Parent)
	       res_dt_collapsed <- res_dt[, lapply(.SD, function(x) {
							   if (length(unique(x)) == 1) {
								   as.character(x[1])
							   } else {
								   paste(x, collapse = ",")
							   }
						}), by = id] %>%
	       dplyr::mutate(trait_N = trait %>% str_split(',') %>% sapply(length),
	       		     SNPID_N = SNPID %>% str_split(',') %>% sapply(length) ) %>%
	       dplyr::arrange(chr, start)



	       # Save
	       message(genes_overlapped_dt %>% str)
	       message(res_dt %>% str)
	       message(sigSNPs %>% str)
	       message(res_dt_collapsed %>% str)

	       res_dt %>%
		       fwrite(paste0(sigSNPs %>% gsub('.tsv$', '', .) %>% gsub('_sigSNPs', '', .), '_genesAround', DISTANCE_str, '.tsv'),
			      col.names = T, sep = '\t') # should be based on input file
	       
	       res_dt_collapsed %>%
		       fwrite(paste0(sigSNPs %>% gsub('.tsv$', '', .) %>% gsub('_sigSNPs', '', .), '_genesAround', DISTANCE_str, '_collapsed.tsv'),
			      col.names = T, sep = '\t') # should be based on input file
					 })

#_sigSNPs_bonf_0.05_genesAround.tsv
