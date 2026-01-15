library(topr) 
packageVersion('topr') # It's important to use > 2.0.0 version, previous versions doesn't work with non-human genomes
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
library(magrittr)
#library(cowplot)
library(GenomicRanges)
library(qvalue)
library(egg)
library(gridExtra)
library(grid)
library(svglite)
args = commandArgs(trailingOnly=TRUE)

###################################### Arguments
ASSOC_TABLE = args[1] # pvalues table for all traits
GFF_TOPR = args[2] # gff transformed into topr annotation format
REGION=args[3] # chr:start-end
ADJUST=args[4] #method_num for threshold line
TRAITS=args[5] %>% str_split(',') %>% unlist # trait names separated by comma
GENES_TO_HIGHLIGHT = args[6] %>% str_split('\\|') %>% unlist # list of genes (sep = '|') to highlight (could be 'all', nameless)
##################################### Functions
# Generalized extraction and concatenation of all digits


extract_and_concatenate <- function(chrom_names) {
  # Extract all digits and concatenate them
  sapply(chrom_names, function(x) {
    paste0(unlist(regmatches(x, gregexpr("[0-9]+", x))), collapse = "") %>% as.numeric %>% unname
  })
}

is_valid_region <- function(region_string) {
  # Extract start and end using regex
  start <- str_extract(region_string, "(?<=:)[0-9]+")
  end   <- str_extract(region_string, "(?<=-)[0-9]+")

  # Check that both start and end are not NA and are numeric
  if (is.na(start) || is.na(end)) return(FALSE)

  start <- as.numeric(start)
  end   <- as.numeric(end)

  # Final check: start and end are numbers and start <= end
  return(!is.na(start) && !is.na(end) && start <= end)
}
##################################### MAIN

# Parse arguments

## Check regi
if(!is_valid_region(REGION)){ stop('The region positions are not numbers, check it!') }

# Pvalues
assoc <- fread(ASSOC_TABLE)
# Threshold
method = ADJUST %>% str_extract('(.*)_', group = 1)
signif = ADJUST %>% str_extract('_(.*)', group = 1) %>% as.numeric

if(method == 'bonf'){ pval = signif / nrow(assoc) }
if(method == 'qval'){
	qvals = qvalue(assoc$pvalue)$qvalues
	pval = max(assoc$pvalue[qvals <= signif], na.rm = TRUE)
} # here signif will be desirable threshold for qvalue corrected pvalue, we need to determine the corresponding pvalue number
if(method == 'custom'){ pval = signif }
if(method == 'top'){ pval = sort(assoc$pvalue)[min(signif, nrow(assoc))] } # here signif will be top 100 or 1000 and etc
# required format: SNP     CHROM     POS     P
# current: SNPID   chr     pos     bio_1 

assoc_lst <-
	lapply(TRAITS, function(trait){
		assoc %>% 
			dplyr::select(SNPID, chr, pos, !!trait) %>% 
			setNames(c('SNP', 'CHROM', 'POS', 'P')) %>%
			dplyr::mutate(CHROM = extract_and_concatenate(CHROM) %>% unname)
  	}) %>% setNames(TRAITS)
	
#TODO one of the problem is chr suffixes so it should be solved before running the pipeline or changed here somehow (hard to generalize)
# One of the solution could provide option for create new labels for chromosomes 


# Parse GFF
gff_table <-
	fread(GFF_TOPR) %>%
		dplyr::mutate(chrom = extract_and_concatenate(chrom) %>% unname )
if(!GENES_TO_HIGHLIGHT %in% c('nameless', 'all')) { # keep only names in list
	gff_table <-
		gff_table %>%	
			dplyr::filter(gene_symbol %in% GENES_TO_HIGHLIGHT) 
}
if(GENES_TO_HIGHLIGHT == 'nameless'){ # erase names but plot gene models
	gff_table <-
		gff_table %>%
			dplyr::mutate(gene_symbol = '')
}								

if(nrow(gff_table) == 0) { message("No genes found in region.") }

chr = REGION %>% str_extract('(.*):', group = 1)
startend = REGION %>% str_extract(':(.*)', group = 1)

chrnum = chr %>% extract_and_concatenate() %>% unname
REGION_chrnum = paste0(chrnum, ':', startend)

#mclapply(1:length(regions_expanded_lst), function(i){
# Sanitize REGION if needed (e.g., replace ":" and "-" with "_")
REGION_SAFE <- gsub("[:\\-]", "_", REGION)
# Construct filename
filename_base <- paste0('plots/regionPlot_',
                   TRAITS %>% paste(collapse = '_'), 
		   '_',
                   REGION_SAFE, '_',
		   ADJUST
		   )

#TEMP change color
#colors = c('#4197d8', '#f8c120', 'black',
#  setdiff(topr::get_topr_colors(), c("darkblue", "#E69F00", "#00AFBB"))
#)
	 									message(gff_table %>% str)
										message(assoc_lst %>% str)
										message(pval)
										message(REGION_chrnum)
										message(REGION_SAFE)
										message(filename_base)
										#message(colors)

png(paste0(filename_base, '.png'), width = 12, height = 5, units = 'in', res = 300) #TODO export width,heigth to config
regionplot(assoc_lst, # https://rdrr.io/cran/topr/man/regionplot.html
	   build=gff_table, #TODO doesnt' work with new version
	   region = REGION_chrnum, # chr:start-end 
	   sign_thresh = pval, #TODO
	   max.overlaps = 999, 
	   show_gene_names = F,
	   show_genes = F,

	   unit_overview = 1,
	   unit_main = 2.5,
	   unit_gene = 1,
	   unit_ratios = 1:2.5:1,

	   #show_gene_legend = F,
	   #gene_label_size = 2,
	   # nudge_y = 1,
	   # nudge_x = 20,

	   axis_text_size = 16,
	   axis_title_size = 18,
	   title_text_size = 18,
	   legend_position = 'right', #TODO keep legend!
	   scale = 1,
	   segment.size = 3,
	   legend_labels = TRAITS,
	   sign_thresh_label_size = 0,
	   extract_plots = F
	   #,color = colors
)
dev.off()

plot_list <-
	regionplot(assoc_lst, # https://rdrr.io/cran/topr/man/regionplot.html
	   build=gff_table,  #TODO doesn't work with new version
	   region = REGION_chrnum, # chr:start-end 
	   sign_thresh = pval, #TODO
	   max.overlaps = 999, 
	   show_gene_names = F,
	   show_genes = F,
	
	   unit_overview = 1,
	   unit_main = 2.5,
	   unit_gene = 1, #TODO figured out the problem with genes
	   unit_ratios = 1:2.5:1,

	   #show_gene_legend = F,
	   #gene_label_size = 2,
	   # nudge_y = 1,
	   # nudge_x = 20,

	   axis_text_size = 16,
	   axis_title_size = 18,
	   title_text_size = 18,
	   legend_position = 'right', #TODO keep legend!
	   scale = 1,
	   segment.size = 3,
	   legend_labels = TRAITS,
	   sign_thresh_label_size = 0,
	   extract_plots = T
	   #,color = colors
	)
# Save each plot individually
ggsave(paste0(filename_base, '_', 'main.svg'), plot = plot_list$main_plot, device = svglite, width = 10, height = 6)
ggsave(paste0(filename_base, '_', 'overview.svg'), plot = plot_list$overview_plot, device = svglite, width = 10, height = 2)
ggsave(paste0(filename_base, '_', 'genes.svg'), plot = plot_list$gene_plot, device = svglite, width = 12, height = 2)
