library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
###################
METADATA=args[1] %>% fread(colClasses = c("site" = "character", "sample" = "character")) # METADATA raw table
SAMPLES_VCF=args[2] %>% fread(header = F, colClasses = "character")
message(SAMPLES_VCF)
message(METADATA %>% str)
OUTPUT=args[3]
###################

dt <-
	METADATA %>%
		dplyr::filter(sample %in% !!SAMPLES_VCF$V1) %>%
		dplyr::arrange(match(sample, !!SAMPLES_VCF$V1)) 

# Save filtred data
dt %>%
	fwrite(OUTPUT, sep = '\t')

