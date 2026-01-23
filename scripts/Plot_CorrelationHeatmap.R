library(ggplot2)
library(data.table)
library(dplyr)
library(qs)
library(stringr)
library(reshape2)  # For melting data
library(scales)    # For better color scales
library(ggcorrplot)
args = commandArgs(trailingOnly=TRUE)
#################################
CLIMATE = args[1] # climate factors from sites
SAMPLES = args[2] # samples with traits if exists
# Plots settings
plot_theme <-
	        theme_classic(base_size = 8, base_family = 'Helvetica') +
		theme(plot.title = element_text(size = 14, face = 'bold'),
		      axis.text.y = element_text(face = "bold", size = 10),
		      axis.text.x = element_text(face = 'bold', size = 10, angle = -90),
		      axis.title = element_blank(),
		      legend.text = element_text(size = 12),
                      legend.title = element_text(face = "bold", size = 12)
									                    )
#################################

samples = fread(SAMPLES)
climate = fread(CLIMATE)

if(ncol(samples) > 4){ traits = cbind(climate, 
				     samples %>% dplyr::select(-site, -sample, -latitude, -longitude)
				     )
} else { traits = climate }

#TEMP
#traits = traits %>% dplyr::select(-bio_14, -bio_10, -bio_11, -bio_5, -bio_6, -bio_8, -bio_9, -bio_13, -bio_16, -bio_17, -bio_18, -bio_19)

# Remove constant columns
std_devs <- apply(traits, 2, sd, na.rm = TRUE)
							
non_constant_traits = names(std_devs)[std_devs != 0]
traits <- traits %>% dplyr::select(!!non_constant_traits)
							

# Compute correlation matrix with removing all rows with NAs
cor_matrix <- cor(traits, use = 'pairwise.complete.obs')
							
							

#cor_matrix %>%
#	fwrite('intermediate/CorMatrix_TEST.tsv', sep = '\t')
# Check for problematic values
#print(any(is.na(cor_matrix)))    # Check for NA
#print(any(is.nan(cor_matrix)))   # Check for NaN
#print(any(is.infinite(cor_matrix)))  # Check for Inf		

#cor_matrix[is.na(cor_matrix)] <- 0   # Replace NA
#cor_matrix[is.nan(cor_matrix)] <- 0 # Replace NaN
#cor_matrix[is.infinite(cor_matrix)] <- 0 # Replace Inf

# Build the heatmap using ggplot2
gHM <-
	ggcorrplot(cor_matrix,
		   #hc.order = TRUE,
		   type = "lower",
		   lab = TRUE,
		   digits = 1,
		   lab_size = 4
		   ) +
	ggtitle('Correlogram of traits') +
		plot_theme

# Determine size based on number of traits
scale_factor <- 0.7  # Adjust this value as needed
width <- ncol(traits) * scale_factor
height <- ncol(traits) * scale_factor

ggsave('plots/CorrelationHeatmap.png', gHM, width = width, height = height, units = "in")
qsave(gHM, 'intermediate/CorrelationHeatmap.qs')
