library(ggplot2)
library(data.table)
library(dplyr)
library(qs)
library(ggcorrplot)

args = commandArgs(trailingOnly=TRUE)
#################################
CLIMATE = args[1]    # climate factors from sites
SAMPLES = args[2]    # samples with traits if exists
PLOT_DIR = args[3]
INTER_DIR = args[4]
#################################

# Plot theme
plot_theme <-
  theme_classic(base_size = 8, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = 'bold', size = 10, angle = -90),
        axis.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 12))

#################################

samples = fread(SAMPLES)
climate = fread(CLIMATE)

# Combine climate with phenotypic traits if available
if (ncol(samples) > 4) {
  traits = cbind(climate,
                 samples %>% dplyr::select(-site, -sample, -latitude, -longitude))
} else {
  traits = climate
}

# Remove constant columns
std_devs <- apply(traits, 2, sd, na.rm = TRUE)
non_constant_traits = names(std_devs)[std_devs != 0]
traits <- traits %>% dplyr::select(all_of(non_constant_traits))

# Compute correlation matrix
cor_matrix <- cor(traits, use = 'pairwise.complete.obs')

# Build the heatmap using ggcorrplot
gHM <- ggcorrplot(cor_matrix,
                  type = "lower",
                  lab = TRUE,
                  digits = 1,
                  lab_size = 4) +
  ggtitle('Correlogram of traits') +
  plot_theme

# Determine size based on number of traits
scale_factor <- 0.7
width <- ncol(traits) * scale_factor
height <- ncol(traits) * scale_factor

# Save
ggsave(paste0(PLOT_DIR, 'CorrelationHeatmap.png'), gHM, width = width, height = height, units = "in")
ggsave(paste0(PLOT_DIR, 'CorrelationHeatmap.svg'), gHM, width = width, height = height, units = "in",
       device = svglite::svglite, bg = 'transparent')
qsave(gHM, paste0(INTER_DIR, 'CorrelationHeatmap.qs'))
