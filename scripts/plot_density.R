library(ggplot2)
library(data.table)
library(dplyr)
library(qs)

args = commandArgs(trailingOnly=TRUE)
#################################
CLIMATE = args[1]    # climate site values TSV
BIO = args[2]        # bio variable name (e.g., "bio_1")
PLOT_DIR = args[3]
INTER_DIR = args[4]
#################################

# Plot theme
plot_theme <-
  theme_classic(base_size = 12, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(face = 'bold', size = 12),
        axis.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = "bold", size = 15))

# Load climate data
climate <- fread(CLIMATE)

# Check if bio variable exists
if (!BIO %in% colnames(climate)) {
  stop(paste0("ERROR: Variable '", BIO, "' not found in climate data. ",
              "Available: ", paste(colnames(climate), collapse = ", ")))
}

# Create density plot
gDens <- climate %>%
  ggplot(aes(x = .data[[BIO]])) +
  geom_density(fill = 'blue') +
  labs(title = BIO, x = "Values", y = "Density") +
  plot_theme

# Save
ggsave(paste0(PLOT_DIR, 'DensityPlot_', BIO, '.png'), gDens)
ggsave(paste0(PLOT_DIR, 'DensityPlot_', BIO, '.svg'), gDens,
       device = svglite::svglite, bg = 'transparent')
qsave(gDens, paste0(INTER_DIR, 'DensityPlot_', BIO, '.qs'))
