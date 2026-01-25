library(ggplot2)
library(dplyr)
library(data.table)
library(qs)

args = commandArgs(trailingOnly=TRUE)
#################################
CLUSTERS = args[1]
K = args[2] %>% as.numeric
PLOT_DIR = args[3]
INTER_DIR = args[4]
#################################

# Colors
my.colors = c("#f3c300", "#a1caf1", "#be0032", "#8db600", "#e68fac", "#0067a5", "#f38400", "#222222",
              "#b3446c", "#2b3d26", "#dcd300", "#875692", "#848482", "#c2b280", "#882d17", "#654522",
              "#e25822", "#008856", "#f99379", "#604e97", "#f6a600")

# Load clusters
clusters <- fread(CLUSTERS)

# Make wide table
IsrPopKs <- clusters %>%
            dplyr::select(-site) %>%
            melt.data.table(id.vars = 'sample')

# Plot
gStructure <-
  ggplot(data = IsrPopKs, aes(y = value, x = sample, fill = variable)) +
    geom_bar(show.legend = TRUE, stat = "identity", position = "fill") +
    ylab("Proportion of assignment") +
    xlab("Accessions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values = my.colors[1:K]) +
    theme(strip.text = element_text(size = 44, angle = 0),
          axis.text.x = element_blank()) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20, face = 'bold')) +
    theme(strip.text.x = element_text(size = 16)) +
    theme(strip.text.y = element_text(size = 16)) +
    guides(fill = guide_legend(title = "Clusters"))

# Save
ggsave(paste0(PLOT_DIR, 'structure_K', K, '.png'), gStructure, width = 3 * 6.4, height = 4.8)
ggsave(paste0(PLOT_DIR, 'structure_K', K, '.svg'), gStructure, width = 3 * 6.4, height = 4.8,
       device = svglite::svglite, bg = "transparent")
qsave(gStructure, paste0(INTER_DIR, 'structure_K', K, '.qs'))
