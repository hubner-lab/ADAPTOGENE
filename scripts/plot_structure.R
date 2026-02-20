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

# Colors — colorblind-safe (Wong 2011 + Paul Tol), visible on yellow/orange climate rasters
my.colors = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9",
              "#332288", "#882255", "#44AA99", "#AA4499", "#999933",
              "#E69F00", "#661100", "#88CCEE", "#CC6677", "#117733",
              "#DDDDDD", "#855C75", "#D9AF6B", "#736F4C", "#526A83", "#625377")

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
