library(LEA)
library(ggplot2)
library(dplyr)
library(tibble)
library(qs)

args = commandArgs(trailingOnly=TRUE)
#################################
SNMF_PROJECT = args[1]
K_START = args[2] %>% as.numeric
K_END = args[3] %>% as.numeric
PLOT_DIR = args[4]
INTER_DIR = args[5]
#################################

# Plot theme
plot_theme <-
  theme_classic(base_size = 12, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 17, face = 'bold'),
        axis.text = element_text(face = "bold", size = 18),
        axis.title = element_text(face = "bold", size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(face = "bold", size = 22))

# Load SNMF project
snmf_res <- load.snmfProject(SNMF_PROJECT)

# Calculate cross-entropy for each K
message(paste0("INFO: Run Cross Entropy test on K", K_START, '-', K_END))

gCrossEntropy <-
  sapply(K_START:K_END, function(K){ min(cross.entropy(snmf_res, K)) }) %>%
  setNames(K_START:K_END) %>%
  enframe %>%
  dplyr::arrange(as.numeric(name)) %>%
  ggplot(aes(x = paste0('K', name), y = value)) +
    geom_point() +
    plot_theme +
    labs(x = 'Number of ancestral populations',
         y = 'Cross-entropy')

# Save outputs
ggsave(paste0(PLOT_DIR, 'cross_entropy_K', K_START, '-', K_END, '.png'), gCrossEntropy)
ggsave(paste0(PLOT_DIR, 'cross_entropy_K', K_START, '-', K_END, '.svg'), gCrossEntropy,
       device = svglite::svglite, bg = 'transparent')
qsave(gCrossEntropy, paste0(INTER_DIR, 'cross_entropy_K', K_START, '-', K_END, '.qs'))
