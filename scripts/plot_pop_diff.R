library(LEA)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(qs)

args = commandArgs(trailingOnly=TRUE)
#################################
SNMF_PROJECT = args[1]
K = args[2] %>% as.numeric
PLOIDY = args[3] %>% as.numeric
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

# Population differentiation tests
p <- snmf.pvalues(snmf_res,
                  entropy = TRUE,
                  ploidy = PLOIDY,
                  K = K)
pvalues.df <- data.frame(pvalues = p$pvalues)

# Histogram plot
gHist <-
  ggplot(pvalues.df, aes(x = pvalues)) +
    geom_histogram(aes(y = ..density..), color = 'black', fill = 'lightblue') +
    geom_density(alpha = 0.2, fill = '#FF6666') +
    plot_theme

# Manhattan-style plot
gPval <-
  ggplot(pvalues.df, aes(y = -log10(pvalues), x = 1:nrow(pvalues.df))) +
    geom_point(color = 'blue') +
    labs(x = 'Index') +
    plot_theme

# Combine plots
gGrid <- ggarrange(gHist, gPval, nrow = 2)

# Save
ggsave(paste0(PLOT_DIR, 'pop_diff_K', K, '.png'), gGrid, width = 2 * 6.4, height = 9.6)
ggsave(paste0(PLOT_DIR, 'pop_diff_K', K, '.svg'), gGrid, width = 2 * 6.4, height = 9.6,
       device = svglite::svglite, bg = "transparent")
qsave(gGrid, paste0(INTER_DIR, 'pop_diff_K', K, '.qs'))
