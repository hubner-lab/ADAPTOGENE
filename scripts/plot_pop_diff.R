library(LEA)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(qs)
library(scattermore)

args = commandArgs(trailingOnly=TRUE)
#################################
SNMF_PROJECT = args[1]
K = args[2] %>% as.numeric
PLOIDY = args[3] %>% as.numeric
PLOT_DIR = args[4]
INTER_DIR = args[5]
SCATTERMORE_THRESHOLD = args[6] %>% as.numeric
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
pvalues.df$idx <- 1:nrow(pvalues.df)
pvalues.df$log10p <- -log10(pvalues.df$pvalues)

# Decide rendering method based on SNP count
n_snps <- nrow(pvalues.df)
use_scattermore <- n_snps > SCATTERMORE_THRESHOLD
if (use_scattermore) {
    message(paste0('INFO: Using scattermore for fast rendering (', n_snps, ' > ', SCATTERMORE_THRESHOLD, ')'))
} else {
    message(paste0('INFO: Using standard geom_point (', n_snps, ' <= ', SCATTERMORE_THRESHOLD, ')'))
}

# Histogram plot
gHist <-
  ggplot(pvalues.df, aes(x = pvalues)) +
    geom_histogram(aes(y = after_stat(density)), color = 'black', fill = 'lightblue') +
    geom_density(alpha = 0.2, fill = '#FF6666') +
    plot_theme

# Manhattan-style plot
gPval <- ggplot(pvalues.df, aes(x = idx, y = log10p))

if (use_scattermore) {
    # pixels should match output: 12.8in x 4.8in (half height) @ 300dpi = ~3840x1440
    gPval <- gPval +
        geom_scattermore(color = 'blue', pointsize = 12, pixels = c(3840, 1440),
                         interpolate = FALSE)
} else {
    gPval <- gPval +
        geom_point(color = 'blue')
}

gPval <- gPval +
    labs(x = 'Index', y = expression(-log[10](p-value))) +
    plot_theme

# Combine plots
gGrid <- ggarrange(gHist, gPval, nrow = 2)

# Save
ggsave(paste0(PLOT_DIR, 'pop_diff_K', K, '.png'), gGrid, width = 2 * 6.4, height = 9.6)
ggsave(paste0(PLOT_DIR, 'pop_diff_K', K, '.svg'), gGrid, width = 2 * 6.4, height = 9.6,
       device = svglite::svglite, bg = "transparent")
qsave(gGrid, paste0(INTER_DIR, 'pop_diff_K', K, '.qs'))
