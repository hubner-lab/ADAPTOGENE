library(LEA)
library(ggplot2)
library(dplyr)
library(data.table)
library(qs)

args = commandArgs(trailingOnly=TRUE)
#################################
LFMM = args[1]
SAMPLES = args[2]
PLOT_DIR = args[3]
INTER_DIR = args[4]
#################################

# Load samples info
samples.df <- fread(SAMPLES,
                    colClasses = c("site" = "character", 'sample' = 'character')) %>%
              dplyr::select(sample, site)

# Run PCA
pc <- LEA::pca(LFMM, scale = TRUE)

# Preprocess for plotting
			message(pc$projections %>% str)
PCA <- as.data.frame(pc$projections) %>%
       setNames(paste0('PC', 1:ncol(.)))
isrPCA <- cbind(samples.df, PCA)

var.explained <- c((pc$eigenvalues[1,] / sum(pc$eigenvalues)) * 100,
                   (pc$eigenvalues[2,] / sum(pc$eigenvalues)) * 100)

# PCA plot
gPCA <- ggplot(isrPCA, aes(x = PC1, y = PC2,
                           label = as.factor(site),
                           color = as.factor(site))) +
        geom_label() +
        xlab("PC1") +
        ylab("PC2") +
        guides(color = guide_legend(ncol=2)) +
        theme_bw()

ggsave(paste0(PLOT_DIR, "pca.png"), gPCA)
ggsave(paste0(PLOT_DIR, "pca.svg"), gPCA, device = svglite::svglite, bg = 'transparent')
qsave(gPCA, paste0(INTER_DIR, "pca.qs"))

# Tracy-Widom test
tw <- tracy.widom(pc)

gTW <- ggplot(data = tw, aes(x = N, y = percentage)) +
       geom_point() +
       xlim(1, 25) +
       theme_classic() +
       ggtitle("Tracy Widom") +
       theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(PLOT_DIR, "tracy_widom.png"), gTW)
ggsave(paste0(PLOT_DIR, "tracy_widom.svg"), gTW, device = svglite::svglite, bg = 'transparent')
qsave(gTW, paste0(INTER_DIR, "tracy_widom.qs"))
