library(LEA)
library(ggplot2)
library(dplyr)
library(data.table)
library(scatterpie)
library(qs)

args = commandArgs(trailingOnly=TRUE)
#################################
LFMM = args[1]
CLUSTERS = args[2]
EIGENVECTORS = args[3]
EIGENVALUES = args[4]
K = args[5] %>% as.numeric
PLOT_DIR = args[6]
INTER_DIR = args[7]
#################################

# Colors
my.colors = c("#f3c300", "#a1caf1", "#be0032", "#8db600", "#e68fac", "#0067a5", "#f38400", "#222222",
              "#b3446c", "#2b3d26", "#dcd300", "#875692", "#848482", "#c2b280", "#882d17", "#654522",
              "#e25822", "#008856", "#f99379", "#604e97", "#f6a600")

# Plot theme
plot_theme <-
  theme_classic(base_size = 12, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 17, face = 'bold'),
        axis.text = element_text(face = "bold", size = 18),
        axis.title = element_text(face = "bold", size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(face = "bold", size = 22))

# Load clusters
clusters <- fread(CLUSTERS)

		    message(clusters %>% str)
# Preprocess - use paste0 for naming, not gsub
# Load PCA
eigenvectors <- fread(EIGENVECTORS, header = F) %>%
           setNames(paste0('PC', 1:ncol(.)))
		    message(eigenvectors %>% str)
eigenvalues <- fread(EIGENVALUES, header = F)$V1
    message(eigenvalues %>% str)

# Validate dimensions
if (nrow(clusters) != nrow(eigenvectors)) {
    stop(paste0("ERROR: Row count mismatch! Clusters has ", nrow(clusters),
                " rows but PCA projections has ", nrow(eigenvectors), " rows."))
}

df <- cbind(clusters, eigenvectors[, 1:4])

var.explained <- c((eigenvalues[1] / sum(eigenvalues)) * 100,
                   (eigenvalues[2] / sum(eigenvalues)) * 100)

# Plot PCA with SNMF clusters
gPCA <- ggplot(df, aes(PC1, PC2, label = as.factor(site))) +
  geom_scatterpie(data = df,
                  aes(x = PC1, y = PC2),
                  cols = paste0('C', 1:K),
                  color = 'black', alpha = 0.8) +
  scale_fill_manual(values = my.colors[1:K]) +
  xlab(paste0('PC1 (', round(var.explained[1], 1), '%)')) +
  ylab(paste0('PC2 (', round(var.explained[2], 1), '%)')) +
  plot_theme

# Save
ggsave(paste0(PLOT_DIR, 'pca_structure_K', K, '.png'), gPCA)
ggsave(paste0(PLOT_DIR, 'pca_structure_K', K, '.svg'), gPCA,
       device = svglite::svglite, bg = "transparent")
qsave(gPCA, paste0(INTER_DIR, 'pca_structure_K', K, '.qs'))
