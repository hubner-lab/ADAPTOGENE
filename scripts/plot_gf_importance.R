library(qs)
library(ggplot2)
library(ggpubr)
library(tibble)
library(forcats)
library(gradientForest)
library(stringr)
library(svglite)
args = commandArgs(trailingOnly=TRUE)
###############
GF_ADAPTIVE_PATH = args[1]
GF_RANDOM_PATH = args[2]    # path or "NULL"
PLOT_DIR = args[3]
INTER_DIR = args[4]
SUFFIX = args[5]
###############

message('INFO: Plotting Overall Importance')

plot_theme <- theme_classic(base_size = 12, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 17, face = 'bold'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 14))

gf <- qread(GF_ADAPTIVE_PATH)
has_random <- GF_RANDOM_PATH != 'NULL'

# Adaptive importance
imp <- importance(gf, type = 'Weighted') %>%
  enframe() %>%
  dplyr::mutate(name = as.factor(name))

max_val <- max(imp$value) + 0.005

gAdapt <- ggplot(imp, aes(y = fct_reorder(name, value), x = value)) +
  geom_bar(stat = 'identity', fill = 'red', color = 'black') +
  plot_theme +
  labs(x = expression(paste("R"^2, " weighted importance")),
       title = 'Adaptive') +
  xlim(c(0, max_val))

if (has_random) {
  gf_random <- qread(GF_RANDOM_PATH)
  imp_random <- importance(gf_random, type = 'Weighted') %>%
    enframe() %>%
    dplyr::mutate(name = as.factor(name))

  gNeutral <- ggplot(imp_random, aes(y = fct_reorder(name, value), x = value)) +
    geom_bar(stat = 'identity', fill = 'blue', color = 'black') +
    plot_theme +
    labs(x = expression(paste("R"^2, " weighted importance")),
         title = 'Neutral') +
    xlim(c(0, max_val))

  gImp <- ggarrange(gAdapt, gNeutral, ncol = 2)
} else {
  gImp <- gAdapt
}

ggsave(paste0(PLOT_DIR, 'OverallImportance_', SUFFIX, '.png'), gImp)
ggsave(paste0(PLOT_DIR, 'OverallImportance_', SUFFIX, '.svg'), gImp,
       device = svglite::svglite, bg = "transparent", fix_text_size = FALSE)
qsave(gImp, paste0(INTER_DIR, 'OverallImportance_', SUFFIX, '.qs'))

message('INFO: Overall importance plot complete')
