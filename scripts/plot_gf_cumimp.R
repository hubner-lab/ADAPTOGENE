library(qs)
library(ggplot2)
library(ggpubr)
library(gradientForest)
library(stringr)
library(svglite)
args = commandArgs(trailingOnly=TRUE)
###############
GF_ADAPTIVE_PATH = args[1]
GF_RANDOM_PATH = args[2]    # path or "NULL"
PREDICTORS_SELECTED = args[3] %>% str_split(',') %>% unlist
PLOT_DIR = args[4]
INTER_DIR = args[5]
SUFFIX = args[6]
###############

message('INFO: Plotting Cumulative Importance')

plot_theme <- theme_classic(base_size = 12, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 17, face = 'bold'),
        axis.text = element_text(face = NULL, size = 12),
        axis.title.y = element_text(face = NULL, size = 12),
        axis.title.x = element_text(face = NULL, size = 16),
        legend.title = element_blank())

gf <- qread(GF_ADAPTIVE_PATH)
has_random <- GF_RANDOM_PATH != 'NULL'
if (has_random) {
  gf_random <- qread(GF_RANDOM_PATH)
}

gList <- lapply(PREDICTORS_SELECTED, function(bio) {
  gf_cumimp <- cumimp(gf, predictor = bio, type = "Overall")
  df <- data.frame(x = gf_cumimp$x,
                   y = gf_cumimp$y,
                   Model = 'Adaptive')

  if (has_random) {
    gf_random_cumimp <- cumimp(gf_random, predictor = bio, type = "Overall")
    df <- rbind(df,
                data.frame(x = gf_random_cumimp$x,
                           y = gf_random_cumimp$y,
                           Model = 'Neutral'))
  }

  gPlot <- ggplot(df, aes(x = x, y = y, color = Model)) +
    geom_line() +
    labs(x = bio, y = 'Cumulative importance') +
    scale_color_manual(values = c('Adaptive' = 'red', 'Neutral' = 'blue')) +
    plot_theme

  return(gPlot)
})

# Determine grid layout
n <- length(gList)
ncols <- min(n, 3)
nrows <- ceiling(n / ncols)

gCumImp <- ggarrange(plotlist = gList, nrow = nrows, ncol = ncols, common.legend = TRUE)

ggsave(paste0(PLOT_DIR, 'CumulativeImportance_', SUFFIX, '.png'), gCumImp)
ggsave(paste0(PLOT_DIR, 'CumulativeImportance_', SUFFIX, '.svg'), gCumImp,
       device = svglite::svglite, bg = "transparent", fix_text_size = FALSE)
qsave(gCumImp, paste0(INTER_DIR, 'CumulativeImportance_', SUFFIX, '.qs'))

message('INFO: Cumulative importance plot complete')
