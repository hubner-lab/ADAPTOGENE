library(ggplot2)
library(data.table)
library(dplyr)
library(qs)
library(ggcorrplot)

args = commandArgs(trailingOnly=TRUE)
#################################
CLIMATE = args[1]    # climate factors from sites
SAMPLES = args[2]    # samples with traits if exists
OUT_PNG = args[3]
INTER_DIR = args[4]
#################################

# Derive SVG and QS paths from PNG path
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)
OUT_QS  <- paste0(INTER_DIR, sub("\\.png$", ".qs", basename(OUT_PNG)))

# Plot theme
plot_theme <-
  theme_classic(base_size = 8, base_family = 'Helvetica') +
  theme(plot.title = element_text(size = 14, face = 'bold'),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = 'bold', size = 10, angle = -90),
        axis.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 12))

#################################

samples = fread(SAMPLES, colClasses = c("site" = "character", "sample" = "character"))
climate = fread(CLIMATE) %>% dplyr::select(-any_of('sample'))

# Combine climate with phenotypic traits if available
if (ncol(samples) > 4) {
  traits = cbind(climate,
                 samples %>% dplyr::select(-site, -sample, -latitude, -longitude))
} else {
  traits = climate
}

# Remove constant columns (filter NA from sd on non-numeric columns)
std_devs <- apply(traits, 2, sd, na.rm = TRUE)
non_constant_traits = names(std_devs)[!is.na(std_devs) & std_devs != 0]
traits <- traits %>% dplyr::select(all_of(non_constant_traits))

# Compute correlation matrix.
#
# BLOCK-WISE WEIGHTING. `climate` has one row per SAMPLE but its values are a
# property of the SITE (identical for every sample sequenced there), while the
# phenotype columns from `samples` are genuinely per-sample. Correlating the whole
# table over sample rows would weight each site by its sample size, so a site with
# 30 individuals counts 30x against a site with 1 in every climate-vs-climate |r|
# a user reads off this figure when curating predictors — and those |r| would then
# disagree with the site-weighted collinearity screen in pregea_rda_setup.R.
# Each pair is therefore correlated at its natural unit: climate-vs-climate over
# distinct SITES, anything involving a phenotype over SAMPLES. (Deduplicating the
# whole table instead would be a silent no-op the moment phenotype columns exist,
# since those vary within a site.)
climate_cols <- base::intersect(colnames(climate), colnames(traits))
# The site vector is taken positionally from `samples`, which the calling rule
# (climate.smk) guarantees is row-aligned to `climate` — a guarantee that has
# been broken once before. Fail loudly rather than dedupe against the wrong sites.
stopifnot(nrow(samples) == nrow(climate))
site_rows    <- !duplicated(samples$site)

traits_df  <- as.data.frame(traits)   # base-R `[` semantics: traits is a data.table,
                                     # whose `[` would treat a character j as an
                                     # expression, not a column selection
cor_matrix <- cor(traits_df, use = 'pairwise.complete.obs')
if (length(climate_cols) > 1 && sum(site_rows) >= 3) {
  cor_matrix[climate_cols, climate_cols] <-
    cor(traits_df[site_rows, climate_cols, drop = FALSE], use = 'pairwise.complete.obs')
  message('INFO: climate-vs-climate correlations computed over ', sum(site_rows),
          ' distinct sites (', nrow(traits_df), ' samples); phenotype pairs over samples')
} else if (length(climate_cols) > 1) {
  message('WARNING: only ', sum(site_rows), ' site(s) — climate correlations left ',
          'sample-weighted (a site-level correlation needs >= 3 sites)')
}

# Build the heatmap using ggcorrplot
gHM <- ggcorrplot(cor_matrix,
                  type = "lower",
                  lab = TRUE,
                  digits = 1,
                  lab_size = 4) +
  ggtitle('Correlogram of traits') +
  plot_theme

# Determine size based on number of traits
scale_factor <- 0.7
width <- ncol(traits) * scale_factor
height <- ncol(traits) * scale_factor

# Save
ggsave(OUT_PNG, gHM, width = width, height = height, units = "in")
ggsave(OUT_SVG, gHM, width = width, height = height, units = "in",
       device = svglite::svglite, bg = 'transparent')
qsave(gHM, OUT_QS)
