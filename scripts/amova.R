library(dplyr)
library(data.table)
library(poppr)
library(vcfR)
library(adegenet)
library(qs)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
#################################
VCF = args[1]
SAMPLES = args[2]
CPU = args[3] %>% as.numeric
PLOT_DIR = args[4]
TABLES_DIR = args[5]
INTER_DIR = args[6]
#################################

RANDOM_N = 10000
set.seed(42)

#=============================================================================
# SINGLE-SITE GUARD
#=============================================================================
# AMOVA partitions variance BETWEEN and WITHIN populations. With one population
# there is no between-population stratum: df_between = n_pop - 1 = 0, the results
# table below indexes rows that do not exist, and randtest() has nothing to
# permute. Keying on the site count (not on coordinate spread) follows the same
# reasoning as MIN_SITES in pregea_dbmem.R.
#
# Write every declared output and exit 0 rather than stop(): population stats are
# an opt-in extra and must not take the whole structure mode down with them.
site_check <- fread(SAMPLES, colClasses = c("site" = "character", "sample" = "character"))
n_sites <- length(unique(site_check$site))
if (n_sites < 2) {
  message(paste0("WARNING: AMOVA skipped -- ", n_sites,
                 " sampling site(s) in the metadata, need >= 2. ",
                 "There is no between-population stratum to partition."))

  data.table(
    Parameter = c('status', 'n_sites',
                  'Between_populations_var_percent',
                  'Within_populations_var_percent',
                  'Phi-samples-total', 'p-value', 'Nperm',
                  'df_between', 'df_within', 'df_total'),
    Value = c('skipped_single_site', as.character(n_sites),
              rep(NA_character_, 8))
  ) %>%
    fwrite(paste0(TABLES_DIR, 'amova.tsv'), sep = '\t')

  skip_plot <- ggplot() +
    annotate("text", x = 0, y = 0, size = 5, colour = "grey50",
             label = paste0("AMOVA unavailable\n(", n_sites,
                            " sampling site, need >= 2)")) +
    theme_void()
  ggsave(paste0(PLOT_DIR, "amova.png"), skip_plot, width = 7, height = 7, dpi = 150)
  ggsave(paste0(PLOT_DIR, "amova.svg"), skip_plot, width = 7, height = 7,
         device = svglite::svglite, bg = 'transparent')
  qsave(list(AMOVA = NULL, AMOVAsignif = NULL, status = 'skipped_single_site'),
        paste0(INTER_DIR, 'amova_result.qs'))
  qsave(skip_plot, paste0(INTER_DIR, 'amova_plot.qs'))

  quit(status = 0)
}

# Read VCF file and convert to genind format
vcf <- vcfR::read.vcfR(VCF)

# Subset to random SNPs for speed
total_snps <- nrow(vcf@fix)
random_rows <- sample(1:total_snps, size = min(RANDOM_N, total_snps))
vcf <- vcf[random_rows, ]

my_genind <- vcfR2genind(vcf,
                         return.alleles = TRUE,
                         ploidy = 2)

my_genlight <- vcfR2genlight(vcf, n.cores = CPU)

# Read population information
samples.df <- fread(SAMPLES,
                    colClasses = c("site" = "character", 'sample' = 'character')) %>%
  dplyr::select(sample, site)

strata(my_genind) = samples.df %>%
  dplyr::mutate(pop = as.character(site)) %>%
  dplyr::select(sample, pop)
pop(my_genind) = samples.df$site

# Calculate genetic distance and perform AMOVA
gen.dist <- bitwise.dist(my_genind,
                         euclidean = TRUE,
                         threads = CPU)

AMOVA <- poppr.amova(my_genind,
                     ~pop,
                     dist = gen.dist,
                     within = FALSE,
                     squared = FALSE,
                     quiet = FALSE,
                     threads = CPU)

# Test significance
AMOVAsignif <- randtest(AMOVA, nrepet = 999)

# Create comprehensive AMOVA results table
data.table(
  Parameter = c('Between_populations_var_percent',
                'Within_populations_var_percent',
                'Phi-samples-total',
                'p-value',
                'Nperm',
                'df_between',
                'df_within',
                'df_total'),
  Value = c(AMOVA$componentsofcovariance$`%`[1],
            AMOVA$componentsofcovariance$`%`[2],
            AMOVA$statphi$Phi,
            AMOVAsignif$pvalue,
            AMOVAsignif$rep %>% as.numeric,
            AMOVA$results$Df[1],
            AMOVA$results$Df[2],
            AMOVA$results$Df[3])
) %>%
  fwrite(paste0(TABLES_DIR, 'amova.tsv'), sep = '\t')

# Create data frame for the pie chart
variations <- data.frame(
  Source = c("Between pop", "Within pop"),
  Percentage = AMOVA$componentsofcovariance$`%`[1:2]
)

# Create the pie chart
pie_chart <- ggplot(variations, aes(x = "", y = Percentage, fill = Source)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("#FF9999", "#66B2FF")) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "AMOVA Results",
    subtitle = paste("p-value =", AMOVAsignif$pvalue,
                     "\nPhi-samples-total =", round(AMOVA$statphi$Phi, 4))
  )

# Add percentage labels
pie_chart <- pie_chart +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_stack(vjust = 0.5))

# Save outputs
ggsave(paste0(PLOT_DIR, "amova.png"), pie_chart)
ggsave(paste0(PLOT_DIR, "amova.svg"), pie_chart,
       device = svglite::svglite, bg = 'transparent')

# Save AMOVA results as qs
AMOVA_result <- list(AMOVA = AMOVA,
                     AMOVAsignif = AMOVAsignif)
qsave(AMOVA_result, paste0(INTER_DIR, 'amova_result.qs'))
qsave(pie_chart, paste0(INTER_DIR, 'amova_plot.qs'))
