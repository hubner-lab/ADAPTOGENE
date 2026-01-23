library(dplyr)
library(data.table)
library(poppr)
library(vcfR)
library(adegenet)
library(qs)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
#################################
VCF = args[1] #TODO result for imputed and not imputed very differ, so it's better to return for non-imputed vcf as input
SAMPLES = args[2]
CPU = args[3]
RANDOM_N = 10000

#################################
set.seed(42)

# Read VCF file and convert to genind format
vcf <- vcfR::read.vcfR(VCF)

# Subset
total_snps <- nrow(vcf@fix)
random_rows <- sample(1:total_snps, size = min(RANDOM_N, total_snps))
vcf <- vcf[random_rows, ]

my_genind <- vcfR2genind(vcf,
                   return.alleles = TRUE,
		   ploidy = 2)
my_genlight <- vcfR2genlight(vcf, n.cores = CPU)

# Read population information from a file
samples.df <-
        fread(SAMPLES,
              colClasses = c("site" = "character", 'sample' = 'character')) %>%
        dplyr::select(sample, site)

strata(my_genind) = samples.df %>% dplyr::mutate(pop = as.character(site)) %>% dplyr::select(sample, pop)
pop(my_genind) = samples.df$site

# Calculate genetic distance and perform AMOVA
gen.dist <- bitwise.dist(my_genind,
             euclidean = T,
             threads = CPU)

AMOVA <- poppr.amova(my_genind,
                     ~pop,
                     ,dist = gen.dist
                     ,within = FALSE,
                     ,squared=FALSE,
                     quiet = FALSE,
                     threads = CPU)

# Test significants
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
  fwrite('tables/AMOVA.tsv', sep = '\t')

# Create data frame for the pie chart
variations <- data.frame(
  Source = c("Between pop", "Within pop"),
  Percentage = AMOVA$componentsofcovariance$`%`[1:2]
)

# Create the main pie chart
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
    subtitle = paste("p-value =",
                     AMOVAsignif$pvalue,
                     "\nPhi-samples-total =",
                     round(AMOVA$statphi$Phi, 4))
  )

# Add percentage labels
pie_chart <-
  pie_chart +
   geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_stack(vjust = 0.5))

# Saving results
ggsave("plots/AMOVA.png", pie_chart)

# Save AMOVA and AMOVAsignif
AMOVA_result <-
  list(AMOVA = AMOVA,
     AMOVAsignif = AMOVAsignif)
qsave(AMOVA_result, 'intermediate/AMOVA_result.qs')
qsave(pie_chart, 'intermediate/AMOVA_plot.qs')
