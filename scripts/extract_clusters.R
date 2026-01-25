library(LEA)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
#################################
SNMF_PROJECT = args[1]
SAMPLES = args[2]
K = args[3] %>% as.numeric
OUTPUT = args[4]
#################################

# Load SNMF project
snmf_res <- load.snmfProject(SNMF_PROJECT)

# Load samples info
samples.df <- fread(SAMPLES,
                    colClasses = c("site" = "character", 'sample' = 'character')) %>%
              dplyr::select(sample, site)

# Select the best run for K (lowest cross-entropy)
best <- which.min(cross.entropy(snmf_res, K = K))

# Extract Q-matrix
Q.matrix <- as.matrix(Q(snmf_res, K = K, run = best))
colnames(Q.matrix) <- paste0('C', 1:K)

# Combine with sample info
clusters <- cbind(samples.df, Q.matrix) %>%
            as.data.table %>%
            dplyr::mutate(sample = as.factor(sample))

# Save
write.table(clusters, OUTPUT, sep = '\t', row.names = FALSE, quote = FALSE)
