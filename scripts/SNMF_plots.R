library(LEA)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
#library(gridExtra)
library(ggpubr)
library(scatterpie)
library(qs)
library(tibble)
args = commandArgs(trailingOnly=TRUE)
set.seed(42)
# Choose colors
my.colors = c("#f3c300", "#a1caf1", "#be0032", "#8db600", "#e68fac", "#0067a5", "#f38400", "#222222", 
  "#b3446c", "#2b3d26", "#dcd300",  "#875692", "#848482", "#c2b280", "#882d17", "#654522", "#e25822",
  "#008856", "#f99379", "#604e97", "#f6a600") # from Polychrome pallete

####################
SNMF_RES=args[1] ; snmf_res = load.snmfProject(SNMF_RES)
LFMM=args[2]
K_START=args[3] %>% as.numeric
K_END=args[4] %>% as.numeric
PLOIDY=args[5]
SAMPLES=args[6] # df of site and sample columns with header
####################
# Plots settings
plot_theme <-
        theme_classic(base_size=12, base_family = 'Helvetica') +
             theme(plot.title = element_text(size = 17, face = 'bold'),
                   axis.text = element_text(face = "bold", size = 18),
                   axis.title = element_text(face = "bold", size = 22),
                   legend.text = element_text(size = 18),
                   legend.title = element_text(face = "bold", size = 22)
                   )
################# Functions ##########
# FUN plot STRUCTURE plots
FUN_structure <- function(project, Kbest, samples.df, my.colors){
  # select the best run for K = 4
  best = which.min(cross.entropy(project, K = Kbest))
  Q.matrix <- as.matrix(Q(project, K = Kbest, run = best)) # replace the number of K
  colnames(Q.matrix) <- paste0(rep('C', Kbest), 1:Kbest)
  clusters <- cbind(samples.df, Q.matrix) %>% as.data.table %>% dplyr::mutate(sample = as.factor(sample))
  #TODO can add SITE option to dive sample on groups + maybe by the rule of naming site_sample
  # Make wide table
  IsrPopKs <- clusters %>% dplyr::select(-site) %>% melt.data.table(id.vars = 'sample')
  IsrPopKs %>% str
  # Plot
  gStructure <-
          ggplot(data = IsrPopKs, aes(y=value, x= sample, fill= variable)) +
            geom_bar(show.legend = T, stat="identity", position = "fill") +
            ylab("Proportion of assignment") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90)) +
            scale_fill_manual(values=my.colors) + #TODO change colors
            xlab("Accessions") +
          #   facet_nested(~ IsrPopKs$Pop, scales = "free", space = "fixed", nest_line = 1) + # add or remove "IsrPopKs$Region +"
            theme(strip.text = element_text(size=44, angle = 0),
                  axis.text.x = element_blank()) +  #change top strip font size & angle
            theme(axis.text=element_text(size=16),
                  axis.title=element_text(size=20,face="bold"),
                  legend.text = element_text(size = 16),
                  legend.title = element_text(size = 20, face='bold')) +
            theme(strip.text.x = element_text(size = 16)) +
            theme(strip.text.y = element_text(size = 16)) +
            guides(fill=guide_legend(title="Clusters"))
          #  theme(strip.background = element_rect(fill1 = "blue", fill2 = "red", fill3 = "green"))
  ggsave(paste0('plots/SNMF_structure_K', Kbest, '.png'), gStructure,  width = 3 * 6.4, height = 4.8)
  ggsave(paste0('plots/SNMF_structure_K', Kbest, '.svg'), gStructure,  width = 3 * 6.4, height = 4.8,
       device = svglite::svglite,
       bg = "transparent"
 	)  
  qsave(gStructure, paste0('intermediate/SNMF_structure_K', Kbest, '.qs'))

  return(clusters)
}

# Fun plot PCA and TW plot
FUN_pca <- function(LFMM, samples.df){
  # browser()
  # PCA
    pc <- LEA::pca(LFMM, scale = T)
  ## Plot
    ### Preprocess
    PCA <-
      as.data.frame(pc$projections) %>%
        setNames(paste0('PC', 1:ncol(.)))
    isrPCA <- cbind(samples.df, PCA)
    var.explained = c((pc$eigenvalues[1,] / sum(pc$eigenvalues)) * 100,
                      (pc$eigenvalues[2,] / sum(pc$eigenvalues)) * 100)

    gPCA <-
      ggplot(isrPCA, aes(x = PC1, y = PC2,
                               label = as.factor(site),
                               color = as.factor(site))) + # Here Colors depends on number of Pops which we can't control TODO
        #geom_point(size = 4, alpha = 0.6) +
	geom_label() +
        # geom_label_repel(max.overlaps = 9999) +
        xlab("PC1") +
        ylab("PC2") +
        guides(color = guide_legend(ncol=2)) +
        theme_bw()
        #theme(legend.position = 'none')
      ggsave(paste0("plots/PCA.png"), gPCA)
      ggsave(paste0("plots/PCA.svg"), gPCA, device = svglite::svglite, bg = 'transparent')
      qsave(gPCA, 'intermediate/PCA.qs')

  # Tracy widow
  tw <- tracy.widom(pc) # to estimate K (number of subpopulations, subistructures)
    ## Plot
    gTW <-
      ggplot(data = tw, aes(x= N, y= percentage))+
        geom_point()+
        xlim(1, 25)+
        theme_classic()+
        ggtitle("Tracy Widom") +
        theme(plot.title = element_text(hjust = 0.5))

    ggsave("plots/PCA_TracyWidom.png", gTW)
    ggsave("plots/PCA_TracyWidom.svg", gTW, device = svglite::svglite, bg = 'transparent')
    qsave(gTW, 'intermediate/PCA_TracyWidow.qs')

  return(pc)
}

# FUN plot structure PCA
FUN_structure_pca <- function(pc, # pc data object: pc$PCA$data
                              clusters, # site sampleName C1 C2 C3...
                              my.colors) {
  # browser()
  # Preprocess
  PCs.all <- as.data.frame(pc$projections) %>% setNames(gsub('V', 'PC', colnames(.)))
  df <- cbind(clusters, PCs.all[,1:4])
  var.explained = c((pc$eigenvalues[1,] / sum(pc$eigenvalues)) * 100,
                    (pc$eigenvalues[2,] / sum(pc$eigenvalues)) * 100)

  clustN = sum(grepl("^C", colnames(df))) # number of clusters
  # Plot PCA with SNMF clusters
  gPCA = ggplot(df, aes(PC1,PC2, label = as.factor(site) )) +
    geom_scatterpie(data = df,
                    aes(x = PC1, y = PC2
                        # , group = site
                        ),
                        cols = paste0('C', 1:clustN),
                        color = 'black', alpha = 0.8) +
    # geom_label_repel(max.overlaps = 9999) +
    scale_fill_manual(values = my.colors) +

    xlab(paste0('PC1 (', round(var.explained[1], 1), '%)')) +
    ylab(paste0('PC2 (', round(var.explained[2], 1), '%)')) +
    plot_theme
    #theme(legend.position = 'none')
 
 ggsave(paste0('plots/PCA_structurePie_K', clustN, '.png'), gPCA)
 ggsave(paste0('plots/PCA_structurePie_K', clustN, '.svg'), gPCA,
       device = svglite::svglite,
       bg = "transparent"
 	)  
 qsave(gPCA, paste0('intermediate/PCA_structurePie_K', clustN, '.qs'))

  return(gPCA)
}

FUN_diff_test <- function(project, ploidy, K){

    # Population differentiation tests
    p = snmf.pvalues(project,
                     entropy = TRUE,
                     ploidy = ploidy,
                     K = K)
    pvalues.df = data.frame(pvalues = p$pvalues)
    # Plot
    gHist <-
      ggplot(pvalues.df, aes(x=pvalues)) +
        geom_histogram(aes(y = ..density..), color = 'black', fill = 'lightblue') +
        geom_density(alpha = 0.2, fill = '#FF6666') +
        plot_theme

    gPval <-
      ggplot(pvalues.df, aes(y = -log10(pvalues), x = 1:nrow(pvalues.df))) +
        geom_point(color = 'blue') +
        labs(x = 'Index') +
        plot_theme


    gGrid <-
            ggarrange(gHist, gPval,
                           nrow = 2
                           #top = paste0('K', K)
	    )
    ggsave(paste0('plots/SNMF_PopDiffTest_K', K, '.png'), gGrid, width = 2 * 6.4, height = 9.6)
    ggsave(paste0('plots/SNMF_PopDiffTest_K', K, '.svg'), gGrid, width = 2 * 6.4, height = 9.6,
       device = svglite::svglite,
       bg = "transparent"
 	)  
 qsave(gPCA, paste0('intermediate/PCA_structurePie_K', clustN, '.qs'))

  return(gPCA)
}

FUN_diff_test <- function(project, ploidy, K){

    # Population differentiation tests
    p = snmf.pvalues(project,
                     entropy = TRUE,
                     ploidy = ploidy,
                     K = K)
    pvalues.df = data.frame(pvalues = p$pvalues)
    # Plot
    gHist <-
      ggplot(pvalues.df, aes(x=pvalues)) +
        geom_histogram(aes(y = ..density..), color = 'black', fill = 'lightblue') +
        geom_density(alpha = 0.2, fill = '#FF6666') +
        plot_theme

    gPval <-
      ggplot(pvalues.df, aes(y = -log10(pvalues), x = 1:nrow(pvalues.df))) +
        geom_point(color = 'blue') +
        labs(x = 'Index') +
        plot_theme


    gGrid <-
            ggarrange(gHist, gPval,
                           nrow = 2
                           #top = paste0('K', K)
	    )
    
    qsave(gGrid, paste0('intermediate/SNMF_PopDiffTest_K', K, '.qs'))
}
################# MAIN ################
							message('INFO: Load samples info')
samples.df <-
        fread(SAMPLES,
	      colClasses = c("site" = "character", 'sample' = 'character')) %>%
        dplyr::select(sample, site)

# Run pca
                                                        message('INFO: Run PCA')
pc <- FUN_pca(LFMM, samples.df)
# Save Clusters and STRUCTURE/PCA plots for each K
for(K in K_START:K_END){
        # Plot structure and calculate clusters
                                                        message(paste0("INFO: Run structure for K", K))
        clusters <- FUN_structure(snmf_res, K, samples.df, my.colors[1:K])
        # Save clusters info of SNMF
        write.table(clusters, paste0('tables/SNMF_clusters_K', K, '.tsv'),
                       sep='\t', row.names = F, quote = F)
        # Plot PCA structure
                                                        message(paste0("INFO: Run structure PCA for K", K))
        FUN_structure_pca(pc, clusters, my.colors[1:K])
							message(paste0("INFO: Run differentiation test for K", K))
	FUN_diff_test(snmf_res, PLOIDY, K)
}

							message(paste0("INFO: Run Cross Entropy test on K", K_START, '-', K_END))
gCrossEntropy <-
	sapply(K_START:K_END, function(K){min(cross.entropy(snmf_res, K))}) %>% 
	  setNames(K_START:K_END) %>%
	  enframe %>%
	  dplyr::arrange(as.numeric(name)) %>%
	    ggplot(aes(x = paste0('K', name),
	               y = value) ) +
	      geom_point() +
	  plot_theme +
	  labs(x = 'Number of ancestral populations',
	       y = 'Cross-entropy')

ggsave(paste0('plots/SNMF_Cross_entropy_K', K_START, '-', K_END, '.png'), gCrossEntropy)
ggsave(paste0('plots/SNMF_Cross_entropy_K', K_START, '-', K_END, '.svg'), gCrossEntropy,
	device = svglite::svglite,
	bg = 'transparent')
qsave(gCrossEntropy, paste0('intermediate/SNMF_Cross_entropy_K', K_START, '-', K_END, '.qs'))
