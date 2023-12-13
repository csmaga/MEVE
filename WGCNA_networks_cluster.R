
BiocManager::install("WGCNA")

library(WGCNA)
options(stringsAsFactors = FALSE)
#############

################ FOR AP vs WO ############################

#############
# # Set working directory to path where RData workspaces are stored 
setwd("/scratch/crs12448/AP_WO_21/WGCNA")
# 
# # Load the workspace. This file ONLY contains the normalized_counts and gene_weights variables, which are needed for running the analysis
load("WGCNA_vars.RData")
# 
# # Now create several networks, using different parameters. Keep the same random seed for easy comparison of methods. For all, included gene weights, using bicor for correlation (see RStudio script for details), and a signed network. ''
# MergeCutHeight is kept at 0.15 which seems to be good based on Peter (author) from online. Thresholding power of 8 was found good so using that. 

# # # This first run is with default parameters except those above, with min module size set to 20 and power set to 8. 
# network1_20_2<-blockwiseModules(normalized_counts, weights=AP_WO_gene_weights_WGCNA, maxBlockSize = 17612, randomSeed=12345,
#                            power = 8, minModuleSize = 20, corType = "bicor", maxPOutliers = 0.05,
#                            reassignThreshold = 1e-6, mergeCutHeight = 0.15,
#                            numericLabels = FALSE, networkType = "signed",
#                            saveTOMs = TRUE, deepSplit = 2,
#                            saveTOMFileBase = "network1TOM",
#                            pamStage=TRUE,
#                            pamRespectDendro=TRUE, nThreads=6,
#                            verbose = 5)
# # 
# # # Make more sensitive by increasing deepSplit to 4 (max sensitivity)
# network2_20_4<-blockwiseModules(normalized_counts, weights=AP_WO_gene_weights_WGCNA, maxBlockSize = 17612, randomSeed=12345,
#                            power = 8, minModuleSize = 20, corType = "bicor", maxPOutliers = 0.05,
#                            reassignThreshold = 1e-6, mergeCutHeight = 0.15,
#                            numericLabels = FALSE, networkType = "signed",
#                            saveTOMs = TRUE, deepSplit = 4,
#                            saveTOMFileBase = "network2TOM",
#                            pamStage=TRUE,
#                            pamRespectDendro=TRUE, nThreads=6,
#                            verbose = 5)
# 
# # # Set power to 8, deepSplit 2, min module size to 40
# network3_40_2<-blockwiseModules(normalized_counts, weights=AP_WO_gene_weights_WGCNA, maxBlockSize = 17612, randomSeed=12345,
#                            power = 8, minModuleSize = 40, corType = "bicor", maxPOutliers = 0.05,
#                            reassignThreshold = 1e-6, mergeCutHeight = 0.15,
#                            numericLabels = FALSE, networkType = "signed",
#                            saveTOMs = TRUE, deepSplit = 2,
#                            saveTOMFileBase = "network3TOM",
#                            pamStage=TRUE,
#                            pamRespectDendro=TRUE, nThreads=6,
#                            verbose = 5)
# # 
# # # Set power 8, deepSplit 4, min module size 40
# network4_40_4<-blockwiseModules(normalized_counts, weights=AP_WO_gene_weights_WGCNA, maxBlockSize = 17612, randomSeed=12345,
#                            power = 8, minModuleSize = 40, corType = "bicor", maxPOutliers = 0.05,
#                            reassignThreshold = 1e-6, mergeCutHeight = 0.15,
#                            numericLabels = FALSE, networkType = "signed",
#                            saveTOMs = TRUE, deepSplit = 4,
#                            saveTOMFileBase = "network3TOM",
#                            pamStage=TRUE,
#                            pamRespectDendro=TRUE, nThreads=6,
#                            verbose = 5)
# 
# # # Set power to 8, deepSplit 2, min module size to 60
# network3_60_2<-blockwiseModules(normalized_counts, weights=AP_WO_gene_weights_WGCNA, maxBlockSize = 17612, randomSeed=12345,
#                            power = 8, minModuleSize = 60, corType = "bicor", maxPOutliers = 0.05,
#                            reassignThreshold = 1e-6, mergeCutHeight = 0.15,
#                            numericLabels = FALSE, networkType = "signed",
#                            saveTOMs = TRUE, deepSplit = 2,
#                            saveTOMFileBase = "network5TOM",
#                            pamStage=TRUE,
#                            pamRespectDendro=TRUE, nThreads=6,
#                            verbose = 5)
# # 
# # # Set power 8, deepSplit 4, min module size 60
# network4_60_4<-blockwiseModules(normalized_counts, weights=AP_WO_gene_weights_WGCNA, maxBlockSize = 17612, randomSeed=12345,
#                            power = 8, minModuleSize = 60, corType = "bicor", maxPOutliers = 0.05,
#                            reassignThreshold = 1e-6, mergeCutHeight = 0.15,
#                            numericLabels = FALSE, networkType = "signed",
#                            saveTOMs = TRUE, deepSplit = 4,
#                            saveTOMFileBase = "network6TOM",
#                            pamStage=TRUE,
#                            pamRespectDendro=TRUE, nThreads=6,
#                            verbose = 5)
# # 
# # 
#save.image(file="WGCNA_networks.RData")

#Site-specific co-expression groups
# AP
AP_normalized_counts<-normalized_counts[c(1,2,5,7,9,11,14,15,17,19,20,21),]
AP_weights<-AP_WO_gene_weights_WGCNA[c(1,2,5,7,9,11,14,15,17,19,20,21),]

AP_network<-blockwiseModules(AP_normalized_counts, weights=AP_weights, maxBlockSize = 17612, randomSeed=12345,
                           power = 8, minModuleSize = 20, corType = "bicor", maxPOutliers = 0.05,
                           reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                           numericLabels = FALSE, networkType = "signed",
                           saveTOMs = TRUE, deepSplit = 2,
                           saveTOMFileBase = "network_AP_TOM",
                           pamStage=TRUE,
                           pamRespectDendro=TRUE, nThreads=6,
                           verbose = 5)
# WO
WO_normalized_counts<-normalized_counts[-c(1,2,5,7,9,11,14,15,17,19,20,21),]
WO_weights<-AP_WO_gene_weights_WGCNA[-c(1,2,5,7,9,11,14,15,17,19,20,21),]

WO_network<-blockwiseModules(WO_normalized_counts, weights=WO_weights, maxBlockSize = 17612, randomSeed=12345,
                             power = 8, minModuleSize = 20, corType = "bicor", maxPOutliers = 0.05,
                             reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                             numericLabels = FALSE, networkType = "signed",
                             saveTOMs = TRUE, deepSplit = 2,
                             saveTOMFileBase = "network_WO_TOM",
                             pamStage=TRUE,
                             pamRespectDendro=TRUE, nThreads=6,
                             verbose = 5)

save.image(file="WGCNA_networks_AP_WO.RData")