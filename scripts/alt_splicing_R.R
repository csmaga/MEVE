library(DEXSeq)
library(BiocParallel)

load("/scratch/crs12448/MEVE/Splicing/MEVE_gonad_splicing.RData")

se<-exon_counts

metadata<-read.csv("/scratch/crs12448/MEVE/Splicing/MEVE_RNASEQ_META_2.csv")

colData(se)$condition<-metadata$temp
colData(se)$pop<-metadata$pop

dxd = DEXSeqDataSetFromSE(se, design= ~ sample + exon + condition:exon)
colData(dxd)
head( counts(dxd), 5 )
dxd = estimateSizeFactors(dxd)

BPPARAM = MultiCoreParam(8)
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM )

dxd = testForDEU( dxd,BPPARAM=BPPARAM )
dxd = estimateExonFoldChanges( dxd,BPPARAM=BPPARAM, fitExpToVar="condition")


dxr1 = DEXSeqResults( dxd )
dxr1

dxr1.df<-as.data.frame(dxr1)
dxr1.df$exon<-rownames(dxr1.df)

save.image("/scratch/crs12448/MEVE/Splicing/analysis_results.RData")
