 # if (!require("BiocManager", quietly = TRUE))
 #     install.packages("BiocManager")
 # 
 # BiocManager::install("GenomicFeatures")
 # if (!require("BiocManager", quietly = TRUE))
 #     install.packages("BiocManager")
 # 
 # BiocManager::install("GenomicAlignments")


library(GenomicAlignments)
library(GenomicFeatures)

# Change directory to where the BAM files are 

setwd("/scratch/crs12448/MEVE/TimeSeries/Alignment/HISAT2/BAM")

all_samples <- c("DRR048609.bam", "DRR048610.bam", "DRR048611.bam", "DRR048612.bam ", "DRR048613.bam", "DRR048614.bam", "DRR048615.bam", " DRR048616.bam", "  DRR048617.bam ",
               "DRR048618.bam ", "DRR048619.bam", "DRR048620.bam", "DRR048621.bam", "DRR048622.bam", "DRR048623.bam", "DRR048624.bam", "DRR048625.bam", "DRR048626.bam", "DRR048627.bam", "DRR048628.bam",
              "DRR048629.bam",  "DRR048630.bam", "DRR048631.bam", "DRR048632.bam", "DRR048633.bam", "DRR048634.bam", "DRR048635.bam", "DRR048636.bam")
summary(all_samples)

# Create a new object for BAM files in a single list
BAM_files<- BamFileList(all_samples)
#
# # Create a storage object for GTF annotations
txdb <- makeTxDbFromGFF("/scratch/crs12448/MEVE/Genome/Amiss.annot.2022.gff", circ_seqs = character())
txdb
#keytypes(txdb)

#
# # Extract all exons grouped within genes - code from MDH and SLB
exons_by_genes <- exonsBy(txdb, by="gene")
exons_by_genes
#write.csv(exons_by_genes, "exons_by_genes.csv")

#transcript_lengths <- transcripts(txdb)
#transcript_lengths

gene_counts <- summarizeOverlaps(features=exons_by_genes, reads=BAM_files, mode="Union", singleEnd=FALSE, ignore.strand=FALSE, fragments=TRUE)
gene_counts

counts_matrix <- assays(gene_counts)
counts_annotations <- rowRanges(gene_counts)

dim(counts_matrix)
dim(counts_annotations)

save.image()

write.csv(counts_matrix, '/scratch/crs12448/MEVE/TimeSeries/ReadCounts/time_series_read_counts.csv')
write.csv(counts_annotations, '/scratch/crs12448/MEVE/TimeSeries/ReadCounts/time_series_read_counts_annotation.csv')
