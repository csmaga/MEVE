library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)

setwd("/scratch/crs12448/MEVE/Alignment/HISAT2/BAM")

all_samples <- c("S231.bam", "S242.bam", "S246.bam", "S247.bam", "S252.bam", "S256_2.bam", "S263.bam",
              "S266_2.bam", "S280.bam", "S295.bam", "S302.bam", "S303.bam", "S314.bam", "S316.bam", "S317.bam", "S319.bam", "S328.bam", "S336.bam",
             "S337.bam",  "S338.bam", "S344.bam", "S345.bam", "S348.bam", "S350.bam", "S353.bam", "S357.bam", "S359.bam", "S367.bam", "S376.bam", "S380.bam", "S384.bam",
            "S388.bam", "S391.bam", "S392.bam", "S393.bam", "S406.bam", "S407.bam", "S408.bam", "S416.bam", "S420.bam", "S421.bam", "S422.bam", "S425.bam",
           "S426.bam", "S427.bam", "S432.bam", "S433.bam", "S435.bam")
summary(all_samples)

BAM_files<- BamFileList(all_samples)

txdb <- makeTxDbFromGFF("/scratch/crs12448/GATOR_GENOME/Amiss.annot.2022.gff", circ_seqs = character())
txdb

flattenedAnnotation = exonicParts( txdb, linked.to.single.gene.only=TRUE )
names(flattenedAnnotation) =
  sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)
          
          
exon_counts <- summarizeOverlaps(features=flattenedAnnotation, reads=BAM_files, singleEnd=FALSE, ignore.strand=FALSE, fragments=TRUE)
exon_counts

save.image(file="/scratch/crs12448/MEVE/Splicing/MEVE_gonad_splicing.RData")

          