#!/bin/bash
#SBATCH --job-name=MEVE_alignment
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH --time=72:00:00
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j
#SBATCH --mail-user=crs12448@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-50 

sample=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/crs12448/MEVE/Data/Raw/sample_names)

## make project directory + make directory for ref genome
OUTDIR="/scratch/crs12448/MEVE"

# This is based off of the GATK "Best Practices" for RNAseq short variant discovery (SNPs + Indels). Found here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-

# First, mark duplicate reads from the bam files.

ml Java/17.0.6
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17

# Run MarkDuplicates to mark duplicate reads from RNAseq data. For each file in the DD, perform the mark duplicates function. Calling on JAVA with 64GB and 8 cores

mkdir $OUTDIR/GATK/MarkDuplicates

 gatk -jar picard.jar MarkDuplicates \
     I=$OUTDIR/Alignment/HISAT2/BAM/${sample}.bam \
     O=$OUTDIR/GATK/MarkDuplicates/${sample}_mark_dup.bam} \
     M=$OUTDIR/GATK/MarkDuplicates/${sample}_mark_dup_metrics.txt}

# SplitNCigar reads to split reads at "N"s (new read to the left and to the right of the N). Not sure exactly why this is necessary....
#  mkdir $OUTDIR/GATK/SplitNCigarReads

#  gatk SplitNCigarReads \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    -I $OUTDIR/GATK/MarkDuplicates/${sample}_mark_dup.bam  \
#    -O $OUTDIR/GATK/SplitNCigarReads/${sample}_cigar.bam

# # Recalibration of base quality
# mkdir $OUTDIR/GATK/Recalibration

# gatk BaseRecalibrator \
#    -I $OUTDIR/GATK/SplitNCigarReads/${sample}_cigar.bam \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    -O $OUTDIR/GATK/Recalibration/recal_data.table

# gatk ApplyBQSR \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    -I $OUTDIR/GATK/SplitNCigarReads/${sample}_cigar.bam \
#    --bqsr-recal-file $OUTDIR/GATK/Recalibration/recal_data.table \
#    -O $OUTDIR/GATK/Recalibration/${sample}_recal.bam
 
# gatk AnalyzeCovariates \
#   -bqsr $OUTDIR/GATK/Recalibration/recal_data.table \
#   -plots AnalyzeCovariates.pdf

# # Variant calling 

# #Make directory for vcf files
# mkdir mkdir $OUTDIR/GATK/HaplotypeCaller

# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#   -I $OUTDIR/GATK/Recalibration/${sample}_recal.bam \
#   -O $OUTDIR/GATK/HaplotypeCaller/${sample}.g.vcf.gz \
#   -ERC GVCF


