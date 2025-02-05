#!/bin/bash
#SBATCH --job-name=bam_subset
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j
#SBATCH --mail-user=crs12448@uga.edu
#SBATCH --mail-type=END,FAIL


cd /scratch/crs12448/MEVE/Alignment/HISAT2/BAM
ml SAMtools/1.18-GCC-12.3.0

# for file in *.bam; do
#   samtools view ${file} "NW_017708899.1:4382450-4395344" > /scratch/crs12448/MEVE/Alignment/HISAT2/BAM/gene_bams/aig1_${file}
# done



for file in *.bam; do
  samtools view ${file} "NW_017713108.1:30798-33207" > /scratch/crs12448/MEVE/Alignment/HISAT2/BAM/gene_bams/cyp1b1/cyp1b1_${file}
done