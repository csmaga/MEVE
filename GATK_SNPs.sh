#!/bin/bash
#SBATCH --job-name=GATK_SNP                 # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=8                           # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=48gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/GATK_prep.o    # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=/scratch/crs12448/MEVE/Logs/GATK_prep.e
#SBATCH --mail-user=christopher.smaga@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

# Directorty where aligned and sorted BAM files are located
DD="/scratch/crs12448/MEVE/Alignment/HISAT2/BAM"

# Set ouput directory
OD="/scratch/crs12448/MEVE/GATK/MarkDuplicates"

#Load the Genome Anlysis Toolkit
ml  GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8

# Run MarkDuplicates to mark duplicate reads from RNAseq data. For each file in the DD, perform the mark duplicates function
cd $DD
for i in *.bam;
do
gatk --java-options "-Xmx48G -XX:+UseParallelGC -XX:ParallelGCThreads=8" MarkDuplicates \
      I=$i \
      O=$OD/${i/.bam/}_mark_dup.bam \
      M=$OD/${i.bam/}_mark_dup_metrics.txt
done

############################################################################################################################################3

cd $OD
OD_2="/scratch/crs12448/MEVE/GATK/SplitNCigarReads"

# Run SplitNCigarReads to split reads that span introns into separate reads
for i in *.bam;
do
 gatk SplitNCigarReads \
      -R /home/crs12448/ALL_METH_PROJ/Amiss.ref.2022.fna \
      -I $i \
      -O $OD_2/${i/_mark_dup.bam/_cigar.bam}
done

############################################################################################################################################

# Recalibrate base quality score
 gatk BaseRecalibrator \
   -I my_reads.bam \
   -R reference.fasta \
   --known-sites sites_of_variation.vcf \
   --known-sites another/optional/setOfSitesToMask.vcf \
   -O recal_data.table