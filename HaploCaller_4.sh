#!/bin/bash
#SBATCH --job-name=GATK_SNP                 # Job name
#SBATCH --partition=highmem_p	                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=4                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/GATK_haplo_4.o    # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=/scratch/crs12448/MEVE/Logs/GATK_haplo_4.e
#SBATCH --mail-user=christopher.smaga@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

OD_4="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/GVCF"
cd /scratch/crs12448/MEVE/GATK/FixBam

for i in S388 S391 S392 S393 S406 S432;
do
 gatk --java-options "-Xmx120g -XX:+UseParallelGC -XX:ParallelGCThreads=4" HaplotypeCaller  \
   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
   -I ${i}_cigar_fix.bam \
   -O $OD_4/${i/_cigar_fix.bam/.g.vcf.gz} \
   -ERC GVCF
done
