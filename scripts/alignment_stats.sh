#!/bin/bash
#SBATCH --job-name=GAlignment_stats		                        # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=1		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=32gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j	    
#SBATCH --mail-user=christopher.smaga@uga.edu               
#SBATCH --mail-type=END,FAIL                            

cd "/scratch/crs12448/MEVE/Alignment/HISAT2/BAM"

ml  SAMtools/1.16.1-GCC-11.3.0

for i in *.bam
do
echo i
samtools flagstat i > /scratch/crs12448/MEVE/Alignment/Stats/{i}_stats
done
