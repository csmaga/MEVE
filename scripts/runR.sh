#!/bin/bash
#SBATCH --job-name=MEVE_count_reads		                        # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=1		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=64gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j	    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=christopher.smaga@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

cd /scratch/crs12448/MEVE/R

ml R/4.3.1-foss-2022a

export R_LIBS=/home/crs12448/R/x86_64-pc-linux-gnu-library/4.3

R CMD BATCH ~/MEVE/scripts/GenomicFeatures_count_reads.R
