#!/bin/bash
#SBATCH --job-name=GATK_SNP                 # Job name
#SBATCH --partition=batch                           # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=8                    # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=48gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/BayeScan.o    # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=/scratch/crs12448/MEVE/Logs/Bayescan.e
#SBATCH --mail-user=christopher.smaga@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

ml  BayeScan/2.1-foss-2019b
cd /scratch/crs12448/MEVE/Bayescan


bayescan AP_WO_SNPS_BAYESCAN.bayescan -snp \
    -od Bayescan_out   \
    -threads 8 \
    -n 5000 \
    -thin 10    \
    -nbp 20 \
    -pilot 5000 \
    -burn 50000 \
    -pr_odds 1000   \
    -out_freq   



