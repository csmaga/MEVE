#!/bin/bash
#SBATCH --job-name=MEVE_QC		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6	                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=60gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/QC.o    # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=/scratch/crs12448/MEVE/Logs/QC.e
#SBATCH --mail-user=christopher.smaga@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

# Set data directory and output directory
DD="/scratch/crs12448/MEVE/Data/Raw"
OD="/scratch/crs12448/MEVE/QC"

#Load FastQC and multiQC
ml FastQC/0.11.9-Java-11
ml multiqc/1.11-GCCcore-8.3.0-Python-3.8.2

# Change directory to the data directory
cd $OD

# Create a for loop that goes into each directory, using the wildcard *, and prints the name of the directoy, then runs FastQC on each .gz file in that directory using 6 cores, then prints "Done"
for dir in $DD/*; do
  cd dir
  echo "$dir"
  fastqc -o $OD -t 6 *.gz
  echo "$dir done"
done

# Use multiQC to put all QC reports together into one .html file
cd $OD
multiqc $OD

module unload FastQC
module unload MultiQC

