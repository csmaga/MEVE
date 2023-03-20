#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-user=christopher.smaga@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=MEVE_trimming
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/crs12448/MEVE/Logs/trim_QC.o
#SBATCH --error=/scratch/crs12448/MEVE/Logs/trim_QC.e


# DD="/scratch/crs12448/MEVE/Data/Raw"
# QC_OD="/scratch/crs12448/MEVE/QC/Trimmed_QC"
# Trim_OD="/scratch/crs12448/MEVE/Trimming/Trimmed_reads"


# module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4

# # Make the current directory the folder with raw data
# cd $DD

# # Create a for loop that goes into each directory, prints that directory, then does trimming on the two reads in that directory and fastqc, outputing fastqc results into QC_DD and trimmed reads into Trim_OD
# # Trim adapter sequences from fastq files, setting strigency to 3 (only trimmed if at least 3 bases match adaptor sequence), which is what was used for STIM2014, SLB's work, and PREE2. 

# for dir in $DD/*; do
#   echo "Starting $dir"
#   cd $dir
#   trim_galore --cores 4 --fastqc --fastqc_args "--outdir $QC_OD" -stringency 3 -o $Trim_OD --paired  *_1.fq.gz *_2.fq.gz
#   echo "Trimming of $dir finished."
# done

module load FastQC/0.11.9-Java-11
module load MultiQC/1.8-foss-2019b-Python-3.7.4

cd $QC_OD
multiqc /scratch/crs12448/MEVE/QC/Trimmed_QC

