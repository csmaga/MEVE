#!/bin/bash
#SBATCH --job-name=SNPGenie
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j
#SBATCH --mail-user=crs12448@uga.edu
#SBATCH --mail-type=END,FAIL


ml SNPGenie/1.0-foss-2022a-Perl-5.34.1

snpgenie.pl --vcfformat=4 --gtffile Amiss.annot.2022.gtf --fastafile Amiss_ref.fasta 

