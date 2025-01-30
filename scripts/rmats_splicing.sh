#!/bin/bash
#SBATCH --job-name=MEVE_splicing
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j
#SBATCH --mail-user=crs12448@uga.edu
#SBATCH --mail-type=END,FAIL

ml rMATS-turbo/4.1.1-foss-2020b

rmats.py --b1 /scratch/crs12448/MEVE/rmats/southern_pops.txt --b2 /scratch/crs12448/MEVE/rmats/northern_pops.txt --gtf /scratch/crs12448/GATOR_GENOME/Amiss.annot.2022.gff -t paired --readLength 150 --libType fr-firststrand --od /scratch/crs12448/MEVE/rmats/results --tmp /scratch/crs12448/MEVE/rmats/temp