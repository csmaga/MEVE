#!/bin/bash
#SBATCH --job-name=Pixy
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j
#SBATCH --mail-user=crs12448@uga.edu
#SBATCH --mail-type=END,FAIL

ml  Miniconda3/23.5.2-0
ml VCFtools
source /apps/eb/Miniconda3/23.5.2-0/etc/profile.d/conda.sh

conda init bash
conda activate /home/crs12448/conda_env

pixy --stats pi fst dxy \
--vcf /scratch/crs12448/MEVE/GATK/GenotypeGVCFs/Filtered/MEVE_sites_pixy.vcf.gz \
--populations /scratch/crs12448/MEVE/PopGen/MEVE_pops.txt \
--n_cores 12 \
--bed_file /scratch/crs12448/MEVE/PopGen/gene_bed_sorted.bed \
--output_folder /scratch/crs12448/MEVE/PopGen \
