#!/bin/bash
#SBATCH --job-name=HISAT2	                        # Job name
#SBATCH --partition=highmem_p	                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=10	                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=250gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/HISAT2.o    # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=/scratch/crs12448/MEVE/Logs/HISAT2.e
#SBATCH --mail-user=christopher.smaga@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

DD="/scratch/crs12448/MEVE/Trimming/Trimmed_reads"

module load BEDTools/2.30.0-GCC-8.3.0
ml HISAT2/2.2.1-foss-2019b

cd /scratch/crs12448/MEVE/Alignment/HISAT2/Genome_Index

# this code is based on MDH's STIMRNAseq github
# using BEDTools to manually create exon and intron tables from the annotation file, which then can be used by hisat2 to build the genome index
# Basically, calling out all exons (from the third column) and their positions, then converting to BED format. Then all genes from the 3rd column are called out and put into BED format.
# To get intronic sights, exons are subtracted from gene coordinates.
#

# First, take the 1, 4, and 5 columns from the annotation file for which the 3rd column is an exon
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4,$5}' /home/crs12448/ALL_METH_PROJ/Amiss.annot.2022.gff > Amiss.ref_exon.bed
echo "Amiss.ref_exon.bed"
head -n 10 Amiss.ref_exon.bed

#Sort exons
sortBed -i Amiss.ref_exon.bed > Amiss.ref_exon_temp.bed

#Rename sorted file
mv -f Amiss.ref_exon_temp.bed Amiss.ref_exon.bed
echo "Amiss.ref_exon.bed"
head -n 10 Amiss.ref_exon.bed

#Merge exons together that overlap
mergeBed -i Amiss.ref_exon.bed > Amiss.ref_exon_merged.bed
echo "Amiss.ref_exon_merged.bed"
head -n 10 Amiss.ref_exon_merged.bed
#
# For lines where 3rd column is gene, extract columns 1,4,5,7 
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4,$5,$7}' /home/crs12448/ALL_METH_PROJ/Amiss.annot.2022.gff > Amiss.ref_gene.bed
echo "Amiss.ref_gene.bed"
head -n 10 Amiss.ref_gene.bed

#Sort
sortBed -i Amiss.ref_gene.bed > Amiss.ref_gene_temp.bed

#Rename
mv -f Amiss.ref_gene_temp.bed Amiss.ref_gene.bed
echo "Amiss.ref_gene.bed"
head -n 10 Amiss.ref.gene.bed

# Subtract exons from gene coordinates to get introns 
subtractBed -a Amiss.ref_gene.bed -b Amiss.ref_exon_merged.bed > Amiss.ref_intron.bed
head -n 10 Amiss.ref_intron.bed

#######################################################################################################################################

module load HISAT2/2.2.1-foss-2019b
# # Build index genome using splice sites and exon files created above
hisat2-build -p 10 --ss /scratch/crs12448/MEVE/Alignment/Genome_Index/Amiss.ref_intron.bed --exon /scratch/crs12448/MEVE/Alignment/Genome_Index/Amiss.ref_exon_merged.bed -f /home/crs12448/ALL_METH_PROJ/Amiss.ref.2022.fna /scratch/crs12448/MEVE/Alignment/Genome_Index/Index/Amiss.index.ref.hisat2

# For each file (S### is prefix), align the two reads resulting from that sample. Output file is SAM format.
cd $DD 
for file in  S392 S319	S280 S249 S388 S337	S343 S359 S391 S302	S263 S231 S316 S246	S432 S376 S406 S317	S266_2 S247 S393 S295 S344 S256_2
do
    echo "Working on $file."
    hisat2 -x /scratch/crs12448/MEVE/Alignment/HISAT2/Genome_Index/Index/Amiss.index.ref.hisat2  -p 10 --rna-strandness FR --dta -q -1 $DD/${file}_1_val_1.fq.gz -2 $DD/${file}_2_val_2.fq.gz -S /scratch/crs12448/MEVE/Alignment/HISAT2/SAM/${file}.sam --summary-file /scratch/crs12448/MEVE/Alignment/HISAT2/HISAT2_alignment_summary
    echo "$file finished."
done

#########################################################################################################################################

# Sort the SAM alignment files and convert to BAM format
ml SAMtools/1.14-GCC-8.3.0
for i in  S392 S319	S280 S249 S388 S337	S343 S359 S391 S302	S263 S231 S316 S246	S432 S376 S406 S317	S266_2 S247 S393 S295 S344 S256_2
do
echo "Working on $i."
samtools sort -@ 10 /scratch/crs12448/MEVE/Alignment/HISAT2/SAM/${i}.sam -o /scratch/crs12448/MEVE/Alignment/HISAT2/BAM/${i}.bam 
echo "$i finished."
done