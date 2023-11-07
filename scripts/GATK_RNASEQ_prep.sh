#!/bin/bash
#SBATCH --job-name=MEVE_GATK
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH --time=72:00:00
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j
#SBATCH --mail-user=crs12448@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-50 

sample=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/crs12448/MEVE/Data/Raw/sample_names)

## make project directory + make directory for ref genome
OUTDIR="/scratch/crs12448/MEVE"

# This is based off of the GATK "Best Practices" for RNAseq short variant discovery (SNPs + Indels). Found here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-

# First, mark duplicate reads from the bam files.

ml Java/17.0.6
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17

# Run MarkDuplicates to mark duplicate reads from RNAseq data. For each file in the DD, perform the mark duplicates function. Calling on JAVA with 64GB and 8 cores

#mkdir $OUTDIR/GATK/MarkDuplicates

#  gatk MarkDuplicates \
#      I=$OUTDIR/Alignment/HISAT2/BAM/${sample}.bam \
#      O=$OUTDIR/GATK/MarkDuplicates/${sample}_mark_dup.bam \
#      M=$OUTDIR/GATK/MarkDuplicates/${sample}_mark_dup_metrics.txt}

# # SplitNCigar reads to split reads at "N"s (new read to the left and to the right of the N). Not sure exactly why this is necessary....
#  mkdir $OUTDIR/GATK/SplitNCigarReads

#  gatk SplitNCigarReads \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    -I $OUTDIR/GATK/MarkDuplicates/${sample}_mark_dup.bam  \
#    -O $OUTDIR/GATK/SplitNCigarReads/${sample}_cigar.bam

# # Recalibration of base quality - This REQUIRES known variable sites. Basically, because quality of base calls is so important in the variantn calling process, and even at 99.9% quality, translates to a lot of error, this function uses machine learning to determine if there are 
# patterns of base calls that might suggest errors. Since I do not have a known-variant file, need to first call variants on the processed BAM files, then filter and use that VCF file as known sites for which this can run. The GATK toolkit suggests doing this iteratively until results are 
# robust. Not quite sure yet what they will look like but will play around with it. 

#mkdir $OUTDIR/GATK/Recalibration

# gatk BaseRecalibrator \
#    -I $OUTDIR/GATK/SplitNCigarReads/${sample}_cigar.bam \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    -O $OUTDIR/GATK/Recalibration/recal_data.table

# gatk ApplyBQSR \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    -I $OUTDIR/GATK/SplitNCigarReads/${sample}_cigar.bam \
#    --bqsr-recal-file $OUTDIR/GATK/Recalibration/recal_data.table \
#    -O $OUTDIR/GATK/Recalibration/${sample}_recal.bam
 
# gatk AnalyzeCovariates \
#   -bqsr $OUTDIR/GATK/Recalibration/recal_data.table \
#   -plots AnalyzeCovariates.pdf

# # Variant calling 

# mkdir $OUTDIR/GATK/BamFix

# gatk --java-options "-Xmx4g" AddOrReplaceReadGroups \
#         I=$OUTDIR/GATK/SplitNCigarReads/${sample}_cigar.bam \
#         O=$OUTDIR/GATK/BamFix/${sample}_cigar_fix.bam \
#         SORT_ORDER=coordinate  RGLB=seq  RGPU=1 RGPL=illumina  RGSM=${sample}.bam  \
#         CREATE_INDEX=True

# #Make directory for vcf files
# mkdir mkdir $OUTDIR/GATK/HaplotypeCaller

# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#   -I $OUTDIR/GATK/BamFix/${sample}_cigar_fix.bam \
#   -O $OUTDIR/GATK/HaplotypeCaller/${sample}.g.vcf.gz \
#   -ERC GVCF

# Filter SNPs (hard filtering to reduce the size of the files). Otherwise, the combined file will be HUGE.
# I don't think this should be bad because I am hard filtering out any variants that have poor read depth
# or quality.

#Change directory to where files are
#cd $OUTDIR/GATK/HaplotypeCaller/

#Create directory for output
#mkdir $OUTDIR/GATK/HaplotypeCaller/Filter_1

#Load modules
#ml VCFtools/0.1.16-GCC-11.2.0

#Filter each sample so the minimum read depth is 20 and quality score is 30. This should reduce the total
# number of sites to a more manageable level before genotyping
#vcftools --gzvcf ${sample}.g.vcf.gz --minQ 30 --min-meanDP 20 --recode --stdout > $OUTDIR/GATK/HaplotypeCaller/Filter_1/${sample}_filtered.g.vcf

# Now that we have done a first pass filter, combine the GVCF files to perform genotyping across all individuals. Here, I am removing the two libraries with limited reads (S256 and S266). The _2 are the new ones with high coverage.

#Combine GVCFs
cd $OUTDIR/GATK/HaplotypeCaller/Filter_1
mkdir $OUTDIR/GATK/GenotypeGVCFs

gatk CombineGVCFs \
   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
 --variant S231_filtered.g.vcf \
 --variant S242_filtered.g.vcf \
 --variant S246_filtered.g.vcf \
 --variant S247_filtered.g.vcf \
 --variant S252_filtered.g.vcf \
 --variant S256_2_filtered.g.vcf \
 --variant S263_filtered.g.vcf \
 --variant S266_2_filtered.g.vcf \
 --variant S280_filtered.g.vcf \
 --variant S295_filtered.g.vcf \
 --variant S302_filtered.g.vcf \
 --variant S303_filtered.g.vcf \
 --variant S314_filtered.g.vcf \
 --variant S316_filtered.g.vcf \
 --variant S317_filtered.g.vcf \
 --variant S319_filtered.g.vcf \
 --variant S328_filtered.g.vcf \
 --variant S336_filtered.g.vcf \
 --variant S337_filtered.g.vcf \
 --variant S338_filtered.g.vcf \
 --variant S344_filtered.g.vcf \
 --variant S345_filtered.g.vcf \
 --variant S348_filtered.g.vcf \
 --variant S350_filtered.g.vcf \
 --variant S353_filtered.g.vcf \
 --variant S357_filtered.g.vcf \
 --variant S359_filtered.g.vcf \
 --variant S367_filtered.g.vcf \
 --variant S376_filtered.g.vcf \
 --variant S380_filtered.g.vcf \
 --variant S384_filtered.g.vcf \
 --variant S388_filtered.g.vcf \
 --variant S391_filtered.g.vcf \
 --variant S392_filtered.g.vcf \
 --variant S393_filtered.g.vcf \
 --variant S406_filtered.g.vcf \
 --variant S407_filtered.g.vcf \
 --variant S408_filtered.g.vcf \
 --variant S416_filtered.g.vcf \
 --variant S420_filtered.g.vcf \
 --variant S421_filtered.g.vcf \
 --variant S422_filtered.g.vcf \
 --variant S425_filtered.g.vcf \
 --variant S426_filtered.g.vcf \
 --variant S427_filtered.g.vcf \
 --variant S432_filtered.g.vcf \
 --variant S433_filtered.g.vcf \
 --variant S435_filtered.g.vcf \
 -O $OUTDIR/GATK/GenotypeGVCFs/all_samples.g.vcf

 # Now we have a single VCF with all samples. We need to genotype them all together now, which can be done using GenotypeGVCFs as below
cd $OUTDIR/GATK/GenotypeGVCFs

 gatk GenotypeGVCFs \
   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
   -V all_samples.g.vcf \
   -O MEVE_variants_1.vcf