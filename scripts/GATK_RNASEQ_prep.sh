#!/bin/bash
#SBATCH --job-name=MEVE_GATK
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=72:00:00
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j
#SBATCH --mail-user=crs12448@uga.edu
#SBATCH --mail-type=END,FAIL
###SBATCH --array=1-50


#sample=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/crs12448/MEVE/Data/Raw/sample_names)

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

# #Combine GVCFs
# cd $OUTDIR/GATK/HaplotypeCaller/Filter_1
# mkdir $OUTDIR/GATK/GenotypeGVCFs

# gatk CombineGVCFs \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#  --variant S231_filtered.g.vcf \
#  --variant S242_filtered.g.vcf \
#  --variant S246_filtered.g.vcf \
#  --variant S247_filtered.g.vcf \
#  --variant S252_filtered.g.vcf \
#  --variant S256_2_filtered.g.vcf \
#  --variant S263_filtered.g.vcf \
#  --variant S266_2_filtered.g.vcf \
#  --variant S280_filtered.g.vcf \
#  --variant S295_filtered.g.vcf \
#  --variant S302_filtered.g.vcf \
#  --variant S303_filtered.g.vcf \
#  --variant S314_filtered.g.vcf \
#  --variant S316_filtered.g.vcf \
#  --variant S317_filtered.g.vcf \
#  --variant S319_filtered.g.vcf \
#  --variant S328_filtered.g.vcf \
#  --variant S336_filtered.g.vcf \
#  --variant S337_filtered.g.vcf \
#  --variant S338_filtered.g.vcf \
#  --variant S344_filtered.g.vcf \
#  --variant S345_filtered.g.vcf \
#  --variant S348_filtered.g.vcf \
#  --variant S350_filtered.g.vcf \
#  --variant S353_filtered.g.vcf \
#  --variant S357_filtered.g.vcf \
#  --variant S359_filtered.g.vcf \
#  --variant S367_filtered.g.vcf \
#  --variant S376_filtered.g.vcf \
#  --variant S380_filtered.g.vcf \
#  --variant S384_filtered.g.vcf \
#  --variant S388_filtered.g.vcf \
#  --variant S391_filtered.g.vcf \
#  --variant S392_filtered.g.vcf \
#  --variant S393_filtered.g.vcf \
#  --variant S406_filtered.g.vcf \
#  --variant S407_filtered.g.vcf \
#  --variant S408_filtered.g.vcf \
#  --variant S416_filtered.g.vcf \
#  --variant S420_filtered.g.vcf \
#  --variant S421_filtered.g.vcf \
#  --variant S422_filtered.g.vcf \
#  --variant S425_filtered.g.vcf \
#  --variant S426_filtered.g.vcf \
#  --variant S427_filtered.g.vcf \
#  --variant S432_filtered.g.vcf \
#  --variant S433_filtered.g.vcf \
#  --variant S435_filtered.g.vcf \
#  -O $OUTDIR/GATK/GenotypeGVCFs/all_samples.g.vcf


#Change directory to where files are
#cd $OUTDIR/GATK/HaplotypeCaller/


# Combine the GVCF files to perform genotyping across all individuals. I am NOT filtering beforehand
# # Here, I am removing the two libraries with limited reads (S256 and S266). The _2 are the new ones with high coverage.
# gatk CombineGVCFs \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#  --variant S231.g.vcf.gz \
#  --variant S242.g.vcf.gz \
#  --variant S246.g.vcf.gz \
#  --variant S247.g.vcf.gz \
#  --variant S252.g.vcf.gz \
#  --variant S256_2.g.vcf.gz \
#  --variant S263.g.vcf.gz \
#  --variant S266_2.g.vcf.gz \
#  --variant S280.g.vcf.gz \
#  --variant S295.g.vcf.gz \
#  --variant S302.g.vcf.gz \
#  --variant S303.g.vcf.gz \
#  --variant S314.g.vcf.gz \
#  --variant S316.g.vcf.gz \
#  --variant S317.g.vcf.gz \
#  --variant S319.g.vcf.gz \
#  --variant S328.g.vcf.gz \
#  --variant S336.g.vcf.gz \
#  --variant S337.g.vcf.gz \
#  --variant S338.g.vcf.gz \
#  --variant S344.g.vcf.gz \
#  --variant S345.g.vcf.gz \
#  --variant S348.g.vcf.gz \
#  --variant S350.g.vcf.gz \
#  --variant S353.g.vcf.gz \
#  --variant S357.g.vcf.gz \
#  --variant S359.g.vcf.gz \
#  --variant S367.g.vcf.gz \
#  --variant S376.g.vcf.gz \
#  --variant S380.g.vcf.gz \
#  --variant S384.g.vcf.gz \
#  --variant S388.g.vcf.gz \
#  --variant S391.g.vcf.gz \
#  --variant S392.g.vcf.gz \
#  --variant S393.g.vcf.gz \
#  --variant S406.g.vcf.gz \
#  --variant S407.g.vcf.gz \
#  --variant S408.g.vcf.gz \
#  --variant S416.g.vcf.gz \
#  --variant S420.g.vcf.gz \
#  --variant S421.g.vcf.gz \
#  --variant S422.g.vcf.gz \
#  --variant S425.g.vcf.gz \
#  --variant S426.g.vcf.gz \
#  --variant S427.g.vcf.gz \
#  --variant S432.g.vcf.gz \
#  --variant S433.g.vcf.gz \
#  --variant S435.g.vcf.gz \
#  -O $OUTDIR/GATK/CombineGVCFs/all_samples_unfiltered.g.vcf.gz
 # Now we have a single VCF with all samples. We need to genotype them all together now, which can be done using GenotypeGVCFs as below
cd $OUTDIR/GATK/GenotypeGVCFs

  gatk GenotypeGVCFs --java-options "-Xmx32g"\
   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
   -V $OUTDIR/GATK/CombineGVCFs/all_samples_unfiltered.g.vcf.gz \
   -L $OUTDIR/PopGen/gene_bed_sorted.bed -ip 100 \
   -all-sites \
   -O MEVE_variants_unfiltered_allgenes.vcf

# Filter SNPs


#vcftools --gvcf MEVE_variants_unfiltered.g.vcf --remove-indels --maf 0.1 --max-missing 0.9 --minQ 30 --min-meanDP 10 --minDP 10 --recode --stdout > $OUTDIR/Filtered/MEVE_SNPs_filter_1.vcf
#vcftools --gvcf MEVE_variants_unfiltered.g.vcf --remove-indels --maf 0.05 --max-missing 0.9 --minQ 30 --min-meanDP 10 --minDP 10 --recode --stdout > Filtered/MEVE_SNPs_filter_2.vcf
#vcftools --gvcf MEVE_variants_unfiltered.g.vcf --remove-indels --maf 0.1 --max-missing 0.95 --minQ 30 --min-meanDP 10 --minDP 10 --recode --stdout > $OUTDIR/Filtered/MEVE_SNPs_filter_3.vcf
#vcftools --gvcf MEVE_variants_unfiltered.g.vcf --remove-indels --maf 0.05 --max-missing 0.95 --minQ 30 --min-meanDP 10 --minDP 10 --recode --stdout > $OUTDIR/Filtered/MEVE_SNPs_filter_4.vcf
#vcftools --gvcf MEVE_variants_unfiltered.g.vcf --remove-indels --maf 0.1 --max-missing 0.85 --minQ 30 --min-meanDP 10 --minDP 10 --recode --stdout > $OUTDIR/Filtered/MEVE_SNPs_filter_5.vcf
#vcftools --gvcf MEVE_variants_unfiltered.g.vcf --remove-indels --maf 0.05 --max-missing 0.85 --minQ 30 --min-meanDP 10 --minDP 10 --recode --stdout > $OUTDIR/Filtered/MEVE_SNPs_filter_6.vcf


# First, filter by depth - only retain sites (SNPs) that have a quality above 30 (all of them in this case) and min mean depth of 10.
#vcftools --gzvcf MEVE_SNPs_filter_depth.vcf --remove-indels --min-meanDP 10 --recode --stdout > Filtered/MEVE_SNPs_filter_depth_sdpth.vcf # 180443 sites remaining


#vcftools --gvcf MEVE_variants_unfiltered.vcf.gz --remove-indels --minQ 30 --minDP 10 --recode --stdout > Filtered/MEVE_SNPs_filter_depth.vcf
#bcftools query -f "A\n" MEVE_SNPs_filter_depth.vcf | wc -l

# 1561991 sites retained

# Now, filter by missing data. Only retain sites present in at least 95%, 90%, 85% of individuals 




#vcftools --gvcf MEVE_SNPs_filter_depth.vcf --remove-indels --max-missing 0.95 --recode --stdout > Filtered/MEVE_SNPs_filter_depth_m95.vcf # 48972 sites remaining
#vcftools --gvcf MEVE_variants_unfiltered.vcf.gz --remove-indels --max-missing 0.90 --recode --stdout > Filtered/MEVE_SNPs_filter_depth_m90.vcf
#vcftools --gvcf MEVE_variants_unfiltered.vcf.gz --remove-indels --max-missing 0.85 --recode --stdout > Filtered/MEVE_SNPs_filter_depth_m85.vcf

# gatk IndexFeatureFile \
#     -I MEVE_variants_unfiltered.vcf.gz


# gatk VariantFiltration -V MEVE_variants_unfiltered.vcf.gz  --filter-expression "QD < 2.0" --filter-name "QD2" \
#  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
#  --filter-expression "FS > 60.0" --filter-name "FS60" \
#  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
#  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#  -O Filtered/MEVE_SNPs.filtered.vcf.gz


#ml  VCFtools/0.1.16-GCC-11.2.0
#cd /scratch/crs12448/MEVE/GATK/GenotypeGVCFs/Filtered

# This filtering further requires that sites are SNPs only, marks genotypes as missing if any individual genotype has a depth of less than 10, and requires that sites being present in at least 90% of inviduals. 
#vcftools --gzvcf MEVE_SNPs.filtered.vcf.gz --remove-indels --minDP 10 --maf 0.05  --max-missing 0.90 --recode --stdout > /scratch/crs12448/MEVE/GATK/GenotypeGVCFs/Filtered/MEVE_SNPs_filtered_013024.vcf

# This gives me 37,434 SNPs. I think this is best approach to take because it filteres for recommended strand biases
# and quality by depth instead of just quality as recommended by GATK Best Practiices for Hard Filtering Variants. I imagine thes
# are the same SNPs as the previous analysis, just slightly fewer in number. 

# gatk VariantFiltration -V MEVE_variants_unfiltered_allsites.vcf  --filter-expression "QD < 2.0" --filter-name "QD2" \
#  --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
#  --filter-expression "FS > 60.0" --filter-name "FS60" \
#  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
#  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#  -O Filtered/MEVE_variants_filtered_allsites.vcf.gz



# # set filters
# MAF=0.05
# MISS=0.9
# QUAL=30
# MIN_DEPTH=10

# cd $OUTDIR/GATK/GenotypeGVCFs/
# # perform the filtering with vcftools
# vcftools --gzvcf MEVE_variants_unfiltered.vcf.gz \
# --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
# --min-meanDP $MIN_DEPTH  \
# --minDP $MIN_DEPTH --recode --stdout > $OUTDIR/GATK/GenotypeGVCFs/Filtered/MEVE_SNPs_filtered_011124.vcf

## This filtering approach gives me 37,434 SNPs















# Trying basequality recalibration to see if it makes any difference
#
#cd $OUTDIR/GATK/Recalibration

#gatk BaseRecalibrator \
#    -I $OUTDIR/GATK/BamFix/${sample}_cigar_fix.bam \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    --known-sites $OUTDIR/GATK/GenotypeGVCFs/Filtered/MEVE_SNPs_filtered_011124.vcf \
#    -O $OUTDIR/GATK/Recalibration/${sample}_recal_data.table

 #gatk ApplyBQSR \
 #   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
 #   -I $OUTDIR/GATK/SplitNCigarReads/${sample}_cigar.bam \
 #   --$OUTDIR/GATK/Recalibration/recal_data.table \
 #   -O $OUTDIR/GATK/Recalibration/${sample}_recal.bam
 
#gatk AnalyzeCovariates \
#   -bqsr $OUTDIR/GATK/Recalibration/S266_recal_data.table \
 #  -plots S266_AnalyzeCovariates.pdf


## Apparently, the read groups need to be added to the new BAM files after correctig them (cannot use HaplotypeCaller on the recal BAM files as is). 
# So, need to re-run AddOrReplaceReadGroups

# gatk --java-options "-Xmx4g" AddOrReplaceReadGroups \
#         I=$OUTDIR/GATK/Recalibration/${sample}_recal.bam \
#         O=$OUTDIR/GATK/BamFix2/${sample}_cigar_fix.bam \
#         SORT_ORDER=coordinate  RGLB=seq  RGPU=1 RGPL=illumina  RGSM=${sample}.bam  \
#         CREATE_INDEX=True

# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#   -I $OUTDIR/GATK/BamFix2/${sample}_cigar_fix.bam \
#   -O $OUTDIR/GATK/HaplotypeCaller2/${sample}.g.vcf.gz \
#   -ERC GVCF

#cd $OUTDIR/GATK/HaplotypeCaller2

# gatk CombineGVCFs \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#  --variant S231.g.vcf.gz \
#  --variant S242.g.vcf.gz \
#  --variant S246.g.vcf.gz \
#  --variant S247.g.vcf.gz \
#  --variant S252.g.vcf.gz \
#  --variant S256_2.g.vcf.gz \
#  --variant S263.g.vcf.gz \
#  --variant S266_2.g.vcf.gz \
#  --variant S280.g.vcf.gz \
#  --variant S295.g.vcf.gz \
#  --variant S302.g.vcf.gz \
#  --variant S303.g.vcf.gz \
#  --variant S314.g.vcf.gz \
#  --variant S316.g.vcf.gz \
#  --variant S317.g.vcf.gz \
#  --variant S319.g.vcf.gz \
#  --variant S328.g.vcf.gz \
#  --variant S336.g.vcf.gz \
#  --variant S337.g.vcf.gz \
#  --variant S338.g.vcf.gz \
#  --variant S344.g.vcf.gz \
#  --variant S345.g.vcf.gz \
#  --variant S348.g.vcf.gz \
#  --variant S350.g.vcf.gz \
#  --variant S353.g.vcf.gz \
#  --variant S357.g.vcf.gz \
#  --variant S359.g.vcf.gz \
#  --variant S367.g.vcf.gz \
#  --variant S376.g.vcf.gz \
#  --variant S380.g.vcf.gz \
#  --variant S384.g.vcf.gz \
#  --variant S388.g.vcf.gz \
#  --variant S391.g.vcf.gz \
#  --variant S392.g.vcf.gz \
#  --variant S393.g.vcf.gz \
#  --variant S406.g.vcf.gz \
#  --variant S407.g.vcf.gz \
#  --variant S408.g.vcf.gz \
#  --variant S416.g.vcf.gz \
#  --variant S420.g.vcf.gz \
#  --variant S421.g.vcf.gz \
#  --variant S422.g.vcf.gz \
#  --variant S425.g.vcf.gz \
#  --variant S426.g.vcf.gz \
#  --variant S427.g.vcf.gz \
#  --variant S432.g.vcf.gz \
#  --variant S433.g.vcf.gz \
#  --variant S435.g.vcf.gz \
#  -O $OUTDIR/GATK/CombineGVCFs2/all_samples_unfiltered.g.vcf.gz

#cd $OUTDIR/GATK/GenotypeGVCFs2

  # gatk GenotypeGVCFs \
  #  -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
  #  -V $OUTDIR/GATK/CombineGVCFs2/all_samples_unfiltered.g.vcf.gz \
  #  -O MEVE_variants_unfiltered.vcf

# # perform the filtering with vcftools
#vcftools --gzvcf M_unfiltered.vcf.gz  --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH  --minDP $MIN_DEPTH --recode --stdout > $OUTDIR/GATK/GenotypeGVCFs/Filtered/MEVE_SNPs_filtered_011124.vcf
# # set filters
# MAF=0.05
# MISS=0.9
# QUAL=30
# MIN_DEPTH=10
