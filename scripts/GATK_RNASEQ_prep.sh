#!/bin/bash
#SBATCH --job-name=MEVE_GATK
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --time=24:00:00
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

#Change directory to where files are
#cd $OUTDIR/GATK/HaplotypeCaller/

#Create directory for output
#mkdir $OUTDIR/GATK/HaplotypeCaller/Filter_1

#Load modules
#ml VCFtools/0.1.16-GCC-11.2.0

# #Combine GVCFs
# cd $OUTDIR/GATK/HaplotypeCaller/Filter_1
# mkdir $OUTDIR/GATK/GenotypeGVCFs

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
#cd $OUTDIR/GATK/GenotypeGVCFs

# gatk GenotypeGVCFs --java-options "-Xmx120g" \
#     -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#     -V $OUTDIR/GATK/CombineGVCFs/all_samples_unfiltered.g.vcf.gz \
#     -L $OUTDIR/PopGen/gene_bed_sorted.bed \
#     -all-sites \
#     -O MEVE_variants_unfiltered_allgenes_02_27.vcf.gz

# Filter SNPs

# gatk IndexFeatureFile \
#     -I MEVE_variants_unfiltered.vcf.gz

# gatk VariantFiltration -V MEVE_variants_unfiltered_allgenes_02_27.vcf.gz  --filter-expression "QD < 2.0" --filter-name "QD2" \
#   --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#   --filter-expression "SOR > 3.0" --filter-name "SOR3" \
#   --filter-expression "FS > 60.0" --filter-name "FS60" \
#   --filter-expression "MQ < 40.0" --filter-name "MQ40" \
#   --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#   --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#   -O Filtered/MEVE_SNPs.filtered.allgenes_02_28.vcf.gz

#ml  VCFtools/0.1.16-GCC-11.2.0
#cd /scratch/crs12448/MEVE/GATK/GenotypeGVCFs

# This filtering further requires that sites are SNPs only, marks genotypes as missing if any individual genotype has a depth of less than 10, and requires that sites being present in at least 90% of inviduals. 
#vcftools --gzvcf MEVE_SNPs.filtered.vcf.gz --remove-indels --minDP 10 --maf 0.05  --max-missing 0.90 --recode --stdout > /scratch/crs12448/MEVE/GATK/GenotypeGVCFs/Filtered/MEVE_SNPs_filtered_013024.vcf
#vcftools --vcf MEVE_variants_filtered_allgenes_PASS_02_29.vcf --minDP 10 --maf 0.05  --max-missing 0.90 --recode --stdout > /scratch/crs12448/MEVE/GATK/GenotypeGVCFs/Filtered/MEVE_variants_filtered_PASS_filt2_02_29.vcf

# This gives me 37,434 SNPs. I think this is best approach to take because it filteres for recommended strand biases
# and quality by depth instead of just quality as recommended by GATK Best Practiices for Hard Filtering Variants. I imagine thes
# are the same SNPs as the previous analysis, just slightly fewer in number. 

##########################################################################################################################
## For the all genes file with non-varient sites too - filter the same way as above
# I am not sure for the strand bias filters whether they should be applied for invariant sites. From what 
# I understand, invariant sites would not have strand bias because there is nothing to compare between
# the reference and alternative alleles (only reads supporting the reference are present)

# Regardless, I can later select variants based on specific filters if need. For example, only selecting variants
# with QD > 2 PASS to use for downstream analyses

# gatk VariantFiltration -V MEVE_variants_unfiltered_allgenes.vcf  --filter-expression "QD < 2.0" --filter-name "QD2" \
#    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
#    --filter-expression "FS > 60.0" --filter-name "FS60" \
#    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
#    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#    -O Filtered/MEVE_variants_filtered_allgenes.vcf.gz

# gatk SelectVariants \
#       -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta  \
#       -V Filtered/MEVE_SNPs.filtered.allgenes_02_28.vcf.gz \
#       --exclude-filtered TRUE \
#       --select-type-to-include NO_VARIATION \
#       --select-type-to-include SNP \
#       -O Filtered/MEVE_variants_filtered_allgenes_PASS_02_29.vcf

#  gatk SelectVariants \
#      -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta  \
#      -V Filtered/MEVE_variants_filtered_allgenes.vcf.gz \
#      --select-type-to-include SNP \
#      --select "QUAL > 30.0" \
#      --select "QD > 2.0" \
#      --exclude-non-variants FALSE \
#      -O Filtered/MEVE_variants_filtered_allgenes_PASS2.vcf



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





# gatk GenotypeGVCFs -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta -V /scratch/crs12448/MEVE/GATK/CombineGVCFs/all_samples_unfiltered.g.vcf.gz -L /scratch/crs12448/MEVE/PopGen/gene_bed_sorted_sub.bed -all-sites -O sub_test.vcf
# # 

# gatk VariantFiltration -V sub_test.vcf  --filter-expression "QD < 2.0" --filter-name "QD2" \
#    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
#    --filter-expression "FS > 60.0" --filter-name "FS60" \
#    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
#    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#    -O sub_test_filtered.vcf.gz


#  gatk SelectVariants \
#      -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta  \
#      -V Filtered/MEVE_variants_filtered_allgenes.vcf.gz \
#      --exclude-non-variants FALSE \
#      --exclude-filtered TRUE \
#      -O sub.vcf

# vcftools --vcf MEVE_variants_filtered_allgenes_PASS_02_29.vcf \
# --max-maf 0 \
# --max-missing 0.90 \
# --minDP 10  \
# --recode --stdout > MEVE_invariant.vcf

# vcftools --vcf MEVE_variants_filtered_allgenes_PASS_02_29.vcf \
# --mac 1 \
# --maf 0.05 \
# --max-missing 0.90 \
# --minDP 10  \
# --recode --stdout > MEVE_invariant.vcf
cd /scratch/crs12448/MEVE/GATK/GenotypeGVCFs/Filtered
ml BCFtools
bcftools view -H MEVE_SNPs.filtered.allgenes_02_28.vcf.gz | wc -l > num_snps
bcftools view -H MEVE_SNPs.filtered.vcf.gz | wc -l > num_snps2
bcftools view -H MEVE_variants_filtered_allgenes.vcf.gz | wc -l > num_snps3
bcftools view -H MEVE_variants_filtered_allsites.vcf.gz| wc -l > num_snps4




