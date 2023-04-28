#!/bin/bash
#SBATCH --job-name=GATK_SNP                 # Job name
#SBATCH --partition=highmem_p                         # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=1                       # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=200gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/GATK2_combine.o    # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=/scratch/crs12448/MEVE/Logs/GATK2_combine.e
#SBATCH --mail-user=christopher.smaga@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

# Directorty where aligned and sorted BAM files are located
#DD="/scratch/crs12448/MEVE/Alignment/HISAT2/BAM"

# Set ouput directory
#OD="/scratch/crs12448/MEVE/GATK/MarkDuplicates2"

#Load the Genome Anlysis Toolkit
ml  GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8

# Run MarkDuplicates to mark duplicate reads from RNAseq data. For each file in the DD, perform the mark duplicates function. 
# cd $DD
# for i in *.bam;
#  do
#  gatk --java-options "-Xmx120G" MarkDuplicates \
#        I=$i \
#        O=$OD/${i/.bam/_mark_dup.bam} \
#        M=$OD/${i/.bam/_mark_dup_metrics.txt}
#  done

############################################################################################################################################3

#This takes a long time, so I split into mutiple scripts. Calling on them here....

# sbatch ~/MEVE/SplitCigar_1.sh
# sbatch ~/MEVE/SplitCigar_2.sh
# sbatch ~/MEVE/SplitCigar_3.sh
# sbatch ~/MEVE/SplitCigar_4.sh

# cd $OD
#OD_2="/scratch/crs12448/MEVE/GATK/SplitNCigarReads2"

# # # Run SplitNCigarReads to split reads that span introns into separate reads
#  for i in *.bam;
#  do
#   gatk --java-options "-Xmx120G -XX:+UseParallelGC -XX:ParallelGCThreads=10" SplitNCigarReads \
#        -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#        -I $i \
#        -O $OD_2/${i/_mark_dup.bam/_cigar.bam}
#  done

############################################################################################################################################
# Call SNPs for each sample individually using haplotype caller in GATK. With the -ERC GVCF option, temporary .gvcf files are created that can then be merged into one, and the genotype calls can be made across all samples.
# From my understanding, this will allow genotypes to be called in samples even if they match the reference, as long as one sample contains a SNP in that position. 

# cd $OD_2

# OD_3="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/GVCF"

# for i in *.bam;
# do
#  gatk --java-options "-Xmx120g -XX:+UseParallelGC -XX:ParallelGCThreads=8" HaplotypeCaller  \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#     -I $i \
#     -O $OD_3/${i/_cigar.bam/.g.vcf.gz} \
#     -ERC GVCF
#  done
 # !!!!!!!!!!!!!!!!!  This gives me an error that I am trying to run on more than one sample. From a little googling, that may be because my BAM files are missing the @RG identifier. Sure enough,
 # when I look at the BAM files, there is no @RG. Luckily, GATK has a tool to add RG to samples with the below code. I don't think the information in RG is useful for my application,  so i only added
 # the required fields, which are RGLB (read group library), RGPU (read group platform unit?), RGPL (read group platform), and RGSM (read group sample - this is the important one). For all samples, I will leave all 
 # identifiers the same except for RG sample name. 

# Fix the read groups......
# cd $OD_2

# for i in *.bam;
#  do
#  gatk --java-options "-Xmx60g" AddOrReplaceReadGroups \
#          I=$i  \
#          O=/scratch/crs12448/MEVE/GATK/FixBam2/${i/_cigar.bam/_cigar_fix.bam}  \
#          SORT_ORDER=coordinate  RGLB=seq  RGPU=1 RGPL=illumina  RGSM=${i/_cigar.bam/}  \
#          CREATE_INDEX=True
#  done



############################################################################################################################################
# Retry haplotype caller.....
#Call SNPs for each sample individually using haplotype caller in GATK. With the -ERC GVCF option, temporary .gvcf files are created that can then be merged into one, and the genotype calls can be made across all samples.
# From my understanding, this will allow genotypes to be called in samples even if they match the reference, as long as one sample contains a SNP in that position. 
# !!!!!!!!!!!!!!!!!!!!!!!!!!
####### I don't think multithreading works for this, which may be why I got such weird results from my first run through the entire pipeline. Trying to run everything again but without multiple threads to see if it fixes the issues.....HOPEFULLY!

#OD_4="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/GVCF"
# cd /scratch/crs12448/MEVE/GATK/FixBam

# for i in *.bam;
# do
#  gatk --java-options "-Xmx120g" HaplotypeCaller  \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    -I $i \
#    -O $OD_4/${i/_cigar.bam/.g.vcf.gz} \
#    -ERC GVCF
# done

# ## AGAIN this takes a LONG time, so split into several scripts
# cd ~/MEVE

# sbatch HaploCaller_1.sh
# sbatch HaploCaller_2.sh
# sbatch HaploCaller_3.sh
# sbatch HaploCaller_4.sh


# The GVCF files are HUGE because they count every base covered regardless of whether it is variable or not (so that genotyping cann be done with multiple individuals). Combining the raw files would take forever and lead to a massive GVCF, so first filter out low quality and low depth SNPs
# Perform the filtering with vcftools

ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

# FILTER_OD="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/FilterGVCF2"

# cd /scratch/crs12448/MEVE/GATK/HaplotypeCaller/GVCF2
# # # This code only keeps sites that have a mean depth > 20, minimum quality of 30. No max depth set because in RNAseq, I expect some depths are really high biologically
#  for i in S231 S242 S246 S247 S252 S256_2 S263 S266_2 S280 S295 S302 S316 S317 S319 S337 S344 S359 S376 S388 S391 S392 S393 S406 S432;
# do
#   vcftools --vcf $i --minQ 30 --min-meanDP 20 --recode --recode-INFO-all --stdout > $FILTER_OD/${i}_filtered.vcf
# done

############################################################################################################################################

                                        # DO NOT NEED TO RUN THIS BLOCK OF CODE #

############################################################################################################################################

# Retry haplotype caller, this time on all samples at once - maybe this will change the high heterozygosity values I saw with the previous pipeline


# OD_5="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/Combined_samples"
# cd /scratch/crs12448/MEVE/GATK/FixBam

# gatk --java-options "-Xmx300g -XX:+UseParallelGC -XX:ParallelGCThreads=12" HaplotypeCaller  \
#  -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#  -I S231_cigar_fix.bam S242_cigar_fix.bam S246_cigar_fix.bam S247_cigar_fix.bam S252_cigar_fix.bam S256_2_cigar_fix.bam S263_cigar_fix.bam S266_2_cigar_fix.bam S280_cigar_fix.bam S295_cigar_fix.bam S302_cigar_fix.bam S316_cigar_fix.bam \
#     S317_cigar_fix.bam S319_cigar_fix.bam S337_cigar_fix.bam S344_cigar_fix.bam S359_cigar_fix.bam S376_cigar_fix.bam S388_cigar_fix.bam S391_cigar_fix.bam S392_cigar_fix.bam S393_cigar_fix.bam S406_cigar_fix.bam S432_cigar_fix.bam \
#  --native-pair-hmm-threads 12 -O $OD_5/AP_WO_combined.vcf 




#########################################################################
# ## AGAIN this takes a LONG time, so split into several scripts
# cd ~/MEVE

# sbatch HaploCaller_1.sh
# sbatch HaploCaller_2.sh
# sbatch HaploCaller_3.sh
# sbatch HaploCaller_4.sh


# The GVCF files are HUGE because they count every base covered regardless of whether it is variable or not (so that genotyping cann be done with multiple individuals). Combining the raw files would take forever and lead to a massive GVCF, so first filter out low quality and low depth SNPs
# Perform the filtering with vcftools
# ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

#FILTER_OD="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/Filter_GVCF"

# cd /scratch/crs12448/MEVE/GATK/HaplotypeCaller/GVCF
# # This code only keeps sites that have a mean depth > 20, minimum quality of 30. No max depth set because in RNAseq, I expect some depths are really high biologically
#  for i in S231 S242 S246 S247 S252 S256_2 S263 S266_2 S280 S295 S302 S316 S317 S319 S337 S344 S359 S376 S388 S391 S392 S393 S406 S432;
#  do
#   vcftools --vcf $i --minQ 30 --min-meanDP 20 --recode --stdout > $FILTER_OD/${i}_filtered.g.vcf
# done

##############################################################################################################################################################################################################################################################################

# Combine the GVCF files into one using the CombineGVCFs function. There is also a GenomicDBImport option that is better for handling more files, but I understand it less and think this might work for only 24 samples (plus they are filtered and have less SNPs):
 
#  FILTER_OD="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/FilterGVCF2"
#  cd $FILTER_OD
 
#  gatk --java-options "-Xmx200g" CombineGVCFs \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    --variant S231_filtered.vcf \
#    --variant S242_filtered.vcf \
#    --variant S246_filtered.vcf \
#    --variant S247_filtered.vcf \
#    --variant S252_filtered.vcf \
#    --variant S256_2_filtered.vcf \
#    --variant S263_filtered.vcf \
#    --variant S266_2_filtered.vcf \
#    --variant S280_filtered.vcf \
#    --variant S295_filtered.vcf \
#    --variant S302_filtered.vcf \
#    --variant S316_filtered.vcf \
#    --variant S317_filtered.vcf \
#    --variant S319_filtered.vcf \
#    --variant S337_filtered.vcf \
#    --variant S344_filtered.vcf \
#    --variant S359_filtered.vcf \
#    --variant S376_filtered.vcf \
#    --variant S388_filtered.vcf \
#    --variant S391_filtered.vcf \
#    --variant S392_filtered.vcf \
#    --variant S393_filtered.vcf \
#    --variant S406_filtered.vcf \
#    --variant S432_filtered.vcf \
#    -O all_samples.vcf

######################################################
# Combine the GVCF files into one using the CombineGVCFs function. There is also a GenomicDBImport option that is better for handling more files, but I understand it less and think this might work for only 24 samples (plus they are filtered and have less SNPs):
 
 FILTER_OD="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/GVCF2"
 cd $FILTER_OD
 
 gatk --java-options "-Xmx200g" CombineGVCFs \
   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
   --variant S231 \
   --variant S242 \
   --variant S246 \
   --variant S247 \
   --variant S252 \
   --variant S256_2 \
   --variant S263 \
   --variant S266_2 \
   --variant S280 \
   --variant S295 \
   --variant S302 \
   --variant S316 \
   --variant S317 \
   --variant S319 \
   --variant S337 \
   --variant S344 \
   --variant S359 \
   --variant S376 \
   --variant S388 \
   --variant S391 \
   --variant S392 \
   --variant S393 \
   --variant S406 \
   --variant S432 \
   -O /scratch/crs12448/MEVE/GATK/HaplotypeCaller/Combined_samples/all_samples.vcf


# Now we have a single VCF with all samples. We need to genotype them all together now, which can be done using GenotypeGVCFs as below
#cd $FILTER_OD

 #gatk --java-options "-Xmx200g" GenotypeGVCFs \
  # -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
   #-V all_samples.g.vcf \
   #-O /scratch/crs12448/MEVE/GATK/GenotypeGVCF/new/AP_WO_SNPs_new.vcf

###################################################

# Filter SNPs
#vcftools --vcf AP_WO_SNPs.vcf --remove-indels --maf 0.1 --max-missing 0.8 --minQ 30 --min-meanDP 20 --minDP 20 --recode --stdout > AP_WO_SNPs_filter.vcf


# Create a matrix of SNPs for exploratory analysis
#bcftools query AP_WO_SNPs_filter.vcf -f '%CHROM\t%POS[\t%GT]\n' > snp_matrix.txt

# Extract sample names
#bcftools query -l AP_WO_SNPs_filter.vcf > sample_names.txt

# Only keep loci with 2 alleles - went from 4871 to 4751 sites
#vcftools --vcf AP_WO_SNPs_filter.vcf --min-alleles 2 --max-alleles 2 --recode --out AP_WO_SNPs_filter_maxtwo.vcf

#vcftools --vcf AP_WO_SNPs_filter_maxtwo.vcf 


#  Transfer those files to local machine for use in R

#vcftools --vcf AP_WO_SNPs.vcf --remove_indels --recode --out AP_WO_SNPs_only.vcf