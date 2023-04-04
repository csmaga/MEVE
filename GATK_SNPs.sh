#!/bin/bash
#SBATCH --job-name=GATK_SNP                 # Job name
#SBATCH --partition=highmem_p                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=8                       # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=200gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/GATK_combine_gvcf.o    # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=/scratch/crs12448/MEVE/Logs/GATK_combine_gvcf.e
#SBATCH --mail-user=christopher.smaga@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

# Directorty where aligned and sorted BAM files are located
DD="/scratch/crs12448/MEVE/Alignment/HISAT2/BAM"

# Set ouput directory
OD="/scratch/crs12448/MEVE/GATK/MarkDuplicates"

#Load the Genome Anlysis Toolkit
ml  GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8

# Run MarkDuplicates to mark duplicate reads from RNAseq data. For each file in the DD, perform the mark duplicates function. Calling on JAVA with 64GB and 8 cores
#cd $DD
# for i in *.bam;
# do
# gatk --java-options "-Xmx64G -XX:+UseParallelGC -XX:ParallelGCThreads=8" MarkDuplicates \
#       I=$i \
#       O=$OD/${i/.bam/_mark_dup.bam} \
#       M=$OD/${i/.bam/_mark_dup_metrics.txt}
# done

############################################################################################################################################3

#This takes a long time, so I split into mutiple scripts. Calling on them here....

# sbatch ~/MEVE/SplitCigar_1.sh
# sbatch ~/MEVE/SplitCigar_2.sh
# sbatch ~/MEVE/SplitCigar_3.sh
# sbatch ~/MEVE/SplitCigar_4.sh

# cd $OD
OD_2="/scratch/crs12448/MEVE/GATK/SplitNCigarReads"

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

#cd $OD_2

#OD_3="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/GVCF"

# for i in *.bam;
# do
#  gatk --java-options "-Xmx120g -XX:+UseParallelGC -XX:ParallelGCThreads=8" HaplotypeCaller  \
#    -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
#    -I $i \
#    -O $OD_3/${i/_cigar.bam/.g.vcf.gz} \
#    -ERC GVCF
# done
 # !!!!!!!!!!!!!!!!!  This gives me an error that I am trying to run on more than one sample. From a little googling, that may be because my BAM files are missing the @RG identifier. Sure enough,
 # when I look at the BAM files, there is no @RG. Luckily, GATK has a tool to add RG to samples with the below code. I don't think the information in RG is useful for my application,  so i only added
 # the required fields, which are RGLB (read group library), RGPU (read group platform unit?), RGPL (read group platform), and RGSM (read group sample - this is the important one). For all samples, I will leave all 
 # identifiers the same except for RG sample name. 

# Fix the read groups......
#cd $OD_2

# for i in *.bam;
# do
# gatk --java-options "-Xmx64g -XX:+UseParallelGC -XX:ParallelGCThreads=4" AddOrReplaceReadGroups \
#         I=$i  O=/scratch/crs12448/MEVE/GATK/FixBam/${i/_cigar.bam/_cigar_fix.bam}  \
#         SORT_ORDER=coordinate  RGLB=seq  RGPU=1 RGPL=illumina  RGSM=${i/_cigar.bam/}  \
#         CREATE_INDEX=True
# done



############################################################################################################################################
# Retry haplotype caller.....
#Call SNPs for each sample individually using haplotype caller in GATK. With the -ERC GVCF option, temporary .gvcf files are created that can then be merged into one, and the genotype calls can be made across all samples.
# From my understanding, this will allow genotypes to be called in samples even if they match the reference, as long as one sample contains a SNP in that position. 


OD_4="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/GVCF"
# cd /scratch/crs12448/MEVE/GATK/FixBam

# for i in *.bam;
# do
#  gatk --java-options "-Xmx120g -XX:+UseParallelGC -XX:ParallelGCThreads=4" HaplotypeCaller  \
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
# ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0

# FILTER_OD="/scratch/crs12448/MEVE/GATK/HaplotypeCaller/Filter_GVCF"

# cd /scratch/crs12448/MEVE/GATK/HaplotypeCaller/GVCF
# # This code only keeps sites that have a mean depth > 20, minimum quality of 30. No max depth set because in RNAseq, I expect some depths are really high biologically
#  for i in S231 S242 S246 S247 S252 S256_2 S263 S266_2 S280 S295 S302 S316 S317 S319 S337 S344 S359 S376 S388 S391 S392 S393 S406 S432;
#  do
#   vcftools --vcf $i --minQ 30 --min-meanDP 20 --recode --stdout > $FILTER_OD/${i}_filtered.g.vcf
# done




# Combine the GVCF files into one using the CombineGVCFs function. There is also a GenomicDBImport option that is better for handling more files, but I understand it less and think this might work for only 24 samples (plus they are filtered and have less SNPs):
 
 cd $FILTER_OD
 
 gatk --java-options "-Xmx200g -XX:+UseParallelGC -XX:ParallelGCThreads=8" CombineGVCFs \
   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
   --variant S231.vcf \
   --variant S242.vcf \
   --variant S246.vcf \
   --variant S247.vcf \
   --variant S252.vcf \
   --variant S256_2.vcf \
   --variant S263.vcf \
   --variant S266_2.vcf \
   --variant S280.vcf \
   --variant S295.vcf \
   --variant S302.vcf \
   --variant S316.vcf \
   --variant S317.vcf \
   --variant S319.vcf \
   --variant S337.vcf \
   --variant S344.vcf \
   --variant S359.vcf \
   --variant S376.vcf \
   --variant S388.vcf \
   --variant S391.vcf \
   --variant S392.vcf \
   --variant S393.vcf \
   --variant S406.vcf \
   --variant S432.vcf \
   -O all_samples.vcf














# OD="/scratch/crs12448/MEVE/GATK/Merge"


#    gatk --java-options "-Xmx250g -Xms200g -XX:+UseParallelGC -XX:ParallelGCThreads=4" GenomicsDBImport \
#       -V $OD_4/S392  \
#       -V $OD_4/S393  \
#       --genomicsdb-workspace-path /scratch/crs12448/MEVE/GATK/Merge/GenomicDB \
#       -L $OD_4/S392 \
#       -L $OD_4/S393


