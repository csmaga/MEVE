#!/bin/bash
#SBATCH --job-name=GATK_SNP                 # Job name
#SBATCH --partition=highmem_p	                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=4                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/GATK_haplo.o    # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=/scratch/crs12448/MEVE/Logs/GATK_haplo.e
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
cd /scratch/crs12448/MEVE/GATK/FixBam

for i in *.bam;
do
 gatk --java-options "-Xmx120g -XX:+UseParallelGC -XX:ParallelGCThreads=4" HaplotypeCaller  \
   -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
   -I $i \
   -O $OD_4/${i/_cigar.bam/.g.vcf.gz} \
   -ERC GVCF
done



