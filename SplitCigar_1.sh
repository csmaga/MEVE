#SBATCH --job-name=GATK_SNP                 # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=2                           # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=32gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/crs12448/MEVE/Logs/GATK_cigar_1.o    # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=/scratch/crs12448/MEVE/Logs/GATK_cigar_1.e
#SBATCH --mail-user=christopher.smaga@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

# Directorty where aligned and sorted BAM files are located
DD="/scratch/crs12448/MEVE/Alignment/HISAT2/BAM"

# Set ouput directory
OD="/scratch/crs12448/MEVE/GATK/MarkDuplicates"

#Load the Genome Anlysis Toolkit
ml  GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
############################################################################################################################################3

cd $OD
OD_2="/scratch/crs12448/MEVE/GATK/SplitNCigarReads"

# # Run SplitNCigarReads to split reads that span introns into separate reads
 for i in S231 S242 S246 S247 S252 S256_2;
 do
  gatk SplitNCigarReads \
       -R /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta \
       -I ${i}_mark_dup.bam  \
       -O $OD_2/${i}_cigar.bam
 done
