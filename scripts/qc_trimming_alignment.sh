#!/bin/bash
#SBATCH --job-name=MEVE_alignment
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

if [ ! -d $OUTDIR/Genome ]
then
    echo 'There is no genome folder.'
fi



echo
echo "******************       BEGIN PIPELINE       **********************"
echo

date

echo 'Project directory = ' $OUTDIR
echo
echo 'Sample: ' $sample
echo
echo 'Raw FastQC'

##load modules for downloading/trimming experiment reads and alignment
echo
echo 'Loading relevant modules...'
echo

ml FastQC/0.11.9-Java-11

echo
echo 'loading modules complete'
echo

#Create directory for raw read fastqc files
echo
echo 'Creating directory for raw read QC...'
echo
if [ ! -d $OUTDIR/RawReadQC ]
then
    mkdir -p $OUTDIR/RawReadQC
fi
echo
echo 'directory created.'
echo

#fastqc -o $OUTDIR/RawReadQC -t 10 $OUTDIR/Data/Raw/${sample}/*.gz

echo 
echo 'Trimming raw reads...'

#Create directory for trimmed reads
echo
echo 'Creating directory for trimmed reads'
echo

if [ ! -d $OUTDIR/TrimmedReads ]
then
    mkdir -p $OUTDIR/TrimmedReads
fi
echo
echo 'directory created.'
echo

#Create directory for trimmed QC files
echo
echo 'Creating directory for trimmed QC...'
echo

if [ ! -d $OUTDIR/TrimmedReadQC ]
then
    mkdir -p $OUTDIR/TrimmedReadQC
fi

echo
echo 'directory created.'
echo

echo
echo 'Load modules for trimming...'
echo

ml Trim_Galore/0.6.7-GCCcore-11.2.0
ml Python/3.10.8-GCCcore-12.2.0
ml pigz/2.7-GCCcore-11.3.0
ml cutadapt/4.5-GCCcore-11.3.0


##trim raw reads

echo
echo 'Trimming raw reads and performing FastQC...'
echo

#trim_galore --cores 4 --fastqc --fastqc_args "--outdir $OUTDIR/TrimmedReadQC" -stringency 3 -o $OUTDIR/TrimmedReads --paired  $OUTDIR/Data/Raw/${sample}/${sample}_1.fq.gz $OUTDIR/Data/Raw/${sample}/${sample}_2.fq.gz
fastqc -o $OUTDIR/TrimmedReadQC -t 10 $OUTDIR/TrimmedReads/*.gz

echo
echo 'trimming complete.'
echo

echo
echo 'Concatenate all FastQC files using MultiQC...'
echo 

multiqc /$OUTDIR/TrimmedQC

date


##Align reads to reference genome

echo
echo "******************        Begin Alignment      ********************"
echo

echo 
echo 'Loading modules...'
echo

ml HISAT2/3n-20201216-gompi-2022a

#Create directory for alignments
if [ ! -d $OUTDIR/Alignment]
then
    mkdir -p $OUTDIR/Alignment
fi

if [ ! -d $OUTDIR/Alignment/HISAT2/SAM]
then
    mkdir -p $OUTDIR/Alignment/HISAT2/SAM
fi

if [ ! -d $OUTDIR/Alignment/HISAT2/Genome_Index]
then
    echo 'Genome index directory not found.'
fi

echo 'Aligning...'

hisat2 -x $OUTDIR/Alignment/HISAT2/Genome_Index/Index/Amiss.index.ref.hisat2  -p 10 --rna-strandness FR --dta -q -1 $OUTDIR/TrimmedReads/${sample}_1_val_1.fq.gz -2 $OUTDIR/TrimmedReads/${sample}_2_val_2.fq.gz -S $OUTDIR/Alignment/HISAT2/SAM/${sample}.sam --summary-file /$OUTDIR/MEVE/Alignment/HISAT2/Stats/${sample}_HISAT2_alignment_summary

echo
echo 'alignment complete'
echo

date

## load samtools and sort bam files
if [ ! -d $OUTDIR/Alignment/HISAT2/BAM]
then
    mkdir $OUTDIR/Alignment/HISAT2/BAM
fi

echo 
echo 'Load modules...'
echo

ml SAMtools/1.16.1-GCC-11.3.0

echo
echo 'Converting SAM to BAM and sorting...'
samtools sort -@ 10 $OUTDIR/Alignment/HISAT2/SAM/${sample}.sam -o $OUTDIR/Alignment/HISAT2/BAM/${sample}.bam

echo 'Sorting and conversion complete.'

date

echo
echo "******************       FINISH PIPELINE        **********************"
echo
