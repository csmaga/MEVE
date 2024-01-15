#!/bin/bash
#SBATCH --job-name=MEVE_characer
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH --time=72:00:00
#SBATCH --output=/scratch/crs12448/MEVE/Logs/log.%j
#SBATCH --mail-user=crs12448@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-19 

#Output directory for all StringTie analyses
#OD="/scratch/crs12448/MEVE/StringTie"

#cd /scratch/crs12448/MEVE/Alignment/HISAT2/BAM

#samples=/scratch/crs12448/MEVE/Alignment/HISAT2/BAM/sample_list

# # load StringTie
#module load StringTie/2.2.1-GCC-11.3.0

#
#
# # the alignment files from HISAT2 have already been sorted with samtools as BAM files
#for i in S231.bam S242.bam S246.bam S247.bam S252.bam S256_2.bam S263.bam S266_2.bam S280.bam S295.bam S302.bam S303.bam S314.bam S316.bam S317.bam S319.bam S328.bam S336.bam S337.bam S338.bam S344.bam S345.bam S348.bam S350.bam S353.bam S357.bam S359.bam S367.bam S376.bam S380.bam S384.bam S388.bam S391.bam S392.bam S393.bam S406.bam S407.bam S408.bam S416.bam S420.bam S421.bam S422.bam S425.bam S426.bam S427.bam S432.bam S433.bam S435.bam;
#do
#  stringtie -p 8 -G /scratch/crs12448/MEVE/Genome/Amiss.annot.2022.gff -o $OD/assemblies/${i/.bam/.gtf} ${i}
#done

#cd /scratch/crs12448/MEVE/StringTie/assemblies
#
# #Merge all annotation files into one file
#stringtie --merge -i -p 8 -G /scratch/crs12448/MEVE/Genome/Amiss.annot.2022.gff -o stringtie_merged.gtf *.gtf

#mkdir /scratch/crs12448/MEVE/StringTie/GFFcompare
#cd /scratch/crs12448/MEVE/StringTie/GFFcompare
#
#module load GffCompare/0.12.6-GCC-11.2.0

# #Compare assembled annotation to the original GTF
#gffcompare /scratch/crs12448/MEVE/StringTie/assemblies/stringtie_merged.gtf -r /scratch/MEVE/Genome/Amiss.annot.2022.gff -G -M -o gtf_comp
# ## M option indicates gffcompare should ignore single-exon transfags and reference transcripts
#
#module load gffread/0.12.7-GCCcore-11.3.0


#mkdir /scratch/crs12448/MEVE/StringTie/sequences
#cd /scratch/crs12448/MEVE/StringTie/sequences
# Extract fasta sequennces from the merged GTF file for each transcript discovered
#gffread -F -w transcript_seqs.fa -g /scratch/crs12448/MEVE/Genome/Amiss_ref.fasta /scratch/crs12448/MEVE/StringTie/assemblies/stringtie_merged.gtf


#module load BLAST+/2.13.0-gompi-2022a
# Originally did this with BLAST module, but the makeblastdb funciton is not there. Trying with BLAST+ which I think is just a newer BLAST.

#cd /scratch/crs12448/MEVE/StringTie/BLAST

## code below downloads the UniprotKB Swiss-Prot database - a high quality, manually annotated and non-redundant protein sequence database
#wget "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
#gunzip uniprot_sprot.fasta.gz

## code below creates a blast database from the previously downloaded Uniprot database
#makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot_database

#cd /scratch/crs12448/MEVE/StringTie/sequences/seq_bins
# ## Break fasta file with assembled transcripts into smaller parts (10000 sequences each)
#awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%10000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < /scratch/crs12448/MEVE/StringTie/sequences/transcript_seqs.fa


# ## the blastx algorithm will translate the sequence in 3 reading frames in the forward direction and 3 reading frames in the reverse direction to generate the amino acid sequences for the search

# Doing this with arrays to speed things up. 
#sequences=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/crs12448/MEVE/StringTie/sequences/seq_bins/seqs)

# cd /scratch/crs12448/MEVE/StringTie/BLAST

# blastx -query /scratch/crs12448/MEVE/StringTie/sequences/seq_bins/${sequences} -db uniprot_sprot_database -out blasted_${sequences} -outfmt 5 -evalue 0.0001 -num_threads 10

## IGNORE below codes - this was before I knew how to use arrays 
# blastx -query myseq10000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx10000 -outfmt 5 -evalue 0.0001 -num_threads 20
# blastx -query myseq20000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx20000 -outfmt 5 -evalue 0.0001 -num_threads 20
# blastx -query myseq30000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx30000 -outfmt 5 -evalue 0.0001 -num_threads 20
# blastx -query myseq40000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx40000 -outfmt 5 -evalue 0.0001 -num_threads 20
# blastx -query myseq50000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx50000 -outfmt 5 -evalue 0.0001 -num_threads 20
# blastx -query myseq60000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx60000 -outfmt 5 -evalue 0.0001 -num_threads 20
# blastx -query myseq70000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx70000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq80000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx80000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq90000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx90000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq100000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx100000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq110000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx110000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq120000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx120000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq130000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx130000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq140000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx140000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq150000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx150000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq160000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx160000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq170000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx170000 -outfmt 5 -evalue 0.0001 -num_threads 12
# blastx -query myseq180000.fa -db uniprot_sprot_database -out Merged_assembly_Blastx180000 -outfmt 5 -evalue 0.0001 -num_threads 12

## put all blastx output files together

#
#
#
# #
#
cd /scratch/crs12448/MEVE/StringTie/BLAST
# #
#module load Python/3.8.2-GCCcore-8.3.0
#module load SciPy-bundle/2021.05-foss-2019b-Python-3.8.2
#module load matplotlib/3.4.1-foss-2019b-Python-3.8.2
#module load Biopython/1.78-foss-2019b-Python-3.8.2
# # # #
# # #
#python blastxml_to_tabular.py Merged_assembly_Blastx > Merged_assembly_Blastx.tsv

# Doing this with arrays to speed things up. 
#sequences_blasted=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /scratch/crs12448/MEVE/StringTie/BLAST/blasted_seqs)
#python blastxml_to_tabular.py /scratch/crs12448/MEVE/StringTie/BLAST/${sequences_blasted} > /scratch/crs12448/MEVE/StringTie/BLAST/blast_tsv/${sequences_blasted/.fa/.tsv}

# Merge all /tsv files together
cat /scratch/crs12448/MEVE/StringTie/BLAST/blast_tsv/*.tsv > /scratch/crs12448/MEVE/StringTie/BLAST/blast_tsv/Merged_assembly_Blastx.tsv

#cat Merged_assembly_Blastx0.tsv Merged_assembly_Blastx10000.tsv Merged_assembly_Blastx20000.tsv Merged_assembly_Blastx30000.tsv Merged_assembly_Blastx40000.tsv Merged_assembly_Blastx50000.tsv Merged_assembly_Blastx60000.tsv Merged_assembly_Blastx70000.tsv Merged_assembly_Blastx80000.tsv Merged_assembly_Blastx90000.tsv Merged_assembly_Blastx100000.tsv Merged_assembly_Blastx110000.tsv Merged_assembly_Blastx120000.tsv Merged_assembly_Blastx130000.tsv Merged_assembly_Blastx140000.tsv Merged_assembly_Blastx150000.tsv Merged_assembly_Blastx160000.tsv Merged_assembly_Blastx170000.tsv Merged_assembly_Blastx180000.tsv > Merged_assembly_Blastx.tsv


# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx0 > Merged_assembly_Blastx0.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx10000 > Merged_assembly_Blastx10000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx20000 > Merged_assembly_Blastx20000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx30000 > Merged_assembly_Blastx30000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx40000 > Merged_assembly_Blastx40000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx50000 > Merged_assembly_Blastx50000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx60000 > Merged_assembly_Blastx60000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx70000 > Merged_assembly_Blastx70000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx80000 > Merged_assembly_Blastx80000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx90000 > Merged_assembly_Blastx90000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx100000 > Merged_assembly_Blastx100000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx110000 > Merged_assembly_Blastx110000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx120000 > Merged_assembly_Blastx120000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx130000 > Merged_assembly_Blastx130000.tsv
# #python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx140000 > Merged_assembly_Blastx140000.tsv
# python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx150000 > Merged_assembly_Blastx150000.tsv
# python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx160000 > Merged_assembly_Blastx160000.tsv
# python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx170000 > Merged_assembly_Blastx170000.tsv
# python ~/PREE2/blastxml_to_tabular.py /scratch/crs12448/work/PREE2/StringTie/BLAST/Merged_assembly_Blastx180000 > Merged_assembly_Blastx180000.tsv
#
# cat Merged_assembly_Blastx0.tsv Merged_assembly_Blastx10000.tsv Merged_assembly_Blastx20000.tsv Merged_assembly_Blastx30000.tsv Merged_assembly_Blastx40000.tsv Merged_assembly_Blastx50000.tsv Merged_assembly_Blastx60000.tsv Merged_assembly_Blastx70000.tsv Merged_assembly_Blastx80000.tsv Merged_assembly_Blastx90000.tsv Merged_assembly_Blastx100000.tsv Merged_assembly_Blastx110000.tsv Merged_assembly_Blastx120000.tsv Merged_assembly_Blastx130000.tsv Merged_assembly_Blastx140000.tsv Merged_assembly_Blastx150000.tsv Merged_assembly_Blastx160000.tsv Merged_assembly_Blastx170000.tsv Merged_assembly_Blastx180000.tsv > Merged_assembly_Blastx.tsv
# #
# # # NOTE - I found out that the gffread code needed to be adjusted to include the -F option which retains the attributes for transcripts in the defline
# # ## This allowed me to link LOC gene_id's to rna_... id's, so the Blastx results could be reconciled with the corresponding LOC gene
# #
# # # ---Initial approach brainstorm --- #
# # # FILES NEEDED: genehits_count_matrix.csv (target), Merged_assembly_Blastx.tsv (search field)
# # # create new column then for each value in first column of genehits_count_matrix.csv, check if it starts with LOC
# # # if not, enter value of column one for new column
# # # if it does start with LOC, search Merged_assembly_Blastx.tsv first column for vlaue matching 'gene-<match LOC string>'
# # ## of those matches, find one with lowest value in 11th column (e-value) then enter value of column 2 for new column in genehits_count_matrix.csv
# #
# # # the problem is that the counts file used for edgeR analysis just has the LOC gene name, but the Blast hits have various forms of id's based on the StringTie results - sometimes its an MSTRG, sometimes its an XM_ id, sometimes it's a gene-LOC id
# # # I need to be able to link the annotation to the LOC id in the counts file, and at this point, I'm not sure how to do that - it is straightforward for the 'gene-LOC' entries in the Blast hits, but there are only 529 of those and there are more than 10,000 LOC genes in the counts file
# #
# # # ---Actual approach--- #
# # ## first, may need to create new file which only includes the best hits of Merged_assembly_Blastx.tsv for each unique value in first column
# # # code below will do this but if a unique value in column 1 has multiple rows with the minimum value in the 11th column (e-value), the final file includes all rows with this minimum value (i.e., the transcript does not have a single best hit)


cd /scratch/crs12448/MEVE/StringTie/BLAST/BLAST_results
awk '
     NR == FNR {
         if (!($1 in min) || $11 < min[$1])
             min[$1] = $11
         next
     }
     $11 == min[$1]
 ' /scratch/crs12448/MEVE/StringTie/BLAST/blast_tsv/Merged_assembly_Blastx.tsv /scratch/crs12448/MEVE/StringTie/BLAST/blast_tsv/Merged_assembly_Blastx.tsv > Blastx_tophits.tsv

head -n 20 Blastx_tophits.tsv
# # #
# # # # grab all the deflines from the file containing the transcript fasta sequences
# # # # the first value in the defline corresponds to the Blastx output, and the following values (separated by spaces) have information about the corresponding LOC gene
# # # # the new file with the 'index' for each transcript is called transcript_seqs_index.txt
grep '>' /scratch/crs12448/MEVE/StringTie/sequences/transcript_seqs.fa > transcript_seqs_index.txt
head -n 20 transcript_seqs_index.txt
# # # #
# # # # # grab those lines of the index file that have gene_name in them - these are the ones we'll need for the LOC genes
# grep gene_name transcript_seqs_index.txt > transcript_seqs_genename.txt
# head -n 20 transcript_seqs_genename.txt
# # # #
# # # # # modify the index file so the first column doesn't contain the '>' from the fasta file, this way it can be used to match with the first column of the blastX output
# awk -F " " '{print substr($1,2), $2, $3}' transcript_seqs_genename.txt > transcript_names.txt
# head -n 20 transcript_names.txt
# # #
# # # # # if the first column of the transcript name matches the first column of the blastx output, then append that row from the blastx output to the corresponding record in the transcript names file
# join <(sort transcript_names.txt) <(sort Blastx_tophits.tsv) > transcript_matches.txt
# head -n 20 transcript_matches.txt
# # # #
# # # # # modify the previous output so the first column in the gene_name
# awk -F " " '{print substr($2,11), $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' transcript_matches.txt > transcript_matches_id.txt
# head -n 20 transcript_matches_id.txt
# # # #
# # # # At the end, we have a file the contains the different aliases of each transcript along with its best blastx hit
# # # # for example,
# # # #LOC109286440 rna-XR_002094475.1 gene_name=LOC109286440 ref_gene_id=gene-LOC109286440 sp|Q8BHG9|CGBP1_MOUSE 34.118 85 54 2 644 892 78 162 9e-11 45.8
# # # #LOC109286440 rna-XR_002094475.1 gene_name=LOC109286440 ref_gene_id=gene-LOC109286440 sp|Q8BHG9|CGBP1_MOUSE 37.097 62 36 2 454 633 21 81 9e-11 43.5
# # # #LOC109286440 rna-XR_002094476.1 gene_name=LOC109286440 ref_gene_id=gene-LOC109286440 sp|Q8BHG9|CGBP1_MOUSE 34.118 85 54 2 616 864 78 162 9e-11 45.8
# # # #LOC109286440 rna-XR_002094476.1 gene_name=LOC109286440 ref_gene_id=gene-LOC109286440 sp|Q8BHG9|CGBP1_MOUSE 37.097 62 36 2 426 605 21 81 9e-11 43.5
# # # #LOC109286488 rna-XR_002094546.1 gene_name=LOC109286488 ref_gene_id=gene-LOC109286488 sp|Q86TG7|PEG10_HUMAN 32.609 46 31 0 1025 1162 429 474 7e-08 36.2
# # # #LOC109286488 rna-XR_002094546.1 gene_name=LOC109286488 ref_gene_id=gene-LOC109286488 sp|Q86TG7|PEG10_HUMAN 42.424 66 34 2 1203 1400 350 411 7e-08 44.7
# # # #LOC106739010 rna-XR_002094576.1 gene_name=LOC106739010 ref_gene_id=gene-LOC106739010 sp|Q99871|HAUS7_HUMAN 33.708 89 48 2 21 254 93 181 7e-05 47.0
# # # #LOC102576687 rna-XR_002094583.1 gene_name=LOC102576687 ref_gene_id=gene-LOC102576687 sp|Q3T0C8|PDLI2_BOVIN 64.634 82 29 0 1535 1780 1 82 5e-27 116
# # # #LOC102576687 rna-XR_002094584.1 gene_name=LOC102576687 ref_gene_id=gene-LOC102576687 sp|Q3T0C8|PDLI2_BOVIN 64.634 82 29 0 87 332 1 82 4e-31 115
# # # #LOC109286521 rna-XR_002094589.1 gene_name=LOC109286521 ref_gene_id=gene-LOC109286521 sp|A0A1B0GUS0|CS085_HUMAN 47.244 127 57 4 67 441 18 136 9e-11 65.1
# # #
# # # #Only select LOC genes
# # grep LOC transcript_matches_id.txt > transcript_match_id_loc.txt
# # head -n 20 transcript_match_id_loc.txt
# # #
# # # # I will have to associate only one annotation with each gene
# # # # since the top hits from Blastx include multiple annotations if their e-values are equal, I will then select the hit with the greatest value for column 3 = pident (Percentage of identical matches)
# # # #
# # awk '
# #       NR == FNR {
# #           if (!($1 in max) || $3 > max[$1])
# #               max[$1] = $3
# #           next
# #       }
# #       $3 == max[$1]
# #   ' Blastx_tophits.tsv Blastx_tophits.tsv > Blastx_tophits_perc.tsv
# # # # #
#
# join <(sort transcript_names.txt) <(sort Blastx_tophits_perc.tsv) > transcript_match_perc.txt
#
# head -n 20 transcript_match_perc.txt
# # # #
# # # # #
# awk -F " " '{print substr($2,11), $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' transcript_match_perc.txt > transcript_match_perc_id.txt
# head -n 20  transcript_match_perc_id.txt
# # # # #
# # # # # how many unique genes are there? 29,090
# awk -F " " '{print $1}' transcript_match_perc_id.txt | uniq | wc -l
# # # # # how many total lines are there? 144093
# wc -l transcript_match_perc_id.txt
# # #
# # # # there are still multiple annotations for a single gene because there were multiple transcripts per gene from StringTie
# # # # I need to select the best annotation per gene
# # # ## to do this, I will select the annotation with the lowest e-value, then if those are equal, the highest percentage of identical matches, for each gene
# # #
# # # ## first converted this into a csv file
#cp transcript_match_perc_id.txt ./transcript_match_perc_id.csv
#sed -i 's/"//g ; s/\s/,/g' transcript_match_perc_id.csv

# it turns out that the transcript_match_perc_id.txt has different numbers of columns for different types of annotation, notably, the transcripts starting with rna_XM... have an extra column, this causes problems when prioritizing annotations based on the value in column number X
# I downloaded it, and then I removed the extra column in the rows that had it
# then I reuploaded to cluster as transcript_match_perc_id_cor.csv
# #
# awk -F "," '{print $13}' transcript_match_perc_id_cor.csv
# awk -F "," '
#      NR == FNR {
#          if (!($1 in min) || $13 < min[$1])
#              min[$1] = $13
#          next
#      }
#      $13 == min[$1]
#  ' transcript_match_perc_id_cor.csv transcript_match_perc_id_cor.csv > pergene_tophits.csv
#
# head -n 20 pergene_tophits.csv
# # # how many unique genes are
# awk -F "," '{print $1}' pergene_tophits.csv | uniq | wc -l
# # # how many total lines are there?
# wc -l pergene_tophits.csv
# # #38,157
# #
# awk -F "," '
#      NR == FNR {
#          if (!($1 in max) || $5 > max[$1])
#              max[$1] = $5
#          next
#      }
#      $5 == max[$1]
#  ' pergene_tophits.csv pergene_tophits.csv > pergene_tophits_perc.csv
# # # how many unique genes are there? 19,559
# awk -F "," '{print $1}' pergene_tophits_perc.csv | uniq | wc -l
# # # how many total lines are there?
# wc -l pergene_tophits_perc.csv
# #26,875

#I guess it's possible that multiple annotations can have the same e-value and percent identical matches
#I checked and this is exactly what's happening
#I downloaded it, added a new column which is just the row number
#Now, I sort for the minimum row number which will arbitrarily select the first annotation for each gene if it has multiple best hits (from what I've seen, it looks like it's the same match anyway)

# awk -F "," '
#      NR == FNR {
#          if (!($2 in min) || $1 < min[$2])
#              min[$2] = $1
#          next
#      }
#      $1 == min[$2]
#  ' pergene_tophits_perc.csv pergene_tophits_perc.csv > pergene_tophits_perc_final.csv
# #
# # # how many unique genes are there? 19,440
# awk -F "," '{print $1}' pergene_tophits_perc_final.csv | uniq | wc -l
# # # how many total lines are there?
# wc -l pergene_tophits_perc_final.csv
# # #19,440
