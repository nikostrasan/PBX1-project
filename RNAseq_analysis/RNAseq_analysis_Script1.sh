#!/bin/bash

#Author: Dr. Nikolaos Trasanidis 
######################################################################
################    RNAseq analysis - scripts   ######################
######################################################################


#PBS -l walltime=70:30:00
                #PBS -l select=1:mem=32gb:ncpus=4
                #PBS -j oe
                #PBS -o $RNAseq_no2
                
#Dependencies:
#STAR 2.6.1d https://github.com/alexdobin/STAR
#STAR Genome [STAR_genome_GRCh38.97] 
#Genome annotation file [Homo_sapiens.GRCh38.97.gtf] 

export STAR_genome="./STAR_genome_GRCh38.97/" # Build from scratch using STAR tool
export Genome_annot_file="Homo_sapiens.GRCh38.97.gtf" # Download from Ensembl website https://www.ensembl.org/info/data/ftp/index.html
export fastq_files="$Fastq_files_directory"

#1. Replace the dots in filenames with underscores & replace back the ".fastq.gz"
cd ./fastq_files/
rename - _ *.fastq.gz
rename - _ *.fastq.gz
rename - _ *.fastq.gz
rename - _ *.fastq.gz
rename - _ *.fastq.gz
rename - _ *.fastq.gz
cd ..

#2. Read all fastq.gz file names, write them to a table file and reformat 
for zipfile in ./fastq_files/*.fastq.gz; do
file=$(basename $zipfile .gz)
gunzip $zipfile
echo $file >> RNAseq_sample_names.txt
done

sed 'n; d' RNAseq_sample_names.txt > k1.txt
sed '1d; n; d' RNAseq_sample_names.txt > k2.txt
paste k1.txt k2.txt | pr -t > RNAseq_sample_names.table


#3. run STAR for fastq files pairs, as designated in the "names" file
mkdir STAR_output

while read -r qfile1 qfile2 ; do
echo $qfile1 $qfile2 
outname=$(basename "$qfile1" | cut -d. -f1)

STAR runThreadN 4 --genomeDir $STAR_genome --readFilesIn ./fastq_files/$qfile1 ./fastq_files/$qfile2 --twopassMode Basic --outSAMtype BAM SortedByCoordinate  --outFileNamePrefix $outname --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000  --alignMatesGapMax 1000000 limitGenomeGenerateRAM 35019136213 -p 4 --limitBAMsortRAM 35019136213 -n 100000 --outBAMsortingThreadN 4

done < RNAseq_sample_names.table

#6. Write table with all "Aligned.sortedByCoord.out.bam" files for downstream analysis

for bamfile in *Aligned.sortedByCoord.out.bam; do
mv $bamfile ./STAR_output/
bamname=$(basename $bamfile Aligned.sortedByCoord.out.bam)
echo $bamname >> Bamfiles_names.table; done

#7. Run the Rscript for Rsubreads analysis
Rscript RNAseq_analysis_script2.R


