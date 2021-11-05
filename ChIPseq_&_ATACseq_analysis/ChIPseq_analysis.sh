#!/bin/bash

#Author: Dr. Nikolaos Trasanidis 
#######################################################################
################    ATACseq analysis - scripts   ######################
#######################################################################

#Dependencies
#1. fastQC https://github.com/s-andrews/FastQC
#2. bowtie2 v2.1.0
#3. samtools v1.9
#4. bedtools 2.14
#5. java 
#6. picard tools
#7. deeptools 3.1.3
#8. macs2 2.1.0.20150731
#9. MSPC.v.2.1
#10. Homer v4.11 tools  
#11. BETA-plus (http://cistrome.org/BETA/)

#Before starting, build hg19 and hg38 genomes using bowtie2 tools prior to analysis
export fastq_files="$Fastq_files_directory"

#1. Run fastQC 
mkdir fastqc_reports
for file in ./fastq_files/*.fastq.gz; do
	name=$(basename $file .fastq)
	#echo $file; echo $name
fastqc $file -o './fastqc_reports/'
done 


# 2. Run bowtie2 for fastq samples
mkdir Analysis
for fastqfile in ./fastq_files/*.fastq.gz; do
	name2=$(basename $fastqfile .fastq)
	#echo $file; echo $name
bowtie2 -x /work/nt112/Genome_files/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -U $fastqfile  -p 8 -S ./Analysis/$name2.sam 2> ./Analysis/$name2.log
done 




##3. Run main analysis
cd ./Analysis/
> Stats_bamfiles_depth.table
for file3 in *.sam; do
	name3=$(basename $file3 .sam)
	echo $file; echo $name
	
	samtools view -bSo $name3.bam $file3 
	samtools sort $name3.bam -o $name3.sorted.bam
	picard MarkDuplicates INPUT=$name3.sorted.bam OUTPUT=$name3.nodup.bam REMOVE_DUPLICATES=true METRICS_FILE=$name3nodup_metrics_file.txt
	picard BuildBamIndex INPUT=$name3.nodup.bam
	samtools sort  $name3.nodup.bam -o $name3.nodup.sorted.bam
 	samtools index $name3.nodup.sorted.bam
	samtools view -c $name3.nodup.sorted.bam >> Stats_bamfiles_depth.table
	bamCoverage --bam $name3.nodup.sorted.bam --binSize 1 -e 200 --normalizeUsing RPKM  -o $name3.nodup.sorted_RPKM_bin1.bw
done 




##4. Peak calling  
mkdir peak_calling

while read -r treatment control ; do
echo $treatment $control 
outname=$(basename "$treatment" .nodup.sorted.bam)
#echo $outname

#a. For TFs
macs2 callpeak -t $treatment -c $control -g hs -n ./peak_calling/$outname.vs_Inp_0.05 -B -q 0.01 --verbose 4 --SPMR --call-summits
#b. For Histone Marks
macs2 callpeak -t $treatment -c $control -g hs -n ./peak_calling/$outname.vs_Inp_0.01 -B --broad --broad-cutoff 0.01 --verbose 4 --SPMR 

bamCompare --bamfile1 $treatment --bamfile2 $control --binSize 1 --normalizeUsing RPKM --scaleFactorsMethod None -e 200 --operation ratio -o peak_calling/$outname.FC_over_Input_bin1.bw
done < ChIPseq_table_samples_match.table





#5. Obtain reproducible peaks and signal tracks from replicates
#example for PBX1_ChIPseq_MM1S
mono '/home/MSPC.v.2.1/MSPC.exe' -i ./peak_calling/'PBX1_rep1_ChIPseq_MM1S.narrowPeak' -i ./peak_calling/'PBX1_rep2_ChIPseq_MM1S.narrowPeak' -r biological -s 1E-8 -w 1E-4 –c 2 

samtools merge PBX1_rep12_ChIPseq_MM1S.merged.bam PBX1_rep1_ChIPseq_MM1S.nodup.sorted.bam PBX1_rep2_ChIPseq_MM1S.nodup.sorted.bam 
samtools merge PBX1_rep12_Input_MM1S.merged.bam PBX1_rep1_Input_MM1S.nodup.sorted.bam PBX1_rep2_Input_MM1S.nodup.sorted.bam

bamCompare --bamfile1 PBX1_rep12_ChIPseq_MM1S.merged.bam --bamfile2 PBX1_rep12_Input_MM1S.merged.bam --binSize 1 --normalizeUsing RPKM --scaleFactorsMethod None -e 200 --operation ratio -o peak_calling/$outname.FC_over_Input_bin1.bw




#6. Motif analysis, genomic regions annotation and super-enhancer calling
#Use Homer v4.11 tools for these tasks

findMotifsGenome.pl ./peak_calling/$file.narrowPeak hg38 MotifOutput/ -size given -mask
annotatePeaks.pl ./peak_calling/$file.narrowPeak hg38 > ./peak_calling/$file.narrowPeak.annotated.txt

#Super-enhancer calling, example for H3K27ac analysis
#a. First make Tag Directories for CHIP AND INPUT files you want to use#
makeTagDirectory "MM1S_ChIP_H3K27ac/" "MM1S_H3K27ac.nodup.sorted.bam"
makeTagDirectory "MM1S_ChIP_Input/" "MM1S_Input.nodup.sorted.bam"

#b. Then perform super-enhancer analysis
findPeaks "MM1S_H3K27ac_ChIP" -i MM1S_Input_ChIP -style super -o MM1S_H3K27ac_ChIP.SuperEnhancers.txt -typical MM1S_H3K27ac_ChIP.typical_enhancers.txt -L 0

sed '/^#/ d' MM1S_ChIP_H3K27ac.SuperEnhancers.txt | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$4,$1,".",$6)}' - > MM1S_ChIP_H3K27ac.SuperEnhancers.bed
sed '/^#/ d' MM1S_ChIP_Input.typical_enhancers.txt | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$4,$1,".",$6)}' - > MM1S_ChIP_Input.typical_enhancers.bed




#7. Cistrome-Transcriptome integration using BETA-plus package
#Use reproducible peaks from MSPC output along with knockdown RNAseq data
BETA plus -p $PEAKFILE.bed -e $Diff_expression_file.xls -o $OUTDIR/  -k BSF -g hg38 --gs 'hg38.fa' --reference='ENSEMBL_GRCh38.97_biomart_forBETA.txt' --info 1,2,3 --da 1 --df 0.1 -c 1e-3 











exit


