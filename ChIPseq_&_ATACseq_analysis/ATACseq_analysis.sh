#!/bin/bash

#Author: Dr. Nikolaos Trasanidis 
#######################################################################
################    ATACseq analysis - scripts   ######################
#######################################################################

#Dependencies
#1. fastQC https://github.com/s-andrews/FastQC
#2. bowtie2 v2.1.0
#3. samtools v1.9
#4. TrimGalore-0.4.3
#5. bedtools 2.14
#6. java 
#7. picard tools
#8. deeptools 3.1.3
#9. macs2 2.1.0.20150731
#10. bedops
#11. Homer v4.11 tools  

#Before starting, build hg19 and hg38 genomes using bowtie2 tools prior to analysis
#Also download HOCOMOCOv11 Human motifs from https://hocomoco11.autosome.ru/downloads_v11
export fastq_files="$Fastq_files_directory"

#1. Read all file names and write them in table file
> ATACseq_sample_names.txt
for zipfile in ./fastq_files/*.fastq.gz; do
file=$(basename $zipfile .gz)
gunzip $zipfile
echo $zipfile "unzipped"
echo $file >> ATACseq_sample_names.txt
done

sed 'n; d' ATACseq_sample_names.txt > a1.txt
sed '1d; n; d' ATACseq_sample_names.txt > a2.txt
paste a1.txt a2.txt | pr -t > ATACseq_sample_names.table
rm a1.txt a2.txt

#2.Fastq files trimming 
#Trims nextera reads, all bps with Q<30 AND removes completely reads <20bp
mkdir fastqc_reports
while read -r qfile1 qfile2 ; do
echo $qfile1 $qfile2 
perl trim_galore --paired  -q 30 --nextera $qfile1 $qfile2 -o trimmed_fastq_files/ 
done < ATACseq_sample_names.table

#3.Quality control using fastQC 
for file in ./trimmed_fastq_files/*.fq.gz; do
	name=$(basename $file .fq)
	echo $file; echo $name
fastqc $file -o './fastqc_reports/'
done 

sed 's/fastq/fq/g' ATACseq_sample_names.table > ATACseq_sample_names2.table

#4. Run bowtie2 for fastq files alignment
mkdir Analysis_hg19 #example for hg19 here

while read -r qfile1 qfile2 ; do
echo $qfile1 $qfile2 
outname=$(basename "$qfile1" | cut -d_ -f1)
bowtie2 -x /work/nt112/Genome_files/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ./trimmed_fastq_files/$qfile1 -2 ./trimmed_fastq_files/$qfile2  -p 8 -S ./Analysis_hg19/$outname.sam -k 4 -X 2000 2> ./Analysis_hg19/$outname.log
done <ATACseq_sample_names2.table


#5. Run main analysis
cd Analysis_hg19
> Stats_bamfiles_depth.table
for file3 in *.sam; do
	name3=$(basename $file3 .sam)
	#echo $file; echo $name
    
    echo 'start analysis of sample' $name3
	samtools view -bSo $name3.bam $file3 
	samtools sort $name3.bam -o $name3.sorted.bam
    samtools index $name3.sorted.bam
	samtools view -F 524 -f 2  -b $name3.sorted.bam | samtools view -F 1804 -f 2 -b -o $name3.sorted.v1.bam -
    samtools sort $name3.sorted.v1.bam -o $name3.v1.resorted.bam
    picard MarkDuplicates INPUT=$name3.v1.resorted.bam OUTPUT=$name3.v1.resorted_nodup.bam REMOVE_DUPLICATES=true METRICS_FILE=$name3.nodup_metrics_file.txt
	picard BuildBamIndex INPUT=$name3.v1.resorted_nodup.bam
	samtools view -h $name3.v1.resorted_nodup.bam | grep -v chrM | grep -v chrEBV | grep -v chrUn_* | samtools view -bS - >$name3.v1.resorted_nodup_noChrM.bam 
    samtools sort  $name3.v1.resorted_nodup_noChrM.bam -o $name3.v1_nodup_noChrM_sorted.bam
 	samtools index $name3.v1_nodup_noChrM_sorted.bam
	count=$(samtools view -c $v1_nodup_noChrM_sorted.bam)
    printf "%s\t%s\n" $name $count >> Stats_bamfiles_depth.table
	bamCoverage --bam $name3.v1_nodup_noChrM_sorted.bam --binSize 1 -e 200 --normalizeUsing RPKM  -o $name3.v1_nodup_noChrM_sorted.bw
	  
done 


#6. Subsample bam files based on depth
minimum=$(sort -k2 -n Samples_depth_stats.txt | head -1 | awk '{print $2}')
echo The minimum depth is $minimum 

> Samples_depth_stats2.txt
for bamfile in *.nodup.sorted_chrUnrm.bam; do
    Name=$(basename $bamfile .bam)
    Count=$(samtools view -c $bamfile)	
    seed_cover=$(echo "scale=40; 5 + ( $minimum / $Count )" | bc)
    samtools view -s $seed_cover -b $bamfile > $Name.subsampl.bam
    samtools index $Name.subsampl.bam
	bedtools bamtobed -i $Name.subsampl.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > $Name.subsampl.bed.gz
    printf "%s\t%s\n" $Name $Count >> Samples_depth_stats2.txt
    echo File $bamfile is finished 
done


#7. Peak calling using macs2
mkdir peak_calling
for file4 in *.shifted.bed.gz; do
	name4=$(basename $file4 .v1_nodup_noChrM_sorted.shifted.bed.gz)
	echo $file; echo $name

    macs2 callpeak -t $file4 -f BED -n $./peak_calling/$name4.narrowpeaks_q001 -g hs -q 0.01 --nomodel --shift -100 --extsize 200 -B --SPMR --call-summits
	#macs2 callpeak -t $file4 -f BED -n $./peak_calling/$name4.broadpeaks_q001 -g hs -q 0.01 --nomodel --shift -100 --extsize 200 --broad --broad-cutoff 0.01 --keep-dup all
done



#8. Differential footprinting 
#Example analysis for 1q_pos versus 1q_neg 

#A. Prepare input files
samtools merge chr1q_pos.bam $file1.subsampl.bam $file2.subsampl.bam $file3.subsampl.bam #...
samtools merge chr1q_neg.bam $file4.subsampl.bam $file5.subsampl.bam $file6.subsampl.bam #...

samtools sort chr1q_pos.bam -o chr1q_pos.sorted.bam && samtools index chr1q_pos.sorted.bam
samtools sort chr1q_neg.bam -o chr1q_neg.sorted.bam && samtools index chr1q_neg.sorted.bam

bedops --merge $file1.subsampl.bed $file2.subsampl.bed $file3.subsampl.bed > chr1q_pos.bed
bedops --merge $file4.subsampl.bed $file5.subsampl.bed $file6.subsampl.bed > chr1q_neg.bed
bedops -m chr1q_pos.bed chr1q_neg.bed > chr1q_consensus.bed

#B. Perform Differential footprinting
wellington_bootstrap.py 'chr1q_pos.bam' 'chr1q_neg.bam' 'chr1q_consensus.bed' 'chr1q_pos_specif_FTs.bed' 'chr1q_neg_specif_FTs.bed' -fdr 0.1 -fdrlimit -5 -A 

#C. Filter footprints based on wellington score (10)
awk '$5>10' chr1q_pos_specif_FTs.bed > chr1q_pos_specif_FTs_filtered.bed
awk '$5>10' chr1q_neg_specif_FTs.bed > chr1q_neg_specif_FTs_filtered.bed

#D. Reformatting and homer motifscan
### re-format files before input to homer ### 
mkdir Homer_motifscan
awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4":"$1":"$2":"$3":"$5,$5,$6)}' chr1q_pos_specif_FTs_filtered.bed > chr1q_pos_specif_FTs_filtered.reform.bed  
awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4":"$1":"$2":"$3":"$5,$5,$6)}' chr1q_neg_specif_FTs_filtered.bed > chr1q_neg_specif_FTs_filtered.reform.bed  

perl findMotifsGenome.pl chr1q_pos_specif_FTs_filtered hg38 -size given -mask -find HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif> chr1q_pos_TFmotif_Footprints.txt
perl findMotifsGenome.pl chr1q_neg_specif_FTs_filtered hg38 -size given -mask -find HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif> chr1q_neg_TFmotif_Footprints.txt

#Filter TFs based on expression and proceed with downstream analysis 


exit


