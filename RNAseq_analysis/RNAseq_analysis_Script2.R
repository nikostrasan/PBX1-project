#!/bin/r

#Author: Dr. Nikolaos Trasanidis 
######################################################################
################    RNAseq analysis - scripts   ######################
######################################################################

#A. Script2. 
## Dependencies (install packages if not present)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rsubread")
#biocLite("GenomicRanges")
#biocLite("biomaRt")
#biocLite("Rsamtools")

#1.Load packages and input files
library("Rsubread")
library("GenomicRanges")
library("biomaRt")
library("Rsamtools")
getwd()
setwd(paste(getwd(),'/STAR_output',sep=""))

samples_table_PBX1_MM1S<-read.csv("PBX1_MM1S_KD_samples_file_92.csv", header=TRUE)
bamfiles.names <- read.table("Bamfiles_names.table",header=FALSE)
gtf.file_97<-file.path("Homo_sapiens.GRCh38.97.gtf")
bam.files <- file.path(paste0(bamfiles.names[,1], "Aligned.sortedByCoord.out.bam")) 
bam.list <- BamFileList(bam.files, yieldSize = 10000) 

### Obtain Ensembl Human genome annotation
library("biomaRt")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)
ensembl_table <- getBM(attributes = c("ensembl_gene_id", "description", "hgnc_symbol", "chromosome_name", "band"), mart = ensembl)
 
#2. Compute raw counts and export file
fc_PBX <- featureCounts(bam.files, annot.ext = gtf.file_97, isGTFAnnotationFile = TRUE, isPairedEnd = TRUE) 
write.table(as.data.frame(fc_PBX$counts), file="Rawcounts_output.txt", sep="\t")

##LOAD THE Bamfiles_names file!!##

###############################################################################################

#listMarts()
#ensembl=useMart("ensembl")
#listDatasets(ensembl)
#ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#attributes = listAttributes(ensembl)
#ensembl_table <- getBM(attributes = c("ensembl_gene_id", "description", "hgnc_symbol", "chromosome_name", "band"), mart = ensembl)

#3. Deseq2 analysis 
dds_PBX <- DESeqDataSetFromMatrix(fc_PBX$counts, colData = samples_table_PBX1_MM1S, design = ~description)
library(rafalib)
mypar()
dds_PBX<-estimateSizeFactors(dds_PBX)
sizeFactors(dds_PBX)
plot(sizeFactors(dds_PBX), colSums(counts(dds_PBX)))
loggeomeans <- rowMeans(log(counts(dds_PBX)))
hist(log(counts(dds_PBX)[,1]) -loggeomeans, col="grey", main="", xlab="", breaks=50)
log.norm.counts_PBX <- log2(counts(dds_PBX, normalized=TRUE) +1)
rs_PBX<-rowSums(counts(dds_PBX))

#save image 
pdf("Boxplot.Normalizations.pdf")
mypar(1,2)
boxplot(log2(counts(dds_PBX)[rs_PBX > 0,]+1))
boxplot(log.norm.counts_PBX[rs_PBX >0,])
dev.off()

#save image 
pdf("Similarity plot.pdf")
rld <- rlog(dds_PBX)
plot(assay(rld)[,1], assay(rld)[,2], cex=.1)
dev.off()

library(vsn)
meanSdPlot(log.norm.counts_PBX, ranks=FALSE, xlab = "mean", ylab = "sd")
meanSdPlot(assay(rld), ranks=FALSE, xlab = "mean", ylab = "sd")

#save image
pdf("PCA.plot.pdf")
plotPCA(rld, intgroup="description")
dev.off()

#save image here
library(rafalib)
mypar(2,1)

#save image 
pdf("Hclustering.plot.pdf")
plot(hclust(dist(t(assay(rld)))), labels = colData(dds_PBX)$samples)
dev.off()

#Perform differential expression analysis 
#Compare P31 vs scrambled to obtain FoldChange values
design(dds_PBX)
levels(dds_PBX$description)
dds_PBX$description<-relevel(dds_PBX$description, "PBX_P41")
dds_PBX$description<-relevel(dds_PBX$description, "scrbl")
levels(dds_PBX$description)
DDS_P31_vs_scrbl<-DESeq(dds_PBX)
results_P31_vs_scrbl<-results(DDS_P31_vs_scrbl)
summary(results_P31_vs_scrbl)
results_P31_vs_scrbl_annotated<-merge(ensembl_table,results_P31_vs_scrbl,by="ensembl_gene_id")
write.table(results_P31_vs_scrbl_annotated,file="results_P31_vs_scrbl.txt")


#Compare P11 vs scrambled to obtain FoldChange values
design(dds_PBX)
levels(dds_PBX$description)
dds_PBX$description<-relevel(dds_PBX$description, "PBX_P31")
dds_PBX$description<-relevel(dds_PBX$description, "scrbl")
levels(dds_PBX$description)
DDS_P11_vs_scrbl<-DESeq(dds_PBX)
results_P11_vs_scrbl<-results(DDS_P11_vs_scrbl)
summary(results_P11_vs_scrbl)
results_P11_vs_scrbl_annotated<-merge(ensembl_table,results_P11_vs_scrbl,by="ensembl_gene_id")
write.table(results_P11_vs_scrbl_annotated,file="results_P11_vs_scrbl.txt")


#Obtain q-values from LTR test across scrambled, P11 and P31. 
dds_LRT<- DESeq(dds_PBX, test="LRT", full=~description, reduced=~1) 
resdds_LRT<-as.data.frame(results(dds_LRT)) 
summary(results(resdds_LRT))
results_LTR_annotated<-merge(ensembl_table,resdds_LRT,by="ensembl_gene_id")
write.table(resdds_LRT,file="results_LRT_test.txt")










