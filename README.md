# Systems medicine dissection of chr1q-amp reveals a novel PBX1-FOXM1 axis for targeted therapy in multiple myeloma (Trasanidis N et al., 2021)
### Available preprint: https://www.biorxiv.org/content/10.1101/2021.11.17.469031v2 


![alt text](https://github.com/nikostrasan/PBX1_project/blob/main/Graphical%20abstract_Blood.png)

# Brief Descriptions of Analysis Scripts
## Multiomic clinical data analysis across chr1q
An R script for this analysis is provided in the "chr1q_multivariate_analysis" folder. 
In addition, all input files required for this analysis were downloaded from the MMRF Data Portal (IA15 release)
( https://research.themmrf.org/population/#/ds/611e776c06f8c500013f9448/downloads  ) and reformatted accordingly. 


## RNAseq analysis
Two scripts are provided for RNAseq analysis. 
Script1: Fastq files alignment
Script2: Count Read, normalize and perform differential expression analysis

**Dependencies**:

-STAR 2.6.1d ( https://github.com/alexdobin/STAR )

-STAR Genome [STAR_genome_GRCh38.97] -> Build from scratch using STAR tool

-Genome annotation file [Homo_sapiens.GRCh38.97.gtf] -> Download from Ensembl website https://www.ensembl.org/info/data/ftp/index.html

-biocLite("Rsubread")

-biocLite("GenomicRanges")

-biocLite("biomaRt")

-biocLite("Rsamtools")


## ChIPseq analysis
One bash script is provided to perform ChIPseq analysis, along with a template for input samples.

**Dependencies**:

-fastQC https://github.com/s-andrews/FastQC

-bowtie2 v2.1.0

-samtools v1.9

-bedtools 2.14

-java 

-picard tools

-deeptools 3.1.3

-macs2 2.1.0.20150731

-MSPC.v.2.1

-Homer v4.11 tools  

-BETA-plus (http://cistrome.org/BETA/)


## ATACseq analysis
One bash script is provided to perform ATACseq analysis.

**Dependencies**:

-1. fastQC https://github.com/s-andrews/FastQC

-2. bowtie2 v2.1.0

-3. samtools v1.9

-4. TrimGalore-0.4.3

-5. bedtools 2.14

-6. java 

-7. picard tools

-8. deeptools 3.1.3

-9. macs2 2.1.0.20150731

-10. bedops 3.2.9. v2.4.31

-11. Homer v4.11 tools  




_Data availability_
High-throughput sequencing data generated during this study have been deposited to the Gene Expression Omnibus repository (GEO): 
MMCL ChIP-seq and RNA-seq files (GSE165060); primary MM ATAC-seq and RNA-seq files (GSE153381).


### Please cite: 
_Systems medicine dissection of chromosome 1q amplification reveals oncogenic regulatory circuits and informs targeted therapy in cancer_
Nikolaos Trasanidis, Alexia Katsarou, Kanagaraju Ponnusamy, Yao-An Shen, Ioannis V Kostopoulos, Bien Bergonia, Keren Keren, Paudel Reema, Xiaolin Xiao, Richard M Szydlo, Pierangela MR Sabbattini, Irene AG Roberts, Holger W Auner, Kikkeri N Naresh, Aristeidis Chaidos, Tian-Li Wang, Luca Magnani, Valentina S Caputo, Anastasios Karadimitris
bioRxiv 2021.11.17.469031; doi: https://doi.org/10.1101/2021.11.17.469031

