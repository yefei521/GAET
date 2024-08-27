
---

#  Genome Annotation by in-depth Epigenetic and Transcriptomic profiling (GAET)

![20240827145514](https://github.com/user-attachments/assets/0a26dd8d-4710-4f71-a97c-919ac66c4c7d)


This repository contains the scripts, tools, and workflows for the comprehensive gene annotation of *Tetrahymena thermophila* using various types of omics data, including transcriptomic, epigenetic, and sequencing data. The pipeline integrates gene model optimization, UTR annotation, alternative splicing (AS) isoform annotation, and manual curation to produce high-quality gene models.

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Workflow](#pipeline-workflow)
3. [Dependencies](#dependencies)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Output](#output)
7. [Contributing](#contributing)
8. [License](#license)

## Overview

This pipeline is designed to optimize and annotate gene models in *Tetrahymena thermophila*. The process involves several steps, including RNA-seq and epigenetic data integration, UTR and regulatory element analysis, AS isoform annotation, and manual genome curation. The final outputs are refined gene models and protein function annotations.

## Pipeline Workflow

### A. Gene Optimization with Transcriptomic Data

1. **Input Data:**
   - RNA-seq 
   - Nanopore DRS
   - Strand-specific RNA-seq

```bash
git clone https://github.com/username/tetrahymena-gene-annotation.git
cd tetrahymena-gene-annotation

    ~/biosoft/TrimGalore-0.6.6/trim_galore -q 20 -length 20 --fastqc -j 4  --paired $j $R2 -o 02Clean_data/ >> 02Clean_data/Trim_Galore.log 2>&1

    if [ ! -d "03Hisat2_align" ];then   mkdir -p 03Hisat2_align ; fi
    hisat2 -p 6 --dta --score-min L,5,-0.4 --sp 1,0 --mp 3,1  --summary-file 03Hisat2_align/${name}_summary.txt -x Transcriptome_Assemblly-221208/BamFile/00Hisat2_index/1-Genome_assembly_new.fasta -S 03Hisat2_align/${name}.sam  -1 02Clean_data/${name}_1_val_1.fq.gz -2 02Clean_data/${name}_2_val_2.fq.gz  2>&1 >>03Hisat2_align/error.ht2.log
    java -jar /home/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar  SortSam -I 03Hisat2_align/${name}.sam -O 03Hisat2_align/${name}_sorted.bam -SORT_ORDER coordinate > 03Hisat2_align/${name}_sort.log 2>&1 
    java -jar /home/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar  MarkDuplicates -I 03Hisat2_align/${name}_sorted.bam -O 03Hisat2_align/${name}_rmdup.bam -REMOVE_DUPLICATES true -M 03Hisat2_align/${name}_marked_dup_metrics.txt > 03Hisat2_align/${name}_rmdup.log 2>&1
     rm 03Hisat2_align/${name}_sorted.bam
     rm 03Hisat2_align/${name}.sam
     perl 01Check_gff-file_format.pl
     perl 02Check_ATG2TGA.pl
     perl 2Merge_orf1_allORF_to_UTR-gff.pl
     perl 2Merge_orf1_to_UTR-gff.pl
     perl 03Remove_error_gene.pl
     perl 04_2Check_CDS_overlap.pl
     perl 04_3Remove_Error_CDS-overlap-gene.pl
     perl 04Check_gene_overlap.pl
     perl 4Filter_main_Nanogff.pl
     perl 05_1Remove_intron2long.pl
     perl 05Find_new_Gap-gene_RNAseq.pl
     perl 5Sort_gff.pl
     perl 6Remove_AS_transcript.pl
     perl 6Sort_gff.pl
     perl 7Merge_4mainGFF.pl
     perl 7Sort_gff.pl
     perl 8Find_new_Gap-gene_RNAseq.pl
     perl 9Merge_4GapNewGene2GFF.pl
     perl 9SeleceGFF_fromGTF.pl
     perl 9Sort_gff.pl
     perl 10Merge_7-9_sortdGFF.pl
     perl Change_OriError.pl
     perl Check_gene_overlap.pl
     perl Remove_LongGapExon_gene.pl
```

2. **Process:**
   - Assemble transcriptome from RNA-seq data.
   - Compare intron-exon junctions across data sources.
   - Perform iterative gene optimization leading to Draft Gene Model v1 and v2.

### B. Gene Optimization with Epigenetic Data

1. **Input Data:**
   - ATAC-seq
   - Epigenetic markers (H3K4me3, H2A.Z, 6mA)

2. **Process:**
   - Predict transcription start sites (TSS) using Random Forest models.

```bash
#Cluster mark 

library("fpc")  #eps: nucleosome + linker = 200，或 half of linker + half of nucleosome = 100

#H2A.Z
data<-read.delim("H2A.Z_group.gff3",sep="\t", header = F)
#model2 <- dbscan(data, eps=147, MinPts=5)
for (i in unique(data$V1)){  a<-subset(data,V1==i); a1<-cbind(a$V4,V2=1); colnames(a1)<-c("V1","V2"); model3 <- dbscan(a1, eps=200, MinPts=3); a1<-as.matrix(a1); a2<-cbind(model3$cluster,i,a1[,1]); write.csv(a2,file=paste("Z:/H2A.Z/",i,".csv"，collapse = NULL, sep=""), sep="\t") }  
data<-read.delim("Z:/H2A.Z_group.gff3",sep="\t", header = F)

#nucleosome
data<-read.delim("Z:/Nucleosome.srt.rmdup.bed",sep="\t", header = F) 
for (i in unique(data$V1)){  a<-subset(data,V1==i); a1<-cbind(a$V2,V2=1); colnames(a1)<-c("V1","V2"); model3 <- dbscan(a1, eps=200, MinPts=3); a1<-as.matrix(a1); a2<-cbind(model3$cluster,i,a1[,1]); write.csv(a2,file=paste("Z:/nucleosome/",i,".csv"，collapse = NULL, sep=""), sep="\t") }

#H3K4me3
H3K4me3<-read.delim("Z:/H3K4me3_group.gff3",sep="\t", header = F)
for (i in unique(data$V1)){  a<-subset(data,V1==i); a1<-cbind(a$V4,V2=1); colnames(a1)<-c("V1","V2"); model3 <- dbscan(a1, eps=200, MinPts=3); a1<-as.matrix(a1); a2<-cbind(model3$cluster,a1[,1]); write.csv(a2,file=paste("Z:/H3K4me3/",i,".csv"，collapse = NULL, sep=""), sep="\t") }

#6mA
data<-read.delim("Z:/veg_6mA.sorted.gff",sep="\t", header = F) 
for (i in unique(data$V1)){  a<-subset(data,V1==i); a1<-cbind(a$V4,V2=1); colnames(a1)<-c("V1","V2"); model3 <- dbscan(a1, eps=200, MinPts=3); a1<-as.matrix(a1); a2<-cbind(model3$cluster,i,a1[,1]); write.csv(a2,file=paste("Z:/6mA/",i,".csv"，collapse = NULL, sep=""), sep="\t") }

#Tranining and prediction using Random Forest

library(randomForest)
library(caret)
library(pROC)
set.seed(1234)
colnames(TSS)<-c("6mA","H3K4me3","H2A.Z","nucleosome","type")

TSS<-read.csv("TSS-input.txt",sep = "\t",header = 1, row.names = 1)
TES<-read.csv("TES-input.txt",sep = "\t",header = 1, row.names = 1)
TSS_TES<-rbind(TSS,TES)
trains <- createDataPartition(y = TSS_TES$type,p = 0.75,list = F)
traindata <- TSS_TES[trains,]
testdata <- TSS_TES[-trains,]
rf.train <- randomForest(as.factor(type)~.,data = traindata,na.action = na.roughfix,importance = TRUE, ntree=500, mtry=1)

plot(rf.train,main="ERROR&TREES")

importance(rf.train)

trainpredprob<-predict(rf.train, newdata = traindata, type="prob")

trainroc<-roc(response=traindata$type, predictor=trainpredprob[,2])

plot(trainroc, print.auc=TRUE, auc.polygon=TRUE, grid=T, max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=T, legacy.axes=T, bty="l")

bestp <- trainroc$threshold[which.max(trainroc$sensitivities + trainroc$specificities -1)]

trainpredlab <- as.factor(ifelse(trainpredprob[,2] > bestp ,1,0))

confusionMatrix(data = as.factor(trainpredlab),  reference = as.factor(traindata$type), positive = '1', mode = 'everything')

testpredprob <- predict(rf.train, newdata = testdata,type = 'prob')

testpredlab <- as.factor(ifelse(testpredprob[,2] > bestp ,1,0))

confusionMatrix(data = as.factor(testpredlab),  reference = as.factor(testdata$type),  positive = '1', mode = 'everything') 

testroc <- roc(response = testdata$type,predictor = testpredprob[,2]) 

plot(trainroc,print.auc = T,grid = c(0.1,0.2),auc.polygon = F,max.auc.polygon = T,main = "random Forest ROC", grid.col = c("green","red"))
plot(testroc,print.auc = T, print.auc.y = 0.4,add = T, col = 'red')
legend("bottomright",legend = c("traindata","testdata"), col = c(par('fg'),'red'), lwd = 2, cex = 0.9)


rf.test <- predict(rf.train,traindata,type="response")
trn_pred<-ifelse(predict(rf.train, type = "response")> 0.3, 1, 0)
trn_tab <- table(predicted = trn_pred, actual = traindata$type)
tst_pred <- ifelse(predict(rf.train, newdata = testdata, type = "response") > 0.3, 1, 0)
tst_tab <- table(predicted = tst_pred, actual = testdata$type)
calc_class_err(actual = teat$type, predicted = tst_pred)
confusionMatrix(trn_tab, positive = "1")
test_roc <- roc(testdata$type , test_prob, plot = TRUE, print.auc = TRUE)

```
     
   - Validate and optimize gene models.

   - Generate Draft Gene Model v3.

### C. UTR Annotation and Regulatory Element Analysis

1. **Input Data:**
   - ATAC-seq
   - Nanopore RNA data

2. **Process:**
   - Annotate 5' and 3' UTRs.
   - Identify regulatory elements like promoters, TATA boxes, and poly-A signals.
   - Update to Draft Gene Model v4.

### D. Manual Curation and Genome Polish

1. **Input Data:**
   - Whole-genome sequencing data
   - Public protein databases

2. **Process:**
   - Manual curation of genes.
   - Align sequencing data to polish the genome.
   - Finalize gene models in Draft Gene Model v5.

### E. Alternatively Spliced Isoform Annotation

1. **Input Data:**
   - RNA-seq (growth, starvation, conjugation)
   - Nanopore DRS

2. **Process:**
   - Identify and annotate alternative splicing events.
   - Generate updated gene annotations in TGD2024.

## Dependencies

- Python 3.x
- R
- GATK
- Samtools
- Nanopore basecalling software
- [List any other specific tools or libraries used]

## Installation

Clone the repository:

```bash
git clone https://github.com/username/tetrahymena-gene-annotation.git
cd tetrahymena-gene-annotation
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Usage

1. **Prepare the input data**: Ensure that all required data files (RNA-seq, ATAC-seq, etc.) are available in the specified directories.
2. **Run the pipeline**: Execute the main script that coordinates the workflow.

```bash
python main_pipeline.py --input data_directory --output results_directory
```

3. **Monitor progress**: Logs and intermediate files will be saved in the output directory.

## Output

- **Gene Models**: Optimized gene models in GFF3 format.
- **UTR Annotations**: Annotated UTR regions.
- **Alternative Splicing Events**: Detailed AS isoform annotations.
- **Polished Genome**: Updated and polished genome sequence.
- **Protein Function Annotations**: Functionally annotated proteins linked to gene models.

## Tetrahymena_Genome_annotation_V2024
The appendix file of paper --"Comprehensive genome annotation of the model ciliate Tetrahymena thermophila by in-depth epigenetic and transcriptomic profiling." Doi: https://www.biorxiv.org/content/10.1101/2024.01.31.578305v1

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to discuss your changes or suggestions.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---
