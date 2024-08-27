#!/bin/bash
#include<stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
[ -e /tmp/fd1 ] || mkfifo /tmp/fd1
exec 3<>/tmp/fd1
rm -rf /tmp/fd1
for ((i=1;i<=4;i++))
do
        echo >&3
done


for j in Raw_data/WT_*_1.fq.gz;
do
read -u3
{
  R2=${j/1_val_1/2_val_2}
  R2=${R2/1.fq.gz/2.fq.gz}
  i=${j/_val_1/}
  name=${i%_1.fq.gz}
  name=${name/Raw_data\//}

#  ~/biosoft/TrimGalore-0.6.6/trim_galore -q 20 -length 20 --fastqc -j 4  --paired $j $R2 -o 02Clean_data/ >> 02Clean_data/Trim_Galore.log 2>&1

    if [ ! -d "03Hisat2_align" ];then   mkdir -p 03Hisat2_align ; fi
    hisat2 -p 6 --dta --score-min L,5,-0.4 --sp 1,0 --mp 3,1  --summary-file 03Hisat2_align/${name}_summary.txt -x Transcriptome_Assemblly-221208/BamFile/00Hisat2_index/1-Genome_assembly_new.fasta -S 03Hisat2_align/${name}.sam  -1 02Clean_data/${name}_1_val_1.fq.gz -2 02Clean_data/${name}_2_val_2.fq.gz  2>&1 >>03Hisat2_align/error.ht2.log
    java -jar /home/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar  SortSam -I 03Hisat2_align/${name}.sam -O 03Hisat2_align/${name}_sorted.bam -SORT_ORDER coordinate > 03Hisat2_align/${name}_sort.log 2>&1 
    java -jar /home/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar  MarkDuplicates -I 03Hisat2_align/${name}_sorted.bam -O 03Hisat2_align/${name}_rmdup.bam -REMOVE_DUPLICATES true -M 03Hisat2_align/${name}_marked_dup_metrics.txt > 03Hisat2_align/${name}_rmdup.log 2>&1
     rm 03Hisat2_align/${name}_sorted.bam
     rm 03Hisat2_align/${name}.sam

    echo >&3 
   
}&
done

wait
exec 3<&-
exec 3>&-

