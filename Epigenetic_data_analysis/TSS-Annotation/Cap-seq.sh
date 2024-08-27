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


for j in 01Raw_data/*_1.fastq.gz;
do
read -u3
{
  j=${j/01Raw_data\//}
  R2=${j/1_val_1/2_val_2}
  R2=${R2/1.fastq.gz/2.fastq.gz}
  i=${j/_val_1/}
  name=${i%_1.fastq.gz}
  name=${name/02Clean_data\//}
  echo $R2
  echo $name
#  ln -s $j Rename_rawdata/${name}_1.fq.gz
#  ln -s $R2 Rename_rawdata/${name}_2.fq.gz

#  ~/biosoft/TrimGalore-0.6.6/trim_galore -q 20 -length 5 --fastqc -j 4  --paired 01Raw_data/${name}_1.fastq.gz 01Raw_data/${name}_2.fastq.gz -o 02Clean_data/ >> 02Clean_data/Trim_Galore.log 2>&1

  if [ ! -d "03Hisat2_align/$name" ];then   mkdir -p 03Hisat2_align/$name ; fi


    ~/biosoft/hisat2-2.2.1/hisat2  -p 6 --dta  --summary-file 03Hisat2_align/${name}_summary.txt --novel-splicesite-outfile 03Hisat2_align/${name}_splice-site.txt  --min-intronlen 9  -x ../index/TGD2024_Genome.fasta -S 03Hisat2_align/${name}.sam  -1 02Clean_data/${name}_1_val_1.fq.gz -2 02Clean_data/${name}_2_val_2.fq.gz 2>&1 >>03Hisat2_align/mic-error.ht2.log 
    java -jar /home/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar  SortSam -I 03Hisat2_align/${name}.sam -O 03Hisat2_align/${name}_sorted.bam -SORT_ORDER coordinate -TMP_DIR temp > 03Hisat2_align/${name}_sort.log 2>&1
    rm 03Hisat2_align/${name}.sam

    echo >&3 
   
}&
done

wait

   stringtie -p 8 -s 0.1 -o cap_merge_mapped.gtf cap_merge_mapped.bam
   
exec 3<&-
exec 3>&-

