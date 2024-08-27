
# bwa index 1-Genome_assembly_new.fasta

bwa mem -t 12 -M -R '@RG\tID:SB210\tSM:SB210-2\tLB:WES\tPL:Illumina' 1-Genome_assembly_new.fasta  02Clean_data/SB210_1_1_val_1.fq.gz 02Clean_data/SB210_1_2_val_2.fq.gz | samtools sort -@ 12 -o 03BWA_align/SB210_1.bam
bwa mem -t 12 -M -R '@RG\tID:SB210\tSM:SB21--1\tLB:WES\tPL:Illumina' 1-Genome_assembly_new.fasta  02Clean_data/SB210_2_1_val_1.fq.gz 02Clean_data/SB210_2_2_val_2.fq.gz | samtools sort -@ 12 -o 03BWA_align/SB210_2.bam

java -jar ~/biosoft/picard-tools-1.119/SortSam.jar I=03BWA_align/SB210_1.bam O=03BWA_align/SB210_1_sorted.bam SORT_ORDER=coordinate
java -jar ~/biosoft/picard-tools-1.119/SortSam.jar I=03BWA_align/SB210_2.bam O=03BWA_align/SB210_2_sorted.bam SORT_ORDER=coordinate


java -jar ~/biosoft/picard-tools-1.119/MarkDuplicates.jar I=03BWA_align/SB210_2_sorted.bam O=03BWA_align/SB210_2_rmdup.bam REMOVE_DUPLICATES=true  METRICS_FILE=03BWA_align/log
java -jar ~/biosoft/picard-tools-1.119/MarkDuplicates.jar I=03BWA_align/SB210_1_sorted.bam O=03BWA_align/SB210_1_rmdup.bam REMOVE_DUPLICATES=true  METRICS_FILE=03BWA_align/log

samtools flagstat 03BWA_align/SB210_1_rmdup.bam > 03BWA_align/SB210_1_rmdup_deduped_stat.txt
samtools flagstat 03BWA_align/SB210_2_rmdup.bam > 03BWA_align/SB210_2_rmdup_deduped_stat.txt

samtools depth 03BWA_align/SB210_1_rmdup.bam > 03BWA_align/SB210_1_rmdup_deduped_coverage.txt &
samtools depth 03BWA_align/SB210_2_rmdup.bam > 03BWA_align/SB210_2_rmdup_deduped_coverage.txt

samtools view -@ 10 -q 30 03BWA_align/SB210_1_rmdup.bam > 03BWA_align/SB210_1_rmdup_filter.bam


#conda activate GATK

#gatk CreateSequenceDictionary -R 1-Genome_assembly_new.fasta

gatk HaplotypeCaller -R 1-Genome_assembly_new.fasta --emit-ref-confidence GVCF -I 03BWA_align/SB210_2_rmdup.bam -O 03BWA_align/SB210_2_rmdup.gvcf
gatk HaplotypeCaller -R 1-Genome_assembly_new.fasta --emit-ref-confidence GVCF -I 03BWA_align/SB210_1_rmdup.bam -O 03BWA_align/SB210_1_rmdup.gvcf

