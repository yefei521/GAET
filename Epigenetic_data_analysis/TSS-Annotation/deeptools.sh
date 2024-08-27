#for i in *bam; do { alignmentSieve -b $i --minFragmentLength 100 --maxFragmentLength 150 -o ${i}_100_150.bam; }& done

#for i in ./*100_150.bam;do { samtools index $i; }& done 

awk 'BEGIN{OFS="\t"} {j=1;  while(j<=NF){ if($3~/W/){print $1,".","SB210_CCS",$2,$2,".","+",".",".";}else{ print $1,".","SB210_CCS
",$2,$2,".","-",".","."; } j++; } }' ~/learning/S24_6mA/SQII_WT-AMT1KO/SB210_SQ2_IPD_ratio_AT_loci_penetration_withou_singletons | sort -k1.5,1.7n -k4.1,4.9n > SB210_SQ2_IPD_ratio_AT_loci_penetr
ation_withou_singletons.gff3

bedToBam -i SB210_SQ2_IPD_ratio_AT_loci_penetration_withou_singletons.gff3 -g chrome.size.txt > SB210_SQ2_IPD_ratio_AT_loci_penetration_withou_singletons.bam


for i in ./*.bam;do { bamCoverage  --outFileFormat bigwig  --skipNAs --bam $i  --outFileName ${i}.bigwig --normalizeUsing RPKM --minMappingQuality 30  --binSize 10 --ignoreDuplicates -p 6; }& done
wait
for i in ./*bigwig; do { computeMatrix scale-regions -S $i -R 2-upd-Genome-GFF3-latest.bed  -bs 10 --outFileNameMatrix ${i}_computeMatrix.txt --outFileSortedRegions ${i}_sortedRegions.bed --missingDataAsZero -b 1000 -a 1000 -o ${i}_computeMatrix.gz --regionBodyLength  1000 -p 6; }& done
wait
for i in ./*gz;do { plotHeatmap  -m $i --heatmapWidth 6  --colorMap Blues -out ${i}_Heatmap.pdf; }& done
wait

H3K4me3
{
   bigwigCompare --numberOfProcessors 12 --outFileName H3K4me3_104-160.bigwig --bigwig1 H3K4me3-input_140-160.srt.rmdup.bam.bigwig --bigwig2 H3K4me3-ChIP_140-160.srt.rmdup.bam.bigwig --skipNAs --operation second

   computeMatrix scale-regions -S H3K4me3_104-160.bigwig -R SB210_sort_sortedRegions.bed_0.85  -bs 10   -b 1000 -a 1000 -o H3K4me3_10
4-160.bigwig_computeMatrix.gz --regionBodyLength  1000 -p 12 --missingDataAsZero

   plotHeatmap  -m H3K4me3_104-160.bigwig_computeMatrix.gz --heatmapWidth 6  --colorMap Blues  -out H3K4me3_104-160.bigwig_exp-up_Heatmap.pdf --heatmapHeight 10
}


H2A.Z
{
bigwigCompare --numberOfProcessors 6 --outFileName H2A.Z_140-160.bigwig --bigwig1 H2A.Z_104-160_input.srt.rmdup.bam.bigwig --bigwig2 H2A.Z_104-160_ChIP.srt.rmdup.bam.bigwig                   

computeMatrix scale-regions -S H2A.Z_140-160.bigwig -R 2-upd-Genome-GFF3-latest-2.gff3_gene_exp0.bed  -bs 10   -b 1000 -a 1000 -o H2A.Z_104-160.bigwig_exp0_computeMatrix.gz --regionBodyLength  1000 -p 12 --missingDataAsZero 

plotHeatmap  -m H2A.Z_104-160.bigwig_exp0_computeMatrix.gz --heatmapWidth 6 --heatmapHeight 10  --colorMap Blues -out H2A.Z_104-160.bigwig_exp0_Heatmap.pdf --zMax 20 --yMax 6


}

#for i in *input_rmdup.bam_100_150.bam.bigwig;do { j=${i/input/IP}; n=${i%input*}; bigwigCompare --numberOfProcessors 6 --outFileName ${n}_log2ratio.bw --bigwig1 $i --bigwig2 $j ; }& done

#for i in ./*bw; do { computeMatrix scale-regions -S $i -R 2-upd-Genome-GFF3-latest.bed  -bs 10 --outFileNameMatrix ${i}_computeMatrix.txt --outFileSortedRegions ${i}_sortedRegions.bed --missingDataAsZero -b 1000 -a 1000 -o ${i}_computeMatrix.gz --regionBodyLength  1000 -p 6; }& done

#for i in ./*bw_computeMatrix.gz;do { plotHeatmap  -m $i --heatmapWidth 6  --colorMap Blues --yMax 1.5  -out ${i}_Heatmap.pdf; }& done



#computeMatrix scale-regions -S veg_6mA.sorted.bigwig H3K4me3-ChIP_140-160_log2ratio.bw H2A.Z_104-160_ChIP_log2ratio.bw H3K4me3-ChIP_140-160.srt.rmdup.bam.bigwig  -R 2-upd-Genome-GFF3-latest-2.gff3_gene_exp-up  -bs 10 --outFileNameMatrix exp-up-gene_computeMatrix.txt --outFileSortedRegions exp-up-gene_sortedRegions.bed --missingDataAsZero -b 1000 -a 1000 -o exp-up-gene_computeMatrix.gz --regionBodyLength  1000 -p 6 &



#multiBigwigSummary bins -b SB210_2_IP_rmdup.bam_100_150.bam.bigwig SB210_2_input_rmdup.bam_100_150.bam.bigwig SB210_3_IP_rmdup.bam_100_150.bam.bigwig SB210_3_input_rmdup.bam_100_150.bam.bigwig SB210_4_IP_rmdup.bam_100_150.bam.bigwig SB210_4_input_rmdup.bam_100_150.bam.bigwig A3KO_1_IP_rmdup.bam_100_150.bam.bigwig A3KO_1_input_rmdup.bam_100_150.bam.bigwig A3KO_2_IP_rmdup.bam_100_150.bam.bigwig A3KO_2_input_rmdup.bam_100_150.bam.bigwig A3KO_3_IP_rmdup.bam_100_150.bam.bigwig A3KO_3_input_rmdup.bam_100_150.bam.bigwig A4KO_1_IP_rmdup.bam_100_150.bam.bigwig A4KO_1_input_rmdup.bam_100_150.bam.bigwig A4KO_2_IP_rmdup.bam_100_150.bam.bigwig A4KO_2_input_rmdup.bam_100_150.bam.bigwig A4KO_3_IP_rmdup.bam_100_150.bam.bigwig A4KO_3_input_rmdup.bam_100_150.bam.bigwig A34KO_1_IP_rmdup.bam_100_150.bam.bigwig A34KO_1_input_rmdup.bam_100_150.bam.bigwig A34KO_2_IP_rmdup.bam_100_150.bam.bigwig A34KO_2_input_rmdup.bam_100_150.bam.bigwig A34KO_3_IP_rmdup.bam_100_150.bam.bigwig A34KO_3_input_rmdup.bam_100_150.bam.bigwig -o 24sample_bigwig.npz

#plotCorrelation --corData 24sample_bigwig.npz --corMethod spearman --whatToPlot heatmap --plotFile 24sample_bigwig.svg --plotTitle "AMT34KO MeRIP-seq"  --colorMap PuBuGn


