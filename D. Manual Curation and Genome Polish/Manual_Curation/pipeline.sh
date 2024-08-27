 delete this gene: T00161861

perl 01Origanized_multiGFF3_groups.pl > merge.gtf

gffread merge.gtf -g ../../TGD2024_Genome.fasta -w merge.fasta

mv merge.* Manual-Check_RawGFF/

perl 02Split_CDS-mRNA-AS-NAT.pl

#perl 03ORFfind.pl Split_CDS-mRNA-AS-NAT/mRNA_need_coding.fast
 predict orf using orffinder with the one line fasta format

 Get the *fasta_orf1 file

perl 03Select_ATG2TGA_cds.pl Split_CDS-mRNA-AS-NAT/mRNA_need_coding.fasta_orf1

mv Split_CDS-mRNA-AS-NAT/mRNA_need_coding.fasta_orf* Predict_mRNA2CDS_GFF/

-------------------------
#check mrna number

#less Predict_mRNA2CDS_GFF/mRNA_need_coding.fasta_orf1_ATG2TGA_longest | grep '>' | cut -f1 > temp.ATG2TGA_name.txt

#grep '>' Split_CDS-mRNA-AS-NAT/mRNA_need_coding.fasta > temp.needCoding_name.txt

#cat temp.needCoding_name.txt temp.ATG2TGA_name.txt | sort | uniq -u | awk '{match($1,/>(.+)/,a); print a[1]}' > temp_uniq.txt


#awk 'NR==FNR{a[$1]++;} NR>FNR && $3~/mRNA/{match($9,/ID=(.+);Man/,n); if(n[1] in a){print $0}} NR>FNR && $3~/exon/{match($9,/Parent=(.+);/,n); if(n[1] in a){print $0;}}' temp_uniq.txt Manual-Check_RawGFF/merge.gtf > temp_uniq_need_check.gtf
--------------------


------------
#
perl 04Merge_orf1_to_UTR-gff.pl Split_CDS-mRNA-AS-NAT/mRNA_need_coding.gtf Predict_mRNA2CDS_GFF/mRNA_need_coding.fasta_orf1_ATG2TGA_longest


gffread Predict_mRNA2CDS_GFF/mRNA_need_coding.fasta_orf1_ATG2TGA_longest_UTR.gff3 -g ../../TGD2024_Genome.fasta -x mRNA_need_coding.fasta_orf1_ATG2TGA_longest_UTR.gff3_CDS.fasta

##Compare the CDS length in UTR.gff and nedd_coding.gtf 
awk 'NR==FNR && $1~/>/{match($1,/>(.+)/,a); name=a[1];} NR==FNR && $1!~/>/{L[name]=L[name]""$1;} NR>FNR && $1~/>/{match($1,/>(.+).t1/,a); name=a[1];} NR>FNR && $1!~/>/{R[name]=R[name]""$1;} END{for(i in R){lenL=length(L[i]); lenR=length(R[i]); if(lenL!=lenR){print i;}}}' Predict_mRNA2CDS_GFF/mRNA_need_coding.fasta_orf1_ATG2TGA_longest mRNA_need_coding.fasta_orf1_ATG2TGA_longest_UTR.gff3_CDS.fasta | wc -l

##Select the gene "not ATG2TGA" in UTR.gff3
less mRNA_need_coding.fasta_orf1_ATG2TGA_longest_UTR.gff3_CDS.fasta |awk '$1~/>/{name=$1;} $1!~/>/{a[name]=a[name]""$1} END{for(i in a){ if(a[i]~/^ATG.+TGA$/){print i; print a[i];}else{print i,"error";}}}' | grep 'error'  | awk '{match($1,/>(.+).t1/,a) ; print a[1];}'  > Predict_mRNA2CDS_GFF/Not_ATG2TGA-in-UTR-GFF.txt

awk 'NR==FNR && $1~/>/{match($1,/>(.+)/,a); name=a[1];} NR==FNR && $1!~/>/{L[name]=L[name]""$1;} NR>FNR{N[$1]++;} END{for(i in N){print ">"i; print L[i];} }' Predict_mRNA2CDS_GFF/mRNA_need_coding.fasta_orf1_ATG2TGA_longest Predict_mRNA2CDS_GFF/Not_ATG2TGA-in-UTR-GFF.txt  > Predict_mRNA2CDS_GFF/Not_ATG2TGA-in-UTR-GFF_CDS.fasta

###wgsim -1 70 -2 70 -r 0 -R 0 -X 0 -S 1 -e 0 -d 0 ./Predict_mRNA2CDS_GFF/Not_ATG2TGA-in-UTR-GFF_CDS.fasta ./Predict_mRNA2CDS_GFF/Not_ATG2TGA-in-UTR-GFF_CDS_sim_1.fq ./Predict_mRNA2CDS_GFF/Not_ATG2TGA-in-UTR-GFF_CDS_sim_2.fq


##awk 'NR==FNR && $1~/>/{match($1,/>(.+)/,a); name=a[1];} NR==FNR && $1!~/>/{L[name]=L[name]""$1;} NR>FNR && $1~/>/{match($1,/>(.+).t1/,a); name=a[1];} NR>FNR && $1!~/>/{R[name]=R[name]""$1;} END{for(i in R){lenL=length(L[i]); lenR=length(R[i]); if(lenL!=lenR){print i,lenL,lenR;}}}' Predict_mRNA2CDS_GFF/mRNA_need_coding.fasta_orf1_ATG2TGA_longest mRNA_need_coding.fasta_orf1_ATG2TGA_longest_UTR.gff3_CDS.fasta

Find 33 gene has not correct ATG or TGA due to the only 1 base in another exon;
#awk 'NR==FNR{a[$1]++;} NR>FNR && $3~/gene/{match($9,/ID=(.+);Name/,b); if(b[1] in a){print}}  NR>FNR && $3~/mRNA/{match($9,/ID=(.+).t1/,b); if(b[1] in a){print}}  NR>FNR && $9~/Parent/{match($9,/Parent=(.+).t1/,b); if(b[1] in a){print}}' Predict_mRNA2CDS_GFF/Not_ATG2TGA-in-UTR-GFF.txt  Predict_mRNA2CDS_GFF/mRNA_need_coding.fasta_orf1_ATG2TGA_longest_UTR.gff3_checked  >  Predict_mRNA2CDS_GFF/Not_ATG2TGA-in-UTR-GFF.txt_checked.gff3

Rename gene
###

mkdir Rename_gene && cp Predict_mRNA2CDS_GFF/mRNA_need_coding.fasta_orf1_ATG2TGA_longest_UTR.gff3_checked Rename_gene/TGD2024_need_rename_gene.gff3

gffread TGD2024_need_rename_gene.gff3 -g ../../../TGD2024_Genome.fasta -w TGD2024_need_rename_gene.fasta
~/biosoft/ncbi-blast-2.10.0+/bin/blastn -db 3-upd-cds-fasta-2021.fasta -query TGD2024_need_rename_gene.fasta -outfmt 7 -out TGD2024blastTGD2021_otf7 -num_threads 20 &

perl 06Rename_gene.pl

cp Rename_gene/TGD2024_duplicated_nedd_check.txt Rename_gene/TGD2024_duplicated_nedd_check.txt_checked 

perl 06Check_duplicated-gene.pl

less -S Interpro_anno/*tsv | grep -v 'Coil' | grep -v 'MobiDBLite' |  grep -v 'TMHMM' | grep -v 'Phobius' | cut -f1,4,13  | awk '$3!~/^-/' | uniq  > Interpro_anno/TGD2024_need_rename_gene_pep.fasta_annotation.txt

perl 06Rename_gene.pl

gffread Rename_gene/TGD2024_renamed.gff3 -g ../../TGD2024_Genome.fasta -w Rename_gene/TGD2024_renamed.gff3_gene.fasta
gffread Rename_gene/TGD2024_renamed.gff3 -g ../../TGD2024_Genome.fasta -x Rename_gene/TGD2024_renamed.gff3_CDS.fasta
awk 'NR==FNR{a[$1]++} NR>FNR && $1~/>/{match($1,/>(.+)/,b); name=b[1]} NR>FNR && $1!~/>/ && name in a{print ">"name; print $0}' Rename_gene/LowQuality_Tname2gname.txt  Rename_gene/3-upd-cds-fasta-2021.fasta > LowQuality_CDS.fasta   
awk 'NR==FNR{a[$1]++} NR>FNR && $1~/>/{match($1,/>(.+)/,b); name=b[1]} NR>FNR && $1!~/>/ && name in a{print ">"name; print $0}' Rename_gene/LowQuality_Tname2gname.txt  Rename_gene/4-upd-Protein-fasta-2021.fasta > LowQuality_pep.fasta




less Rename_gene/TGD2024_renamed.gff3_gene.fasta | awk '$1~/>/{name=$1;} $1!~/>/{seq[name]=seq[name]""$1;} END{for(i in seq){if(i~/t1/){print i; print seq[i];}}}' > Rename_gene/TGD2024_renamed.gff3_mRNA.fasta

wk 'NR==FNR && $1~/>/{match($1,/>(.+)/,a); name=a[1];} NR==FNR && $1!~/>/{seq[name]=seq[name]""$0;}  NR>FNR{if($4~/Unknow/){}else{match($4,/TGDAnno=(.+)/,a);Anno1[$1]=a[1];} if($5!~/Anno=$/){match($5,/Anno=(.+)/,a);Anno2[$1]=a[1];}  Tn[$1]=$2;} END{for(i in seq){print ">"i"|"Tn[i],Anno1[i],Anno2[i]; print seq[i];}}' Rename_gene/TGD2024_renamed.gff3_pep.fasta TGD2024_Tname2Annotation.xlxs > Rename_gene/TGD2024_renamed.gff3_pep_addAnno.fasta









































