#perl 3Scan_gene_H3K4me3.pl H3K4me3_group.gff3 7Merged_RMoverlap_gene.gff3_sorted_9Merged_RMoverlap_newgene.gff3_sorted_sorted
use strict;
my(@H3K4me3, @gene, @H3K4me3_gene, @temp, @sort_temp);
my$dir=$ENV{'PWD'};
open(IN1,"<$dir/chrname") or die;
open(IN2,"<$ARGV[0]") or die; @H3K4me3=<IN2>; close IN2;
open(IN3,"<$ARGV[1]") or die; @gene=<IN3>; close IN3;
open(OUT1,">$ARGV[0]_H3K4me3_ori") or die;
foreach my$chr(<IN1>){chomp $chr; print"$chr\n";
   foreach my$gene(@gene){ if($gene=~/$chr.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);/){ my$gene_start=$1; my$gene_end=$2; my$gene_name=$4; undef@H3K4me3_gene; my$ori=$3;
       foreach my$H3K4me3(@H3K4me3){  if($H3K4me3=~/$chr.+?H3K4me3\s+(\d+)\s+(\d+)\s/){ if($1>=$gene_start && $1<=$gene_end){ push@H3K4me3_gene,$1; } } }
       my$sum_ratio=0; my$ave_ratio=0 ; my$n=@H3K4me3_gene; undef@temp; my$gene_len=$gene_end-$gene_start; 
       foreach my$H3K4me3(@H3K4me3_gene){ $sum_ratio+=($H3K4me3-$gene_start)/($gene_end-$gene_start); my$ratio=($H3K4me3-$gene_start)/($gene_end-$gene_start); push@temp,$ratio; }
       if(@H3K4me3_gene>0){ $ave_ratio=$sum_ratio/@H3K4me3_gene; }else{ $ave_ratio=0; } @sort_temp=sort{$a<=>$b}@temp;
       if($sort_temp[0]<0.5){ print OUT1"$gene_name\t$ori\t+\t$sort_temp[0]\t$sort_temp[-1]\t$ave_ratio\t$gene_len\n";}else{print OUT1"$gene_name\t$ori\t-\t$sort_temp[0]\t$sort_temp[-1]\t$ave_ratio\t$gene_len\n"; }
    }  
  }
}
