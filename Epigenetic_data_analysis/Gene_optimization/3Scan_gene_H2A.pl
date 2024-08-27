#perl 3Scan_gene_H2A.pl H2A.Z_group.gff3 7Merged_RMoverlap_gene.gff3_sorted_9Merged_RMoverlap_newgene.gff3_sorted_sorted
use strict;
my(@H2A, @gene, @H2A_gene, @temp, @sort_temp);
my$dir=$ENV{'PWD'};
open(IN1,"<$dir/chrname") or die;
open(IN2,"<$ARGV[0]") or die; @H2A=<IN2>; close IN2;
open(IN3,"<$ARGV[1]") or die; @gene=<IN3>; close IN3;
open(OUT1,">$ARGV[0]_H2A_ori") or die;
foreach my$chr(<IN1>){chomp $chr; print"$chr\n";
   foreach my$gene(@gene){ if($gene=~/$chr.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);/){ my$gene_start=$1; my$gene_end=$2; my$gene_name=$4; undef@H2A_gene; my$ori=$3;
       foreach my$H2A(@H2A){  if($H2A=~/$chr.+?H2A.Z\s+(\d+)\s+(\d+)\s/){ if($1>=$gene_start && $1<=$gene_end){ push@H2A_gene,$1; } } }
       my$sum_ratio=0; my$ave_ratio=0 ; my$n=@H2A_gene;  undef@temp; my$gene_len=$gene_end-$gene_start+1;
       foreach my$H2A(@H2A_gene){ $sum_ratio+=($H2A-$gene_start)/($gene_end-$gene_start); my$ratio=($H2A-$gene_start)/($gene_end-$gene_start); push@temp,$ratio;  }
       if(@H2A_gene>0){ $ave_ratio=$sum_ratio/@H2A_gene; }else{ $ave_ratio=0; } @sort_temp=sort{$a<=>$b}@temp;
       if($ave_ratio<0.5){ print OUT1"$gene_name\t$ori\t+\t$sort_temp[0]\t$sort_temp[-1]\t$ave_ratio\t$gene_len\n";}else{print OUT1"$gene_name\t$ori\t-\t$sort_temp[0]\t$sort_temp[-1]\t$ave_ratio\t$gene_len\n"; }
    }  
  }
}
