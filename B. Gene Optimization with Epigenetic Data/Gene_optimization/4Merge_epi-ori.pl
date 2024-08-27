#perl 4Merge_epi-ori.pl veg_6mA.sorted.gff_m6A_ori H2A.Z_group.gff3_H2A_ori H3K4me3_group.gff3_H3K4me3_ori
use strict;
my(%m6A, %H2A, %H3K4me3, %gene_no_6mA, %gene_no_H2A, %gene_no_H3K4, %gene_mid_6mA, %gene_mid_H2A, %gene_mid_H3K4, %gene_ne_6mA, %gene_ne_H2A, %gene_ne_H3K4);
my$dir=$ENV{'PWD'};
open(IN1,"<$dir/chrname");
open(IN2,"<$ARGV[0]"); #6mA
open(IN3,"<$ARGV[1]"); #H2A.Z
open(IN4,"<$ARGV[2]"); #H3K4me3
my@m6A=<IN2>; close IN2;
my@H2A=<IN3>; close IN3;
my@H3K4me3=<IN4>; close IN4;

#star_STRG.2.1_ORF1      +       -       0.577431539187913
#star_STRG.5.1_ORF1      +       +       0.48465855050409
foreach my$m6A(@m6A){  if($m6A=~/(.+?)\s+(.)\s+(.)\s+(.+?)\s+(.+?)\s+(.+?)\s+(\d+)\n/){  my$ratio=(int(100*$6))/100;  $m6A{"$1\t$2"}="m6A:$3($ratio)"; 
                                                                     if($6==0){ $gene_no_6mA{"$1\t$2"}++; }elsif($6<0.6 && $6>0.4){$gene_mid_6mA{"$1\t$2"}++; }
                                                                     if( ($2 ne $3)&&($6>0) ){$gene_ne_6mA{"$1\t$2"}++;  } 
                       }elsif($m6A=~/(.+?)\s+(.)\s+(.)\s+(0)\s+(\d+)\n/){ $gene_no_6mA{"$1\t$2"}++;}
}
foreach my$H2A(@H2A){  if($H2A=~/(.+?)\s+(.)\s+(.)\s+(.+?)\s+(.+?)\s+(.+?)\s+(\d+)\n/){  my$ratio=(int(100*$6))/100;  $H2A{"$1\t$2"}="H2A:$3($ratio)"; 
                                                                     if($6==0){ $gene_no_H2A{"$1\t$2"}++; }elsif($6<0.6 && $6>0.4){$gene_mid_H2A{"$1\t$2"}++; } 
                                                                     if( ($2 ne $3)&&($6>0) ){$gene_ne_H2A{"$1\t$2"}++;  } 
                       }elsif($H2A=~/(.+?)\s+(.)\s+(.)\s+(0)\s+(\d+)\n/){$gene_no_6mA{"$1\t$2"}++;}  
}
foreach my$H3K4me3(@H3K4me3){  if($H3K4me3=~/(.+?)\s+(.)\s+(.)\s+(.+?)\s+(.+?)\s+(.+?)\s+(\d+)\n/){  my$ratio=(int(100*$6))/100;  $H3K4me3{"$1\t$2"}="H3K4me3:$3($ratio)"; 
                                                                     if($6==0){$gene_no_H3K4{"$1\t$2"}++;}elsif($6<0.6 && $6>0.4){$gene_mid_H3K4{"$1\t$2"}++; }
                                                                     if( ($2 ne $3)&&($6>0) ){$gene_ne_H3K4{"$1\t$2"}++;  } 
                         }elsif($H3K4me3=~/(.+?)\s+(.)\s+(.)\s+(0)\s+(\d+)\n/){ $gene_no_6mA{"$1\t$2"}++;  } 
}

open(OUT1,">$dir/4Merge_epi-ori.gff");
open(OUT2,">$dir/4Merge_gene-no-EpiOri.gff");
open(OUT3,">$dir/4Merge_gene-mid-EpiOri.gff");
open(OUT4,">$dir/4Merge_gene-ne-EpiOri.gff");
foreach my$n(keys%m6A){ 
     print OUT1"$n\t$m6A{$n}\t$H2A{$n}\t$H3K4me3{$n}\n";
     if(exists$gene_no_6mA{$n} && exists$gene_no_H2A{$n} && exists$gene_no_H3K4{$n}){print OUT2"$n\t$m6A{$n}\t$H2A{$n}\t$H3K4me3{$n}\n";}
     if(exists$gene_mid_6mA{$n} && exists$gene_mid_H2A{$n} && exists$gene_mid_H3K4{$n}){print OUT3"$n\t$m6A{$n}\t$H2A{$n}\t$H3K4me3{$n}\n";}
     if(exists$gene_ne_6mA{$n} && exists$gene_ne_H2A{$n} && exists$gene_ne_H3K4{$n}){print OUT4 "$n\t$m6A{$n}\t$H2A{$n}\t$H3K4me3{$n}\n"; }
}
