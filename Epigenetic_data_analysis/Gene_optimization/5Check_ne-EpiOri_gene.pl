use strict;
my$dir=$ENV{'PWD'};
my( %EpiOri_ne_gene, %EpiOri_check_gene, @EpiOri_ne_gene_H3K4, @H3K4me3_gene, @H3K4me3, @temp);
open(IN1,"<$dir/4Merge_gene-ne-EpiOri.gff") or die; my@EpiOri_ne_gene=<IN1>; close IN1;
open(IN2,"<$dir/veg_6mA.sorted.gff") or die; my@m6A_gff=<IN2>; close IN2;
open(IN3,"<$dir/H2A.Z_group.gff3") or die; my@H2A_gff=<IN3>; close IN3;
open(IN4,"<$dir/H3K4me3_group.gff3") or die; my@H3K4_gff=<IN4>; close IN4;
open(IN5,"<$dir/7Merged_RMoverlap_gene.gff3_sorted_9Merged_RMoverlap_newgene.gff3_sorted_sorted") or die; my@Gene=<IN5>; close IN5;
#open(IN5,"<$dir/5Checked.gff3") or die; my@Gene=<IN5>; close IN5;
open(OUT1,">$dir/5Checked.gff3"); 
##chr_001 AUGUSTUS        gene    283     989     .       +       .       ID=star_STRG.2.1_ORF1;Name=star_STRG.2.1_ORF1;Note=
my($chr, $start, $end, $Name, %gff, $gene_id, $mRNA_start, $mRNA_end, @m, @n, );
foreach my$gene(@Gene){
       if($gene=~/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+.+?ID=(chr_\d+);/){
       }elsif($gene=~/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5;  $gff{$gene_id}.=$gene;
       }elsif($gene=~/chr_\d+.+?mRNA.+?Parent=$gene_id/){   $gff{$gene_id}.=$gene;
       }elsif($gene=~/chr_\d+.+?UTR.+?Parent=$gene_id/){   $gff{$gene_id}.=$gene;
       }elsif($gene=~/chr_\d+.+?CDS.+?Parent=$gene_id/){   $gff{$gene_id}.=$gene;
       }elsif($gene=~/chr_\d+.+?MANUAL_ENTRY.+?ID=(.+?);Name/){
       }elsif($gene=~/chr_\d+.+?rDNA.+?ID=(.+?);Name/){
       }
    }
my$n=keys%gff; print "$n\n";
my($ori, $chr, @UTR5, @UTR3, @m, @n, @CDS_line, @lines, $y );
foreach my$name(keys%gff){ 
    @lines=split(/\n/,$gff{$name}); undef@UTR5; undef@UTR3; undef@m; undef@n; undef@CDS_line; 
    foreach my$line(@lines){ 
           if($line=~/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=.+?;Name=.+?;.+/){ $ori=$4; $chr=$1;
           }elsif($line=~/(chr_\d+\s.+?mRNA)\s+(\d+)\s+(\d+)\s+(.+)/){ $mRNA_start=$2; $mRNA_end=$3; 
           }elsif($line=~/(chr_\d+\s.+?five_prime_UTR)\s+(\d+)\s+(\d+)\s+(.+)/){  push@UTR5,$2,$3;
           }elsif($line=~/(chr_\d+\s.+?three_prime_UTR)\s+(\d+)\s+(\d+)\s+(.+)/){ push@UTR3,$2,$3;
           }elsif($line=~/(chr_\d+\s.+?CDS)\s+(\d+)\s+(\d+)\s+(.+)/){ push@CDS_line,$line;
           }
    }
           if(@UTR5<3 && @UTR3<3){ print OUT1"$gff{$name}"; next; }
           if($ori eq "+"){
               my$i=2; while($i<@UTR5){ if($UTR5[$i]-$UTR5[$i-1]>1000){ shift@UTR5; shift@UTR5; }else{  $i+=2; }  } @m=@UTR5;  push@m,$mRNA_start if@m==0;
               my$x=$y=2; while($x<@UTR3){ if($UTR3[$x]-$UTR3[$x-1]>1000){ last; }else{ $y=$x+1; } $x+=2; }  
               if(@UTR3==0){push@n,$mRNA_end; }else{ my$x=0; while($x<$y){ push@n,$UTR3[$x],$UTR3[$x+1];  $x+=2; } } 
               print OUT1 "$chr\tAUGUSTUS\tgene\t$m[0]\t$n[-1]\t.\t$ori\t.\tID=$name;Name=$name;Note=\n";
               print OUT1 "$chr\tAUGUSTUS\tmRNA\t$m[0]\t$n[-1]\t.\t$ori\t.\tID=$name;Name=$name;Note=\n";
               $i=1; while ( $i < @m ) { print OUT1 "$chr\tAUGUSTUS\tfive_prime_UTR\t$m[$i-1]\t$m[$i]\t.\t$ori\t.\tID=$name.t1.utr5;Parent=$name\n"; $i+=2; }
               $i=0; while ($i<@CDS_line) {   print OUT1 "$CDS_line[$i]\n"; $i++; }
               $i=1; while ( $i < @n ) { print OUT1 "$chr\tAUGUSTUS\tthree_prime_UTR\t$n[$i-1]\t$n[$i]\t.\t$ori\t.\tID=$name.t1.utr3;Parent=$name\n"; $i+=2; }
           }elsif($ori eq "-"){
               my$i=$y=2; while($i<@UTR5){ if($UTR5[$i]-$UTR5[$i-1]>1000){ last; }else{ $y=$i+2; } $i+=2; } 
               if(@UTR5==0){push@m,$mRNA_end }else{my$x=1; while($x<$y){ push@m,$UTR5[$x-1],$UTR5[$x];  $x+=2; } } # push@m,$mRNA_end if@m==0;
               my$x=2; while($x<@UTR3){ if($UTR3[$x]-$UTR3[$x-1]>1000){ shift@UTR3; shift@UTR3; }else{  $x+=2; }  } @n=@UTR3; push@n,$mRNA_start if@n==0;
               print OUT1 "$chr\tAUGUSTUS\tgene\t$n[0]\t$m[-1]\t.\t$ori\t.\tID=$name;Name=$name;Note=\n";
               print OUT1 "$chr\tAUGUSTUS\tmRNA\t$n[0]\t$m[-1]\t.\t$ori\t.\tID=$name;Name=$name;Note=\n";
               $i=1; while ( $i < @n ){ print OUT1 "$chr\tAUGUSTUS\tthree_prime_UTR\t$n[$i-1]\t$n[$i]\t.\t$ori\t.\tID=$name.t1.utr3;Parent=$name\n"; $i+=2; }
               $i=0; while ($i<@CDS_line){   print OUT1 "$CDS_line[$i]\n"; $i++; }
               $i=1; while ( $i < @m ){ print OUT1 "$chr\tAUGUSTUS\tfive_prime_UTR\t$m[$i-1]\t$m[$i]\t.\t$ori\t.\tID=$name.t1.utr5;Parent=$name\n"; $i+=2; }
           }
}
