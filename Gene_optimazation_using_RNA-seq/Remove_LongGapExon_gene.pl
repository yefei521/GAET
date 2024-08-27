#perl Remove_LongGapExon_gene.pl SB210_DRS_DNA_Q7_maped_sorted.gtf 
use strict;
my$dir=$ENV{'PWD'};
my( %EpiOri_ne_gene, %EpiOri_check_gene, @EpiOri_ne_gene_H3K4, @H3K4me3_gene, @H3K4me3, @temp);
open(IN5,"<$ARGV[0]") or die; my@Gene=<IN5>; close IN5;
open(OUT1,">$ARGV[0]_RmLongGapExon.gtf"); 
##chr_001 AUGUSTUS        gene    283     989     .       +       .       ID=star_STRG.2.1_ORF1;Name=star_STRG.2.1_ORF1;Note=
my($chr, $start, $end, $Name, %gff, $gene_id, $mRNA_start, $mRNA_end, @m, @n, );
foreach my$gene(@Gene){
       if($gene=~/(chr_\d+)\s.+?transcript\s+(\d+)\s+(\d+)\s+\d+\s+(.)\s+.+? transcript_id "(.+?)"; cov/){ $gene_id=$5;  $gff{$gene_id}.=$gene;
       }elsif($gene=~/chr_\d+.+?exon.+?transcript_id "($gene_id)"; exon_numbe/){   $gff{$gene_id}.=$gene;
       }
    }
my$n=keys%gff; print "$n\n";
my($ori, $chr, @UTR5, @UTR3, @m, @n, @CDS_line, @lines, $y );
foreach my$name(keys%gff){ 
    @lines=split(/\n/,$gff{$name}); undef@UTR5; undef@UTR3; undef@m; undef@n; undef@CDS_line;  
    foreach my$line(@lines){ 
           if($line=~/(chr_\d+)\s.+?transcript\s+(\d+)\s+(\d+)\s+\d+\s+(.)\s+.+? transcript_id "(.+?)"; cov/){ $ori=$4; $chr=$1; $mRNA_start=$2; $mRNA_end=$3; $gene_id=$5;
           }elsif($line=~/(chr_\d+\s.+?exon)\s+(\d+)\s+(\d+)\s+.+?transcript_id "($gene_id)"; exon_numbe/){  push@UTR5,$2,$3;
           }
    }
    my$i=2; while($i<@UTR5){ if($UTR5[$i]-$UTR5[$i-1]>1000){ shift@UTR5; shift@UTR5; }else{  $i+=2; }  } @UTR3=@UTR5; 
    my$x=$y=2; while($x<@UTR3){ if($UTR3[$x]-$UTR3[$x-1]>1000){ last; }else{ $y=$x+1; } $x+=2; }  
    my$x=0; while($x<$y){ push@n,$UTR3[$x],$UTR3[$x+1];  $x+=2; }  
    print OUT1 "$chr\tAUGUSTUS\tgene\t$n[0]\t$n[-1]\t.\t$ori\t.\tID=$name;Name=$name;Note=\n";
    print OUT1 "$chr\tAUGUSTUS\tmRNA\t$n[0]\t$n[-1]\t.\t$ori\t.\tID=$name;Name=$name;Note=\n";
    $i=1; while ( $i < @n ) { print OUT1 "$chr\tAUGUSTUS\texon\t$n[$i-1]\t$n[$i]\t.\t$ori\t.\tID=$name.t1.exon;Parent=$name\n"; $i+=2; }
}
#open(IN);
