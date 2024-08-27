use strict;
my$dir=$ENV{'PWD'}; 
my(@GFF, $chr, $gene_id, %gene, @Chr_gff, %temp, $gene_start, $gene_end, $n);
open(IN1,"<$dir/03T-t-manualcheck_220224_rmErrorGene.gff3") or die;
#open(IN1,"<$dir/temp.gff");
@GFF=<IN1>; close IN1;
#chr_001 AUGUSTUS        supercontig     1       1459056 .       .       .       ID=chr_001;Name=chr_001
my$cmd=`rm chr/chr*`;
print "$cmd";
foreach(@GFF){
       if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ 
           $gene_id=$5; $chr=$1;
           open(TXT,">>$dir/chr/$chr.gff3") or die; 
           print TXT"$_"; close TXT; 
           $gene{$_}=$_; 
       }
    }

open(IN2,"<$dir/chrname") or die;
open(OUT1,">$dir/04Overlap_gene.txt");
foreach(<IN2>){ 
    if(/^(chr_\d+)\s.+?supercontig\s+(\d+)\t(\d+)\s+/){ 
        my$chr=$1; my$contig_start=$2; my$contig_end=$3;  
        open(TXT,"<$dir/chr/$chr.gff3") or die;
        @Chr_gff=<TXT>; close TXT; 
        foreach my$line(@Chr_gff){ 
             if($line=~/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_start=$2; $gene_end=$3;  }else{ next; }
             %temp=%gene; delete$temp{$line};
             foreach my$m(keys%temp){
                 if($m=~/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ my$x=$2; my$y=$3;
                   if( ($gene_start<$x) and ($gene_end>$y)  ){ print OUT1"Error\n$line$m"; 
                   }elsif( ($gene_start<$x) and ($gene_end<$y) and ($x<$gene_end)   ){ print OUT1"Error\n$line$m"; 
                   }elsif( ($gene_start>$x) and ($gene_end>$y) and ($y>$gene_start) ){ print OUT1"Error\n$line$m"; 
                   }elsif( ($gene_start>$x) and ($gene_end<$y) ){ print OUT1"Error\n$line$m";   }
                 }
             }
        }
    } 
}

