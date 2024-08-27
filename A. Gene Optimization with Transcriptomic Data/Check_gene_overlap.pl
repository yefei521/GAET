use strict;
my$dir=$ENV{'PWD'}; 
my(%gff, @GFF,%Overlap,$mark, $CDS_len, $chr, $gene_id, %gene, @Chr_gff, %temp, $gene_start, $gene_end, $n, $ID, $gene_len, $max_len, $len, $name);
open(IN1,"<$dir/7Merged_RMoverlap_gene.gff3_sorted_9Merged_RMoverlap_newgene.gff3_sorted_sorted_RmLongGapUTR_RmLong1K-UTR.gff3_OriErrorType1.gff3_sorted") or die;
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
foreach(@GFF){ 
     if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gff{$5}=$_; $ID=$5; $mark=0;  $mark++; $CDS_len=0;
     }elsif(/(chr_\d+)\s.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=($ID)/){ $gff{$ID}.=$_; $CDS_len+=$3-$2+1; $mark++;
     }elsif(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=($ID)/){ $gff{$ID}.=$_; $mark++;
     }elsif(/(chr_\d+)\s.+?UTR\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=($ID)/){ $gff{$ID}.=$_; $mark++;
     }
     if($CDS_len<150 && $mark>4){print "$gff{$ID}"; } 
}

open(IN2,"<$dir/chrname") or die;
open(OUT1,">$dir/04Overlap_gene.txt");
foreach(<IN2>){ 
    if(/^(chr_\d+)\s.+?supercontig\s+(\d+)\t(\d+)\s+/){ 
        my$chr=$1; my$contig_start=$2; my$contig_end=$3;  
        open(TXT,"<$dir/chr/$chr.gff3") or die;
        @Chr_gff=<TXT>; close TXT; 
        foreach my$line(@Chr_gff){ 
             if($line=~/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_start=$2; $gene_end=$3; $max_len=$gene_end-$gene_start+1; $name=$5; 
             %temp=%gene; delete$temp{$line};
             next if  exists$Overlap{$name};
             foreach my$m(keys%temp){
                 if($m=~/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ my$x=$2; my$y=$3; $len=$y-$x+1;
                   if( ($gene_start<$x) and ($gene_end>$y) && $len>$max_len ){   $name=$5; $Overlap{$name}++; 
                   }elsif( ($gene_start<$x) and ($gene_end<$y) and ($x<$gene_end) && ($gene_end-$x)/($gene_end-$gene_start)>0.7 && ($gene_end-$x)/($y-$x)>0.7 && $len>$max_len  ){ $name=$5; $Overlap{$name}++;
                   }elsif( ($gene_start>$x) and ($gene_end>$y) and ($y>$gene_start) && ($y-$gene_start)/($gene_end-$gene_start)>0.7 && ($y-$gene_start)/($y-$x)>0.7 && $len>$max_len ){ $name=$5; $Overlap{$name}++;
                   }elsif( ($gene_start>$x) and ($gene_end<$y) && $len>$max_len ){  $name=$5; $Overlap{$name}++; }
                 }
            } 
            print OUT1"$gff{$name}";
          }
        }
    } 
}

