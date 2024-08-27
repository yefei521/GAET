use strict;
my($before, $after, @GFF, $chr, $gene_id, %gene, @Chr_gff, %temp, $gene_start, $gene_end, $n, %gff, @CDS, );
my$dir=$ENV{'PWD'};
open(IN1,"<$dir/03T-t-manualcheck_220224_rmErrorGene.gff3");
open(GFF_OUT,">$dir/temp.gff3");
my$cmd=`rm chr/chr*`;
print "$cmd";

while(<IN1>){
       if(/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+.+?ID=(chr_\d+);/){ 
       }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5;  $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?CDS.+?Parent=$gene_id/){   $gff{$gene_id}.=$_; 
       }elsif(/chr_\d+.+?MANUAL_ENTRY.+?ID=(.+?);Name/){ 
       }elsif(/chr_\d+.+?rDNA.+?ID=(.+?);Name/){  
       }
    }
foreach(keys%gff){
       my@temp=split(/\n/,$gff{$_}); undef@CDS;
       foreach my$m(@temp){ 
           if($m=~/((chr_\d+)\s.+?gene)\s+(\d+)\s+(\d+)\s+(.\s+(.)\s+.+?ID=(.+?);Name=(.+?);.+)/){
                  $before=$1; $after=$5;
           }elsif($m=~/(chr_\d+)\s.+?CDS\s+(\d+)\s+(\d+)\s+/){ push@CDS,$2,$3; }
       }
       print GFF_OUT"$before\t$CDS[0]\t$CDS[-1]\t$after\n";
}
open(GFF_IN,"<$dir/temp.gff3");
@GFF=<GFF_IN>; close GFF_IN;
foreach(@GFF){
       if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){
           $gene_id=$5; $chr=$1;
           open(TXT,">>$dir/chr/$chr.gff3") or die;
           print TXT"$_"; close TXT;
           $gene{$_}=$_;
       }
    }

open(IN2,"<$dir/chrname") or die;
open(OUT1,">$dir/04Overlap_gene_CDS.txt");
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

my$cmd=`rm temp.gff3`;
print "$cmd";
