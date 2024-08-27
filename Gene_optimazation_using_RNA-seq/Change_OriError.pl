#
#Check gene ori according to Specific-strand ori and select new gff to correct gene
use strict;
my$dir=$ENV{'PWD'};
open(IN1,"<$dir/SB210_DRS_DNA_Q7_maped_sorted.gtf_RmLongGapExon.gtf"); my@IN1=<IN1>; close IN1;#
open(IN2,"<","$dir/2-upd-Genome-GFF3-latest-2.gff3"); #

my(%Old_Ori, %Specific_Ori, $name, %GFF, %gene_start, %gene_end);
foreach(@IN1){ if(/^(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.+?\s+(.)\s.+?/){  $Specific_Ori{"$1"}{"$2\t$3"}=$4;   }   }
foreach(<IN2>){ if(/^(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.+?\s+(.)\s.+?ID=(.+?);/){  $Old_Ori{"$1"}{"$2\t$3"}=$4; $GFF{$_}=$_; $name=$5; $gene_start{$1}{"$2\t$3"}{"start"}=$2; $gene_end{$1}{"$2\t$3"}{"end"}=$3;
                }elsif(/^(chr_\d+)\s.+?\s+(\d+)\s+(\d+)\s+.+?\s+(.)\s.+?ID=($name);/){ $GFF{$_}.=$_; }
   }   

my(@ori, @corrd, , );
foreach my$a(@IN1){
   if($a=~/^(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?);/){ my$chr=$1; my$ori=$4; my$start=$2; my$end=$3; my$name=$5; my$a_corrd="$2\t$3";
      undef@ori;undef@corrd;
      foreach my$corrd(keys%{$Old_Ori{$chr}}){ @corrd=split(/\t/,$corrd); 
         if($corrd[0]>$end || $corrd[1]<$start){ 
         }else{ push@ori,$Old_Ori{$1}{$corrd};  print "$chr\t$corrd\t$ori\t$name\t$Old_Ori{$1}{$corrd}\n";
               # if($Old_Ori{$1}{$corrd}==$ori){ push@ori,$Old_Ori{$1}{$corrd}; }  
         }
      }
      my$mark=0; foreach my$b(@ori){ if($b ne $ori){ $mark++; }else{ $mark=0; last; } }   
      #if($mark>0){  print "$a";   }
   }
}

#open(IN4,"<$dir/");
#foreach(){}
