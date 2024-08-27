#####################################################
#README before use this script
#   1. mkdir chr 
#   2. mkdir chr_filter
#   3.Please remove files under directory named chr & chr_filter
#   4.Change input-file IN1 and remove # in line 13-19 
#   5.Change input gff file IN4 and IN5
####################################################
use strict;
my$dir=$ENV{'PWD'};
my(@chr, @gff, @m6A_corrd, $j, $n, $number, $gene, $chr, $mid, $m, $p, $q, $x, @m6A, %m6A_corrd, %gene);
open(IN1,"<$dir/veg_6mA.sorted.gff");
#chr_001 .       SB210_veg_6mA   689     689     .       +       .       .
foreach my$line(<IN1>){
   if($line=~/^(chr_\d+)\t.\t.+?\t(\d+)\t(\d+)\t/){
       open(TXT,">>$dir/chr/$1.txt");
       print TXT "$line";
       close TXT;
   }
}
print "STEP1: Split m6A into group and import into m6A_gene.txt\nPlease wait ........\n";
@chr=glob"$dir/chr_6mA/*txt";
open(OUT1,">$dir/m6A_gene.txt");
foreach my$line(@chr){
   $line=~/.+\/(chr_\d+)\.txt/;$chr=$1;
   open(TXT,"<$line");
   @gff=<TXT>; close TXT; undef@m6A_corrd;
   foreach my$i(@gff){ if($i=~/^(chr_\d+)\t.\t.+?\t(\d+)\t(\d+)\t/){ push@m6A_corrd,$2;} }
   $j=0;  
    while($j<@m6A_corrd){ $n=$j+1;
       $gene++; 
       while($n<@m6A_corrd){ 
          if($m6A_corrd[$n]-$m6A_corrd[$j]<=50){ $number++; $n++;
          }else{ if($number>=3){ print OUT1"\n$chr\tgene$gene\tgroup1:"; $x=$j; while($x<$n){print OUT1"\t$m6A_corrd[$x]"; $x++;}} last; }
       }
       my$mark=$n-1;
       while($n<@m6A_corrd){ if($m6A_corrd[$n]-$m6A_corrd[$j]<140){ $n++;  }else{ last;} }
       $number=0; $m=$n; $mid=$m6A_corrd[$j]+(($m6A_corrd[$mark]-$m6A_corrd[$j])/2);
       while($n<@m6A_corrd){
          if($m6A_corrd[$n]-$mid>=150 && $m6A_corrd[$n]-$mid<=220){
              $number++; $n++;  
          }else{ if($number>=3){ print OUT1"\tgroup2:";   $x=$m; while($x<$n){print OUT1"\t$m6A_corrd[$n]"; $x++;}} last; }
       }
       $number=0;  $q=$n;
       while($n<@m6A_corrd){ if($m6A_corrd[$n]-$mid<330){$n++; }else{last;} }
       while($n<@m6A_corrd){ 
          if($m6A_corrd[$n]-$mid>=350 && $m6A_corrd[$n]-$mid<=420){
              $number++; $n++;
          }else{ if($number>=3){ print OUT1"\tgroup3:"; $x=$q; while($x<$n){print OUT1"\t$m6A_corrd[$x]"; $x++; }} last; }
       }
       $number=0; $p=$n;
       while($n<@m6A_corrd){ if($m6A_corrd[$n]-$mid<520){$n++; }else{last;} }
       while($n<@m6A_corrd){ 
          if($m6A_corrd[$n]-$mid>=550 && $m6A_corrd[$n]-$mid<=620){
              $number++; $n++;
          }else{ if($number>=3){ print OUT1"\tgroup4:"; $x=$p; while($x<$n){print OUT1"\t$m6A_corrd[$x]"; $x++;}} last; }
       }
       while($n<@m6A_corrd){ if($m6A_corrd[$n]-$mid<700){$n++; print OUT1"\t$m6A_corrd[$n]";}else{last;} }
       $j=$n+1;
   }
}
undef %m6A_corrd;
#--------------------------------------------------
print"STEP2: Filter line  less than 2 m6A corrd and import into m6A_gene_filtered.txt\nPlease wait ...............\n";
open(IN2,"<$dir/m6A_gene.txt");
open(OUT2,">$dir/m6A_gene_filterd.txt");
foreach my$i(<IN2>){
    @m6A=split(/\t/,$i);
    my$n=0;
    foreach my$j(@m6A){ if($j=~/^\d+$/){$n++;} }
    if($n>=5){print OUT2"$i";}
}

#--------------------------------------------------
#chr_001 gene1   group1: 685     689     698
print "STEP3: Start Checking gene model with 6mA in m6A_gene_filtered.txt\nPlease wait .................\n ";
open(IN3,"<$dir/m6A_gene_filterd.txt") or die;
foreach(<IN3>){ 
    if(/^(chr_\d+)\s/){ 
      $m6A_corrd{$_}=$_;
      open(TXT,">>$dir/chr_filter/$1.txt"); 
      print TXT"$_"; 
      close TXT;
    } 
}
my($gene_start, $gene_end, @m6A, $mark, $gene_id);
undef@m6A; $m=0;
#open(IN4,"<$dir/2-upd-Genome-GFF3-latest-2.gff3");
open(IN4,"<$dir/06Merge_Deleted636-81_UTR2901-550_To13gff3.gff3");
open(OUT3,">$dir/m6A_gene_filterd_checked_gene.gff");
#open(OUT3,">$dir/2-upd-Genome-GFF3-latest-2.gff3_6mA_check.gff");
print "STEP4: Check gene with regular 6mA distribution like nuclesome\nPlease wait ................\n";
foreach my$i(<IN4>){ chomp $i;
   if($i=~/^(chr_\d+)\s.+?gene\t(\d+)\t(\d+)\t.\t\+\t.+?ID=(.+?);Name/){
      $gene_start=$2; $gene_end=$3; $gene_id=$4;
      open(TXT,"<$dir/chr_filter/$1.txt");
      foreach my$j(<TXT>){ 
         @m6A=split(/\t/,$j); undef@m6A_corrd;
         foreach(@m6A){if(/^\d+$/){push@m6A_corrd,$_;} }
         $number=0; 
         foreach(@m6A_corrd){ 
            if($gene_end-$gene_start>1000){
              if($_>$gene_start && $_<$gene_start+1000){$number++; $mark=0; }else{ $mark=1;  } 
            }else{
              if($_>$gene_start && $_<$gene_end){$number++; $mark=0; }else{ $mark=1;  }
            }
         }
         print OUT3"$i\n$j" if$number>=3;
         delete $m6A_corrd{$j} if$number>=3; 
         if($number>=3){$gene{$gene_id}=$gene_id; $m++;}
      }
      close TXT;
   }elsif($i=~/^(chr_\d+)\s.+?gene\t(\d+)\t(\d+)\t.\t\-\t.+?ID=(.+?);Name/){
      $gene_start=$2; $gene_end=$3; 
      open(TXT2,"<$dir/chr_filter/$1.txt") or die;
      foreach my$j(<TXT2>){
         @m6A=split(/\t/,$j); undef@m6A_corrd;
         foreach(@m6A){if(/^\d+$/){push@m6A_corrd,$_;} }
         $number=0; 
         foreach(@m6A_corrd){ 
           if($gene_end-$gene_start>1000){
              if($_>($gene_end-1000) && $_<$gene_end){$number++; $mark=0; }else{ $mark=1;  } 
           }else{
              if($_>$gene_start && $_<$gene_end){$number++; $mark=0; }else{ $mark=1;  }
           }
         }
         print OUT3"$i\n$j" if$number>=3;
         delete $m6A_corrd{$j} if$number>=3;
         if($number>=3){$gene{$gene_id}=$gene_id; $m++;}
      }
      close TXT2; 
   }
}
close IN4; close OUT3;
print "$m 6mA group distribute downstream of TSS\n";
my$n=keys%m6A_corrd;
print"number of m6A group misdistribution: $n\n";
#-----------------------------------------------------------
open(IN5,"<$dir/06Merge_Deleted636-81_UTR2901-550_To13gff3.gff3") or die;
open(OUT4,">$dir/06gff3_m6A_unregular_ingene.txt");
open(OUT5,">$dir/06gff3_m6A_unregular_ingene.temp");
print "STEP5: Check gene with unregular 6mA distribution\nPlease wait .............\n";
push my@gene_corrd,0;
foreach(<IN5>){ 
    if(/^(chr_\d+)\s.+?gene\t(\d+)\t(\d+)\t.\t.\t.+?ID=(.+?);Name/){
        open(TXT,">>$dir/chr_gene/$1.txt") or die;
        print TXT "$_"; 
        close TXT;
        push@gene_corrd,$2,$3;
    }
}
my@gene_chr=glob"$dir/chr_gene/*txt";
my($gene_before, $gene_ungelar_6mA); $chr="chr_001"; 
foreach(@gene_chr){
    open(GENE,"<$_");
    $gene_before=0; $mark=0; unshift@gene_corrd,0;
    foreach my$i(<GENE>){
       if($i=~/^(chr_\d+)\s.+?gene\t(\d+)\t(\d+)\t.\t\+\t.+?ID=(.+?);Name/){
          if($1 eq $chr){}else{shift@gene_corrd; shift@gene_corrd; shift@gene_corrd; unshift@gene_corrd,0; }
          if($2!=$gene_corrd[1]){ shift@gene_corrd; }
          $gene_start=$2; $gene_end=$3; $chr=$1; $gene_before=$gene_corrd[0];
          print OUT5"gene_id:$4\tgene_corrd0:$gene_corrd[0] gene_corrd1:$gene_corrd[1] gene_corrd2:$gene_corrd[2]\n$i";
          if(!exists$gene{$4}){
            foreach my$j(sort{$a<=>$b}keys%m6A_corrd){
              if($j=~/$chr\s/){@m6A=split(/\t/,$j);
                  foreach my$m(0 .. $#m6A){
                       if($m6A[$m]=~/^(\d+)$/){ 
                           if($1>$gene_before && $1<$gene_end){
                               my$corrd=$1-$gene_start; $m6A[$m]="$1|$corrd\t"; $number++;
                           }  
                       } $m++; 
                  }
                  if($number>3){$gene_ungelar_6mA++;}
                  print OUT4"$i$gene_before\t$gene_start\t$gene_end\n@m6A\n" if $number>3; $number=0;
              }
            }
          }
          shift@gene_corrd;
       }elsif($i=~/^(chr_\d+)\s.+?gene\t(\d+)\t(\d+)\t.\t\-\t.+?ID=(.+?);Name/){
          if($1 eq $chr){}else{shift@gene_corrd; shift@gene_corrd; shift@gene_corrd; unshift@gene_corrd,0;}
          if($2!=$gene_corrd[1]){ shift@gene_corrd; }
          $gene_start=$2; $gene_end=$3; $chr=$1; $gene_before=$gene_corrd[3];
          print OUT5"gene_id:$4\tgene_corrd3:$gene_corrd[3] gene_corrd1:$gene_corrd[1] gene_corrd2:$gene_corrd[2]\n$i";
          if(!exists$gene{$4}){
            foreach my$j(sort{$a<=>$b}keys%m6A_corrd){
              if($j=~/$chr\s/){@m6A=split(/\t/,$j);
                  foreach my$m(0 .. $#m6A){
                       if($m6A[$m]=~/^(\d+)$/){
                           if($1>$gene_start && $1<$gene_before){
                               my$corrd=$gene_end-$1; $m6A[$m]="$1|$corrd\t"; $number++;
                           }
                       } $m++;
                  }
                  if($number>3){$gene_ungelar_6mA++;}
                  print OUT4"$i$gene_before\t$gene_start\t$gene_end\n@m6A\n" if $number>3; $number=0;
              }
            }
          }
          shift@gene_corrd;
       }
    }
}
print"$gene_ungelar_6mA gene with ungular 6mA distribution\nFINISHED! THANK YOU!! BYE\n__________________________\n";









