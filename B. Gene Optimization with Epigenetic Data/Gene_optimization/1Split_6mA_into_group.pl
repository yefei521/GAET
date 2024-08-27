use strict;
my$dir=$ENV{'PWD'};
my(@chr, @gff, @m6A_corrd, $j, $n, $number, $gene, $chr, $mid, $m, $p, $q, $x, @m6A, %m6A_corrd, %gene);
open(IN1,"<$dir/veg_6mA.sorted.gff");
my$cmd=`rm chr/*`;
print"$cmd";
open(IN2,"<$dir/chrname");
my@m6A=<IN1>;
foreach(<IN2>){ chomp; $chr=$_;
  open(TXT,">$dir/chr/$chr.txt");
  foreach my$line(@m6A){
   if($line=~/^($chr)\t.\t.+?\t(\d+)\t(\d+)\t/){  print TXT "$line"; }
  }
}

print "STEP1: Split m6A into group.\nPlease wait ........\n";
@chr=glob"$dir/chr/*txt";
open(OUT1,">$dir/m6A_group.bed");
foreach my$line(@chr){
   $line=~/.+\/(chr_\d+)\.txt/;$chr=$1;
   open(TXT,"<$line");
   @gff=<TXT>; close TXT; undef@m6A_corrd;
   foreach my$i(@gff){ if($i=~/^(chr_\d+)\t.\t.+?\t(\d+)\t(\d+)\t/){ push@m6A_corrd,$2;} }
   $j=0;
    while($j<@m6A_corrd){ $n=$j+1; $number=0;
       $gene++;
       while($n<@m6A_corrd){
          if($m6A_corrd[$n]-$m6A_corrd[$j]<=75){ $number++; $n++;
          }else{ if($number>=3){  print OUT1"$chr\t.\t6mA"; $x=$j; while($x<$n){print OUT1"\t$m6A_corrd[$x]"; $x++;} print OUT1"\n";  } last; }
       }
       $j=$n;
   }
}
