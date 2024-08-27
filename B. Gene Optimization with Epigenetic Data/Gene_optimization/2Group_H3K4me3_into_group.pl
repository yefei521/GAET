use strict;
my$dir=$ENV{'PWD'};
my(@chr, @gff, @H3K4me3_corrd, $j, $n, $number, $gene, $chr, $mid, $m, $p, $q, $x, @H3K4me3, %H3K4me3_corrd, %gene);
open(IN1,"<$ARGV[0]");
my$cmd=`rm chr/*`;
print"$cmd";
open(IN2,"<$dir/chrname");
my@H3K4me3=<IN1>;
foreach(<IN2>){ chomp; $chr=$_;
  open(TXT,">$dir/chr/$chr.txt");
  foreach my$line(@H3K4me3){
   if($line=~/^($chr)\t.\t.+?\t(.+?)\t(.+?)\s/){  print TXT "$line"; }
  }
}

print "STEP1: Split H3K4me3 into group.\nPlease wait ........\n";
@chr=glob"$dir/chr/*txt";
open(OUT1,">$dir/H3K4me3_group.bed");
foreach my$line(@chr){
   $line=~/.+\/(chr_\d+)\.txt/;$chr=$1;
   open(TXT,"<$line");
   @gff=<TXT>; close TXT; undef@H3K4me3_corrd;
   foreach my$i(@gff){ if($i=~/^(chr_\d+)\t.\tH3K4me3\t(.+?)\t(.+?)\s/){ push@H3K4me3_corrd,$2;} }
   $j=0;
    while($j<@H3K4me3_corrd){ $n=$j+1; $number=0;
       $gene++;
       while($n<@H3K4me3_corrd){
          if($H3K4me3_corrd[$n]-$H3K4me3_corrd[$j]<=75){ $number++; $n++;
          }else{ if($number>=6){  print OUT1"$chr\t.\tH3K4me3"; $x=$j; while($x<$n){print OUT1"\t$H3K4me3_corrd[$x]"; $x++;} print OUT1"\n";  } last; }
       }
       $j=$n;
   }
}
