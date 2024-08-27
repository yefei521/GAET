use strict;
my$dir=$ENV{'PWD'};
my(@chr, @gff, @H2A_corrd, $j, $n, $number, $gene, $chr, $mid, $m, $p, $q, $x, @H2A, %H2A_corrd, %gene);
open(IN1,"<$ARGV[0]"); #
my$cmd=`rm chr/*`;
print"$cmd";
open(IN2,"<$dir/chrname");
my@H2A=<IN1>;
foreach(<IN2>){ chomp; $chr=$_;
  open(TXT,">$dir/chr/$chr.txt");
  foreach my$line(@H2A){
   if($line=~/^($chr)\t.\t.+?\t(.+?)\t(.+?)\s/){  print TXT "$line"; }
  }
}

print "STEP1: Split H2A into group.\nPlease wait ........\n";
@chr=glob"$dir/chr/*txt";
open(OUT1,">$dir/H2A.Z_group.bed");
foreach my$line(@chr){
   $line=~/.+\/(chr_\d+)\.txt/;$chr=$1;
   open(TXT,"<$line");
   @gff=<TXT>; close TXT; undef@H2A_corrd;
   foreach my$i(@gff){ if($i=~/^(chr_\d+)\t.\tH2A.Z\t(.+?)\t(.+?)\s/){ push@H2A_corrd,$2;} }
   $j=0;
    while($j<@H2A_corrd){ $n=$j+1; $number=0;
       $gene++;
       while($n<@H2A_corrd){
          if($H2A_corrd[$n]-$H2A_corrd[$j]<=75){ $number++; $n++;
          }else{ if($number>=6){  print OUT1"$chr\t.\tH2A"; $x=$j; while($x<$n){print OUT1"\t$H2A_corrd[$x]"; $x++;} print OUT1"\n";  } last; }
       }
       $j=$n;
   }
}
