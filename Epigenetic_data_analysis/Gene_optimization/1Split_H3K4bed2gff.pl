use strict;
open(IN1,"<$ARGV[0]");
my$dir=$ENV{'PWD'};
open(IN2,"<$dir/chrname");
my@GFF=<IN1>;
foreach(<IN2>){
  chomp; my$chr=$_;
  open(OUT1,">","$dir/H3K4me3/$chr.txt");
  foreach(@GFF){
   if(/^($chr)\t(\d+)\t(\d+)\t/){
       my$n=($2+$3)/2; print OUT1"$chr\t.\tH3K4me3\t$n\t$n\n";
    }
  }
}
my$cmd=`cat ./H3K4me3/*txt > H3K4me3_group.gff3`;
print"$cmd";
my$cmd=`rm -r ./H3K4me3/*`;
print"$cmd";
