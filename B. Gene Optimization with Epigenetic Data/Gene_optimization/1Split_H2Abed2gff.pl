use strict;
open(IN1,"<$ARGV[0]");
my$dir=$ENV{'PWD'};
open(IN2,"<$dir/chrname");
my@GFF=<IN1>;
foreach(<IN2>){
  chomp; my$chr=$_;
  open(OUT1,">","$dir/H2A.Z/$chr.txt");
  foreach(@GFF){
   if(/^($chr)\t(\d+)\t(\d+)\t/){
       my$n=($2+$3)/2; print OUT1"$chr\t.\tH2A.Z\t$n\t$n\n";
    }
  }
}
my$cmd=`cat ./H2A.Z/*txt > H2A.Z_group.gff3`;
print"$cmd";
my$cmd=`rm -r ./H2A.Z/*`;
print"$cmd";
