$dir=$ENV{'PWD'};
open(IN1,"<","$dir/Whole-genome_Epigenetic-mark_density.txt");
open(IN2,"<","$dir/Whole-genome_Epigenetic-mark_density_RandomForest-prediction.txt");

foreach(<IN1>){
   if(/^(((chr_\d+)-(\d+)-(\d+)-(\d+))\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?))\s+/){ #chr-cluster-start-end   6mA     H3K4me3 H2A.Z   nucleosome  chr_062-90-51994-52236  0       0       0       3
       $pos=$2; $m6A{$2}=$7; $H3K4me3{$2}=$8; $H2A{$2}=$9; $nucleosome{$2}=$10;
       $total{$2}=$1;
   }
}

open(OUT1,">","$dir/Merge-prediction-TSS-density.txt");
print OUT1"chr-cluster-start-end\t6mA\tH3K4me3\tH2A.Z\tnucleosome\ttype\n";
foreach(<IN2>){
   chomp;
   if(/^"(.+?)","(.)"/){ #"chr_062-63-83693-84909","1"
      print OUT1 "$total{$1}\t$2\n";
   }
}
