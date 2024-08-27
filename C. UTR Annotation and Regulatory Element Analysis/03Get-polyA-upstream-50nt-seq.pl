$dir=$ENV{'PWD'};
open(IN1,"<","$dir/WT_mix_seqkit2DNA_MaxIntron2k.sam_polya.txt");
open(OUT1,">","$dir/polyA-up-downstream-50nt-seq.bed");
while(<IN1>){
   if(/^(chr_\d+)\s+(\d+)\s+(\d+)\s+.+?\d+\s+(.)\s+/){ #chr_001 235664  236255  3dc132ab-8754-42da-9d29-c2b542cffdaa    14      -       TCTAAAAAAAAAAAAGA       polya:AAAAAAAAAAAAGA
       if($4 eq '+'){ $start=$3-50; $end=$3+50; print OUT1"$1\t$start\t$end\t.\t.\t+\n"; 
      }elsif($4 eq '-'){ $start=$2-50; $end=$2+50; print OUT1"$1\t$start\t$end\t.\t.\t-\n"; }
   }
}

