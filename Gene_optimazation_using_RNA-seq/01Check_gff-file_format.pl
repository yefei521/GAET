use strict;
my$dir=$ENV{'PWD'};
my($supercontig, $gene_id, %mRNA, %gff, @CDS, $chr, $gene_start, $gene_end, $ori, $gene_name, @UTR, @UTR5, @UTR3);

#open(IN1,"<$dir/Thermophila_tetrahymena_2021_UTR_manualcheck_220224.gff3");
#open(OUT1,">$dir/01Error_gene.txt");
 open(IN1,"<$ARGV[0]");
 open(OUT1,">$ARGV[0]_01Check");

while(<IN1>){
       if(/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+.+?ID=(chr_\d+);/){ $supercontig++;   #print OUT"$_";
       }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5; $mRNA{$gene_id}=$6; $gff{$gene_id}.=$_; #if(exists$hash{$5}){}else{print OUT"$_"; $gene_id=$5;}
       }elsif(/chr_\d+.+?mRNA.+?Parent=$gene_id/){  $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?UTR.+?Parent=$gene_id/){   $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?CDS.+?Parent=$gene_id/){   $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?MANUAL_ENTRY.+?ID=(.+?);Name/){ #print OTR"$_";
       }elsif(/chr_\d+.+?rDNA.+?ID=(.+?);Name/){  #print OUT"$_"; push@gff,$1; $gff{$1}=$_;
       }
    }

print "chromosome number:$supercontig\n";
#chr_001 AUGUSTUS        gene    8269    10024   1755    .       +       .       "ID=g18412;Name=TTHERM_00161850;Note=""triose-phosphate transporter family protein"""
#chr_001 AUGUSTUS        mRNA    8269    10024   1755    .       +       .       ID=g18412.t1;Parent=g18412
#chr_001 AUGUSTUS        five_prime_UTR  8269    8638    369     .       +       .       ID=g18412.t1.utr;Parent=g18412.t1
#chr_001 AUGUSTUS        CDS     8639    9838    1199    .       +       .       ID=g18412.t1.cds;Parent=g18412.t1
#chr_001 AUGUSTUS        three_prime_UTR 9839    10024   185     .       +       .       ID=g18412.t1.utr;Parent=g18412.t1
foreach my$name(keys%gff){
   my@temp=split(/\n/,$gff{$name}); undef@CDS; undef@UTR; undef@UTR5; undef@UTR3; my$n=0;
   foreach my$line(@temp){

       if($line=~/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){
          $chr=$1; $gene_start=$2; $gene_end=$3; $ori=$4; $gene_id=$5; $gene_name=$6;
       }elsif($line=~/$chr.+?mRNA\s+($gene_start)\s+($gene_end)\s+.\s+(.)\s+.+?Parent=$gene_id/){
       }elsif($line=~/$chr.+?five_prime_UTR\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?Parent=$gene_id/){ 
          push@UTR5,$1,$2;
       }elsif($line=~/$chr.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?Parent=$gene_id/){
          push@CDS,$1,$2;
       }elsif($line=~/$chr.+?three_prime_UTR\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?Parent=$gene_id/){
          push@UTR3,$1,$2;
       }else{print "Error:\n\t$line\n";}
    }
    if($ori eq "+"){
          if( @UTR5==0 && @UTR3==0){    if($gene_start==$CDS[0] && $gene_end==$CDS[-1]){ pop@CDS; shift@CDS; while($n<$#CDS){ if($CDS[$n]<$CDS[$n+1]){}else{print OUT1"$ori\t$name\t1\n"; }  $n+=2;} }else{print OUT1"$ori\t$name\t2\n";}
          }elsif(@UTR5>0 && @UTR3==0 ){ if($gene_start==$UTR5[0]&& $gene_end==$CDS[-1] && ($UTR5[-1]+1)==$CDS[0]){ pop@CDS;shift@CDS; while($n<$#CDS){ if($CDS[$n]<$CDS[$n+1]){}else{print OUT1"$CDS[$n]\t$CDS[$n+1]\t$ori\t$name\t3\n"; } $n+=2; } }else{print OUT1"$ori\t$name\t4\n";}
          }elsif(@UTR5==0 && @UTR3>0 ){ if($gene_end==$UTR3[-1]&& $gene_start==$CDS[0] && ($UTR3[0]-1)==$CDS[-1]){ pop@CDS;shift@CDS; while($n<$#CDS){ if($CDS[$n]<$CDS[$n+1]){}else{print OUT1"$ori\t$name\t5\n"; } $n+=2; } }else{print OUT1"$ori\t$name\t6\n";}
          }elsif(@UTR5>0 && @UTR3>0){ if($gene_start==$UTR5[0]&& $gene_end==$UTR3[-1] && ($UTR5[-1]+1)==$CDS[0] && ($UTR3[0]-1)==$CDS[-1]){ pop@CDS;shift@CDS; while($n<$#CDS){ if($CDS[$n]<$CDS[$n+1]){}else{print OUT1"$ori\t$name\t7\n"; }$n+=2; } }else{print OUT1"$ori\t$name\t8\n";}
          }
    }elsif($ori eq "-"){
          @UTR=@UTR5; @UTR5=@UTR3; @UTR3=@UTR;
          if( @UTR5==0 && @UTR3==0){    if($gene_start==$CDS[0] && $gene_end==$CDS[-1]){ pop@CDS; shift@CDS; while($n<$#CDS){ if($CDS[$n]<$CDS[$n+1]){}else{print OUT1"$ori\t$name\t1\n"; }  $n+=2;} }else{print OUT1"$ori\t$name\t2\n";}
          }elsif(@UTR5>0 && @UTR3==0 ){ if($gene_start==$UTR5[0]&& $gene_end==$CDS[-1] && ($UTR5[-1]+1)==$CDS[0]){ pop@CDS;shift@CDS; while($n<$#CDS){ if($CDS[$n]<$CDS[$n+1]){}else{print OUT1"$CDS[$n]\t$CDS[$n+1]\t$ori\t$name\t3\n"; } $n+=2; } }else{print OUT1"$ori\t$name\t4\n";}
          }elsif(@UTR5==0 && @UTR3>0 ){ if($gene_end==$UTR3[-1]&& $gene_start==$CDS[0] && ($UTR3[0]-1)==$CDS[-1]){ pop@CDS;shift@CDS; while($n<$#CDS){ if($CDS[$n]<$CDS[$n+1]){}else{print OUT1"$ori\t$name\t5\n"; } $n+=2; } }else{print OUT1"$ori\t$name\t6\n";}
          }elsif(@UTR5>0 && @UTR3>0){ if($gene_start==$UTR5[0]&& $gene_end==$UTR3[-1] && ($UTR5[-1]+1)==$CDS[0] && ($UTR3[0]-1)==$CDS[-1]){ pop@CDS;shift@CDS; while($n<$#CDS){ if($CDS[$n]<$CDS[$n+1]){}else{print OUT1"$ori\t$name\t7\n"; }$n+=2; } }else{print OUT1"$ori\t$name\t8\n";}
          }
   }
}
