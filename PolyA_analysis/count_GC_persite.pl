use strict;
my$dir=$ENV{'PWD'};
my(%motif, $name, $fa, @sub, $sub, $sub, %count);
open(IN1,"<$ARGV[0]"); #$fasta
open(OUT1,">$ARGV[0]_ratio");

my$n=0;
foreach(<IN1>){
   if(/>/){  #$name=$1; 
    }else{ chomp; @sub=split(//,$_); $n++;
           my$i=0; while($i<@sub){ $count{$i}{$sub[$i]}++; $i++; }  
   }
}

#my$n=20937;
foreach my$a(keys%count){
    print OUT1"$a\t";
    foreach my$b(keys%{$count{$a}}){ my$ratio=$count{$a}{$b}/$n*100; print OUT1"$b\t$ratio\t"; }
    print OUT1 "\n";
}
