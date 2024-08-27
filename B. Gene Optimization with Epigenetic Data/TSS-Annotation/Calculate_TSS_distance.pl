$dir=$ENV{'PWD'};
open(IN1,"<","$dir/Homer_peak_NFR.bed");
open(IN2,"<","$dir/TSS_annotated_gene.txt");
open(IN3,"<","$dir/TGD2024_renamed_TSSanno_addAS-NAT-NC.gtf");
open(IN4,"<","$dir/TGD2024_renamed.gff3");
foreach(<IN2>){
   if(/(TTHERM_\d+)\s/){
     $TSS_annoT{$1}++; #print "$1\n";
   }
}

foreach(<IN4>){
   if(/^(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(g\d+);Name=(TTHERM_\d+);UserName/){
      $T2g{$6}=$5; $TSS_annog{$5}++ if exists $TSS_annoT{$6}; 
   }
}

foreach(<IN3>){
   if(/^(chr_\d+)\s.+?transcript\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+transcript_id "(g\d+).t1";/){
      if(exists$TSS_annog{$5} && $4 eq "+"){ $TSSf{$1}{$5}=$2;   }
      if(exists$TSS_annog{$5} && $4 eq "-"){ $TSSr{$1}{$5}=$3;   }
   }
}

foreach(<IN1>){
   if(/^(chr_\d+)\s+(\d+)\s+(\d+)\s/){ 
       undef@temp;
       foreach $name(keys%{$TSSf{$1}}){
         $posTSS=$TSSf{$1}{$name}; $distance=abs(($3+$2)/2-$posTSS);
         push@temp,$distance;
       }
       @sorttemp=sort{$a<=>$b}@temp;
       $min_dis1=$sorttemp[0];
       undef@temp;
       foreach $name(keys%{$TSSr{$1}}){
         $posTSS=$TSSr{$1}{$name}; $distance=abs(($3+$2)/2-$posTSS);
         push@temp,$distance;
       }
       @sorttemp=sort{$a<=>$b}@temp;
       $min_dis2=$sorttemp[0];
       if($min_dis1<$min_dis2){$min=$min_dis1}else{$min=$min_dis2}
       print "$min\n";
    }
}
