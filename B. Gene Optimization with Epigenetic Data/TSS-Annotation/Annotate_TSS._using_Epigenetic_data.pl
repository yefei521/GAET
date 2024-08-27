##~/miniconda3/envs/homer/bin/annotatePeaks.pl /apps/users/yefei/learning/SB210-ATAC-seq/06Call_peak_100bp/Homer_peak none  -gtf TGD2024_renamed_TSSanno_addAS-NAT-NC.gtf > NFR_peaks.bed
my$dir=$ENV{'PWD'};
open(IN1,"<","$dir/Manual_check-total-gene.gff3_Right_UTR.gff3") or die;
open(IN2,"<","$dir/77803_GGACTCCT_S5_L008_rmdup_100bp.bed_q1_peak_summits.bed") or die;
open(IN3,"<","$dir/77804_TAGGCATG_S6_L008_100bp.bed_peak_summits.bed") or die;
open(IN4,"<","$dir/NFR_peaks.bed") or die;
open(IN5,"<","$dir/Whole-genome_Epigenetic-mark_density_RandomForest-prediction.txt");
foreach(<IN1>){
   if(/^(chr_\d+)\s.+?(gene)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);Name=(.+?);ManInfo=(.+?)\s/){ #chr_088 GSAman  gene    461906  463746 -  ID=YF00013595;Name=YF00013595;ManInfo=Manual-Check
      $gff{$1}{$6}{$2}.=$_; push@{$pos{$1}{$6}{$2}},$3,$4; $num++; $ori{$6}=$5;
   }elsif(/^(chr_\d+)\s.+?(mRNA)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+Parent=(.+?);ID=(.+?)\s/){ #chr_088 GSAman  mRNA    461906  463746 -  Parent=YF00013595;ID=YF00013595.t1
      $gff{$1}{$6}{$2}.=$_; push@{$pos{$1}{$6}{$2}},$3,$4;
   }elsif(/^(chr_\d+)\s.+?(CDS)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+Parent=(.+?);ID=(.+?)\s/){ #chr_088 GSAman  CDS     462005  462413  -   Parent=YF00013595.t1;ID=YF00013595.t1.cds
      $gff{$1}{$6}{$2}.=$_; push@{$pos{$1}{$6}{$2}},$3,$4;
   }elsif(/^(chr_\d+)\s.+?(UTR)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+Parent=(.+?);ID=(.+?)\s/){
      $gff{$1}{$6}{$2}.=$_; push@{$pos{$1}{$6}{$2}},$3,$4;
   }elsif(/^(chr_\d+)\s.+?(exon)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+Parent=(.+?);ID=(.+?)\s/){
      $gff{$1}{$6}{$2}.=$_; push@{$pos{$1}{$6}{$2}},$3,$4;
   }else{print "$_"; } 
}

print "gene number:$num\n";
$num=0;
foreach(<IN2>){
    if(/^(chr_\d+)\s+(\d+)\s+(\d+)\s+(.+?)\s/){ #chr_001 238     239     77803_GGACTCCT_S5_L008_rmdup_100bp.bed_q1_peak_peak_1   3.76082
        $peak{$1}{$2}++; $num++;
    }
}
print "broad peak number:$num\n";
$num=0;

foreach(<IN3>){
    if(/^(chr_\d+)\s+(\d+)\s+(\d+)\s+(.+?)\s/){ #chr_001 238     239     77803_GGACTCCT_S5_L008_rmdup_100bp.bed_peak_peak_1   3.76082
        $high_peak{$1}{$2}++; $num++;
    }
}

print "high-confidence peak number:$num\n";
$num=0;

foreach(<IN4>){
    if(/^(chr_\d+)\s+(\d+)\s+(\d+)\s+(.+?)\s+.+?\s+(.)\s/){ #chr_001 238     239     77803_GGACTCCT_S5_L008   3.76082
       if($5 eq "+"){ $start=$2+4; $end=$3; }else{ $end=$3-5; $start=$2; }
           foreach $i($start .. $end){ $cov{$1}{$i}++;  } 
    }
}

foreach(<IN5>){
   chomp;
   if(/^"(.+?)","(.)"/){ #"chr_062-63-83693-84909","1"
      $TSS_r{$1}=$2;
   }
}


foreach$chr(keys%gff){
    foreach$name(keys%{$gff{$chr}}){ 
         if($ori{$name} eq "+"){
              @temp=@{$pos{$chr}{$name}{"gene"}};
              $start=$temp[0];
              foreach$i($start-500 .. $start+40){ 
                  if(exists$high_peak{$chr}{$i}){
                       $TSS{$name}=$i; #print "$name\t+\t$start\t$i\n"; 
                       last;}
                  if(exists$peak{$chr}{$i}){
                       $TSS{$name}=$i; #print "$name\t+\t$start\t$i\n"; 
                       last;}
              }
          }elsif($ori{$name} eq "-"){
              $start=${$pos{$chr}{$name}{"gene"}}[1]; 
              foreach$i($start-40 .. $start+500){   
                  if(exists$high_peak{$chr}{$i}){
                     $TSS{$name}=$i; #print "$name\t-\t$start\t$i\n"; 
                     last;}
                  if(exists$peak{$chr}{$i}){
                     $TSS{$name}=$i; #print "$name\t-\t$start\t$i\n";
                     last;}
              }
          }
   }
}

open(OUT1,">","$dir/TSS.bed");
open(OUT2,">","$dir/eTSS.bed");
foreach $chr(keys%gff){
   foreach $name(keys%{$gff{$chr}}){
      @temp=@{$pos{$chr}{$name}{"gene"}};
      if(exists$TSS{$name}){
          print OUT1"$chr\t$TSS{$name}\t$TSS{$name}\t$name\t$high_peak{$chr}{$TSS{$name}}\t$ori{$name}\n";
          print OUT2"$chr\t$TSS{$name}\t$TSS{$name}\t$name\t$high_peak{$chr}{$TSS{$name}}\t$ori{$name}\n";
          $num++;
       }else{ if($ori{$name} eq "+"){
              $start=$temp[0]; undef%bin;
              foreach$i($start-300 .. $start+20){ if(exists$cov{$chr}{$i}){ $bin{$i}+=$cov{$chr}{$i}; } }
                  @sort_bin=sort{$bin{$b}<=>$bin{$a}}keys%bin;
                  print OUT1"$chr\t$sort_bin[0]\t$sort_bin[0]\t$name\t.\t$ori{$name}\n" if $bin{$sort_bin[0]}>10;
             }elsif($ori{$name} eq "-"){
              $start=$temp[0]; undef%bin;
              foreach$i($start-300 .. $start+20){ if(exists$cov{$chr}{$i}){ $bin{$i}+=$cov{$chr}{$i};  } }
                  @sort_bin=sort{$bin{$b}<=>$bin{$a}}keys%bin;
                  print OUT1"$chr\t$sort_bin[0]\t$sort_bin[0]\t$name\t.\t$ori{$name}\n" if $bin{$sort_bin[0]}>10;
            }
       }
   }
}

print "Annotated TSS:$num";

