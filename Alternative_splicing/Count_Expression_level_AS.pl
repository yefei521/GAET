$dir=$ENV{'PWD'};
open(IN1,"<","$dir/TGD2024_renamed_addAS-NAT-NC.gtf");
#chr_001 GSAman  transcript      258     1182    .       +       .       gene_id "g18410.t2"; transcript_id "g18410";
#chr_001 GSAman  exon    258     453     .       +       .       gene_id "g18410.t2"; transcript_id "g18410";
#chr_001 GSAman  exon    505     1182    .       +       .       gene_id "g18410.t2"; transcript_id "g18410";
foreach(<IN1>){
    if(/^(chr_\d+)\s+.+?\stranscript\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+gene_id ".+?.(t.)"; transcript_id "(.+?)";/){
         $chr=$1; $start=$2; $end=$3; $ori=$4; $gnameTn=$5; $gname=$6; 
     }elsif(/^(chr_\d+)\s+.+?\sexon\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+gene_id ".+?(t.)"; transcript_id "(.+?)";/){
         $chr=$1; $start=$2; $end=$3; $ori=$4; $gnameTn=$5; $gname=$6;
         push@{$gff{$chr}{$gname}{$gnameTn}},"$start-$end";
         $exon_len{$gname}{$gnameTn}+=$end-$start+1;
     }
}
 $TotalCount{"C10_1"}=70818126; $TotalCount{"C4"}=43101447; $TotalCount{"C5"}= 61052530; $TotalCount{"C6"}= 63088370; $TotalCount{"C8"}=59185127 ; $TotalCount{"C10"}= 61235682; 
 $TotalCount{"C18"}= 69833166; $TotalCount{"veg"}= 75613504; $TotalCount{"S24"}= 48948585; 
#Geneid  Chr     Start   End     Strand  Length  WT-C10  WT-C18  WT_C10  WT_C4   WT_C5   WT_C6   WT_C8   WT_S24  WT_veg
#g18410.t1       chr_001 256     342     +       87      0       0       0       0       0       0       0       0       0
open(IN2,"<","$dir/read.count");
foreach(<IN2>){
   @temp=split(/\t/);
   $chr=$temp[1]; $temp[0]=~/^(.+).(t\d+)/; $gname=$1; $gnameTn=$2;  $exon="$temp[2]-$temp[3]";
   push@{$Count{$chr}{$gname}{$gnameTn}{$exon}},@temp[6..14];
   
}

open(IN3,"<","$dir/TGD2024_renamed_addAS-NAT-NC.gff3");
#chr_001 GSAman  mRNA    51940   54539   .       -       .       ID=g18425.t2;Parent=g18425;ManInfo=AS;
foreach(<IN3>){
   if(/^(chr_\d+)\s+.+?\smRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=.+?(t\d+);Parent=(.+?);ManInfo=(A.+?);/){
      $chr=$1; $start=$2; $end=$3; $ori=$4; $gnameTn=$5; $gname=$6; $info=$7; $info{"$gname.$gnameTn"}=$info;
   }
}
open(OUT,">","$dir/TGD2024_AS_rpkm.txt");
print OUT "chr\tAS_id\tAS_type\tveg\tS24\tC4\tC5\tC6\tC8\tC10\tC10\n";
foreach $chr(keys%Count){
   foreach $gname(keys%{$Count{$chr}}){
       next if keys%{$Count{$chr}{$gname}}<2;
            @mRNA=@{$gff{$chr}{$gname}{"t1"}};
            #undef%ReadCount;
            foreach $exon(@mRNA){ 
                $ReadCount{$gname}{t1}{"C10_1"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[0] ; 
                $ReadCount{$gname}{t1}{"C18"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[1] ;
                $ReadCount{$gname}{t1}{"C10"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[2] ;
                $ReadCount{$gname}{t1}{"C4"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[3] ;
                $ReadCount{$gname}{t1}{"C5"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[4] ;
                $ReadCount{$gname}{t1}{"C6"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[5] ;
                $ReadCount{$gname}{t1}{"C8"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[6] ;
                $ReadCount{$gname}{t1}{"S24"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[7] ;
                $ReadCount{$gname}{t1}{"veg"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[8] ;
            }
            #$ReadCount=      
       foreach $gnameTn(keys%{$Count{$chr}{$gname}}){
            if(exists$info{"$gname.$gnameTn"}){ $info=$info{"$gname.$gnameTn"};
                 if($info=~/ASnIR/){   }    
                 @AS=@{$gff{$chr}{$gname}{$gnameTn}};
                 foreach $exon(@AS){ 
                      $ReadCount{$gname}{$gnameTn}{"C10_1"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[0] ;
                      $ReadCount{$gname}{$gnameTn}{"C18"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[1] ;
                      $ReadCount{$gname}{$gnameTn}{"C10"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[2] ; 
                      $ReadCount{$gname}{$gnameTn}{"C4"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[3] ;
                      $ReadCount{$gname}{$gnameTn}{"C5"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[4] ; 
                      $ReadCount{$gname}{$gnameTn}{"C6"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[5] ;
                      $ReadCount{$gname}{$gnameTn}{"C8"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[6] ;
                      $ReadCount{$gname}{$gnameTn}{"S24"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[7] ;
                      $ReadCount{$gname}{$gnameTn}{"veg"}+=${$Count{$chr}{$gname}{$gnameTn}{$exon}}[8] ;
                  }
                 $ReadCount{$gname}{$gnameTn}{"C10_1"}=$ReadCount{$gname}{$gnameTn}{"C10_1"}-$ReadCount{$gname}{"t1"}{"C10_1"};
                 $ReadCount{$gname}{$gnameTn}{"C18"}=$ReadCount{$gname}{$gnameTn}{"C18"}-$ReadCount{$gname}{"t1"}{"C18"};
                 $ReadCount{$gname}{$gnameTn}{"C10"}=$ReadCount{$gname}{$gnameTn}{"C10"}-$ReadCount{$gname}{"t1"}{"C10"};
                 $ReadCount{$gname}{$gnameTn}{"C4"}=$ReadCount{$gname}{$gnameTn}{"C4"}-$ReadCount{$gname}{"t1"}{"C4"};
                 $ReadCount{$gname}{$gnameTn}{"C5"}=$ReadCount{$gname}{$gnameTn}{"C5"}-$ReadCount{$gname}{"t1"}{"C5"};
                 $ReadCount{$gname}{$gnameTn}{"C6"}=$ReadCount{$gname}{$gnameTn}{"C6"}-$ReadCount{$gname}{"t1"}{"C6"};
                 $ReadCount{$gname}{$gnameTn}{"C8"}=$ReadCount{$gname}{$gnameTn}{"C8"}-$ReadCount{$gname}{"t1"}{"C8"};
                 $ReadCount{$gname}{$gnameTn}{"S24"}=$ReadCount{$gname}{$gnameTn}{"S24"}-$ReadCount{$gname}{"t1"}{"S24"};
                 $ReadCount{$gname}{$gnameTn}{"veg"}=$ReadCount{$gname}{$gnameTn}{"veg"}-$ReadCount{$gname}{"t1"}{"veg"};
                 next if !exists$exon_len{$gname}{$gnameTn};
                 $veg=int(abs($ReadCount{$gname}{$gnameTn}{"veg"}/$TotalCount{"veg"}/(abs($exon_len{$gname}{$gnameTn}-$exon_len{$gname}{"t1"})+1))*1000000000000)/1000;
                 $S24=int(abs($ReadCount{$gname}{$gnameTn}{"S24"}/$TotalCount{"S24"}/(abs($exon_len{$gname}{$gnameTn}-$exon_len{$gname}{"t1"})+1))*1000000000000)/1000;
                 $C4=int(abs($ReadCount{$gname}{$gnameTn}{"C4"}/$TotalCount{"C4"}/(abs($exon_len{$gname}{$gnameTn}-$exon_len{$gname}{"t1"})+1))*1000000000000)/1000;;
                 $C5=int(abs($ReadCount{$gname}{$gnameTn}{"C5"}/$TotalCount{"C5"}/(abs($exon_len{$gname}{$gnameTn}-$exon_len{$gname}{"t1"})+1))*1000000000000)/1000;
                 $C6=int(abs($ReadCount{$gname}{$gnameTn}{"C6"}/$TotalCount{"C6"}/(abs($exon_len{$gname}{$gnameTn}-$exon_len{$gname}{"t1"})+1))*1000000000000)/1000;
                 $C8=int(abs($ReadCount{$gname}{$gnameTn}{"C8"}/$TotalCount{"C8"}/(abs($exon_len{$gname}{$gnameTn}-$exon_len{$gname}{"t1"})+1))*1000000000000)/1000;
                 $C10=int(abs($ReadCount{$gname}{$gnameTn}{"C10"}/$TotalCount{"C10"}/(abs($exon_len{$gname}{$gnameTn}-$exon_len{$gname}{"t1"})+1))*1000000000000)/1000;
                 $C10_1=int(abs($ReadCount{$gname}{$gnameTn}{"C10_1"}/$TotalCount{"C10_1"}/(abs($exon_len{$gname}{$gnameTn}-$exon_len{$gname}{"t1"})+1))*1000000000000)/1000;
                 $C18=int(abs($ReadCount{$gname}{$gnameTn}{"C18"}/$TotalCount{"C18"}/(abs($exon_len{$gname}{$gnameTn}-$exon_len{$gname}{"t1"})+1))*1000000000000)/1000;
                 print OUT "$chr\t$gname.$gnameTn\t$info\t$veg\t$S24\t$C4\t$C5\t$C6\t$C8\t$C10\t$C10_1\n";
            }
       }
   }
}

