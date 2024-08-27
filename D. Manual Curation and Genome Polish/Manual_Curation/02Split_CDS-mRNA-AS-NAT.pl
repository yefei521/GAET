$dir=$ENV{'PWD'};
open(IN1,"<","$dir/Manual-Check_RawGFF/merge.gtf"); #merge.gtf
open(IN2,"<","$dir/Manual-Check_RawGFF/merge.fasta"); #merge.fasta
open(IN3,"<","$dir/Rename_gene/TGD2024_duplicated_nedd_check.txt_checked");
foreach $line(<IN2>){# print "$line";
    if($line=~/>(.+?)\s/){ $chrID=$1;
    }else{chomp $line; $sequence{$chrID}.=$line; }
}
foreach(<IN3>){ chomp; @temp=split(/\t/);
      foreach$line(@temp){ if($line=~/^(.+?):(AS.+)/){ $AS_name{$1}++; }elsif($line=~/^(.+?):(NAT)/){ $NAT_name{$1}++;  }else{   }  }
}
open(CDS,">","$dir/Split_CDS-mRNA-AS-NAT/CDS.gff3");
open(CDSfa,">","$dir/Split_CDS-mRNA-AS-NAT/CDS.fasta");

foreach $line(<IN1>){ #print "$line";
  if($line=~/^(chr_\d+)\s+GSAman\s+(mRNA)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);ManInfo=(.+?);/){
     $chr=$1; $region=$2; $start=$3; $end=$4; $ori=$5; $chrID=$6;  $info=$7;
     if($sequence{$chrID}=~/^ATG.+TGA$/){  
         $CDS_name{$chrID}++;
         print CDSfa "$chrID\n$sequence{$chrID}\n";
     }elsif($line=~/NA[T|t].T(\d+)/ or $line=~/NA[T|t](\d+)/){    $NAT_name{$chrID}++;  
         $line="$chr\tGSAman\t$region\t$start\t$end\t.\t$ori\t.\tID=$chrID;ManInfo=NATnT$1;\n"; 
     }elsif($info=~/AS/){ $AS_name{$chrID}++;
     }elsif($info=~/NC/ or $info=~/non/){ $NC_name{$chrID}++;
        $line="$chr\tGSAman\t$region\t$start\t$end\t.\t$ori\t.\tID=$chrID;ManInfo=NC;\n"; 
     }elsif($info=~/CDS/){ #print "$line"; #corrdite error in CDS line;
         
     }elsif($info=~/gene/ or $info=~/gee/ or $info=~/GENE/ or $info=~/gne/ or $info=~/ge/){
          $mRNA{$chrID}++;
     }elsif($info=~/delete/){ $Delete{$chrID}++; 

     }else{ print "$line"; }
     push@{$gff{$chrID}{$region}},$line;
  }elsif($line=~/^(chr_\d+)\s+GSAman\s+(mRNA)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);ManInfo=;/){ 
     $chr=$1; $region=$2; $start=$3; $end=$4; $ori=$5; $chrID=$6;
     push@{$gff{$chrID}{$region}},$line;
     if($sequence{$chrID}=~/^ATG.+TGA$/){
         $CDS_name{$chrID}++;
         print CDSfa "$chrID\n$sequence{$chrID}\n";
     }else{ $mRNA{$chrID}++;
     }
  }elsif($line=~/^(chr_\d+)\s+GSAman\s+(exon)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+Parent=(.+?);/){
     $chr=$1; $region=$2; $start=$3; $end=$4; $ori=$5; $chrID=$6;
     push@{$gff{$chrID}{$region}},$line;
  }else{  print "$line"; 

  }
}
open(GFFmRNA,">","$dir/Split_CDS-mRNA-AS-NAT/mRNA_need_coding.gtf");
open(RNA,">","$dir/Split_CDS-mRNA-AS-NAT/mRNA_need_coding.fasta");
foreach $chrID(keys%mRNA){
    print RNA ">$chrID\n$sequence{$chrID}\n";
    @temp=@{$gff{$chrID}{"mRNA"}};
    foreach$line(@temp){   print GFFmRNA "$line"; }
    @temp=sort@{$gff{$chrID}{"exon"}};
    foreach$line(@temp){   print GFFmRNA "$line"; }
    
}
open(AS,">","$dir/Split_CDS-mRNA-AS-NAT/AS.gff3");
open(NC,">","$dir/Split_CDS-mRNA-AS-NAT/NC.gff3");
open(NAT,">","$dir/Split_CDS-mRNA-AS-NAT/NAT.gff3");
foreach $chrID(keys%AS_name){ 
    @temp=@{$gff{$chrID}{"mRNA"}};
    foreach$line(@temp){   print AS "$line"; }
    @temp=sort@{$gff{$chrID}{"exon"}};
    foreach$line(@temp){   print AS "$line"; } 
}
$n=keys%AS_name;
print "AS num:$n\n";
foreach $chrID(keys%NC_name){ 
    @temp=@{$gff{$chrID}{"mRNA"}};
    foreach$line(@temp){   print NC "$line"; }
    @temp=sort@{$gff{$chrID}{"exon"}};
    foreach$line(@temp){   print NC "$line"; }
 }
$n=keys%NC_name;
print "NC num:$n\n";
foreach $chrID(keys%NAT_name){ 
    @temp=@{$gff{$chrID}{"mRNA"}};
    foreach$line(@temp){   print NAT "$line"; }
    @temp=sort@{$gff{$chrID}{"exon"}};
    foreach$line(@temp){   print NAT "$line"; }
 }
$n=keys%NAT_name;
print "NAT num:$n\n";




