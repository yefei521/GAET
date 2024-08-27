$dir=$ENV{'PWD'};
#open(IN1,"<","$dir/");
#open(IN2,"<","$dir/");
my@file_raw_gff3=glob("$dir/Manual-Check_RawGFF/*gff3");
foreach $file(@file_raw_gff3){
    open(IN,"<","$file") or die;
    foreach $line(<IN>){ #print "$line"; 
     #chr_001 GSAman  mRNA    255     1182    .       +       .       ID=GSAman0003;Parent=GSAman0003.gene
     #chr_001 GSAman  exon    286     342     .       +       .       Parent=GSAman00161810CDS
       if($line=~/^(chr_\d+)\s+GSAman\s+(mRNA)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);Parent=(.+?);ManInfo=(.+?)\s/){
             $chr=$1; $region=$2; $start=$3; $end=$4; $ori=$5; $ID=$6; $Parent=$7; $info=$8; $chrID="$chr-$ID";
             $line1="$chr\tGSAman\t$region\t$start\t$end\t.\t$ori\t.\tID=$chrID;";
             push@{$Transcripts{$chrID}},$line1;  if($ID=~/(GSAman\d+)(\D+)/){  $chrID2Info{$chrID}=$2;  }else{  }
       }elsif($line=~/(chr_\d+)\s+GSAman\s+(mRNA)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);Parent=(.+?)\s/){
         #chr_001 GSAman  mRNA    258     1182    .       +       .       ID=GSAman00161861AS2IR;Parent=GSAman00161861AS2IR.gene
             $chr=$1; $region=$2; $start=$3; $end=$4; $ori=$5; $ID=$6; $Parent=$7;  $chrID="$chr-$ID";
             $line1="$chr\tGSAman\t$region\t$start\t$end\t.\t$ori\t.\tID=$chrID;";
             push@{$Transcripts{$chrID}},$line1; if($ID=~/(GSAman\d+)(\D+)/){ $chrID2Info{$chrID}=$2;  }else{ }
       }elsif($line=~/^(chr_\d+)\s+GSAman\s+(mRNA)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);ManInfo=(.+?)\s/){
             $chr=$1; $region=$2; $start=$3; $end=$4; $ori=$5; $ID=$6; $info=$7;  $chrID="$chr-$ID";
             $line1="$chr\tGSAman\t$region\t$start\t$end\t.\t$ori\t.\tID=$chrID;";
             push@{$Transcripts{$chrID}},$line1;   $chrID2Info{$chrID}=$info;
       }elsif($line=~/^(chr_\d+)\s+GSAman\s+(exon)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+Parent=(.+?)\s/){
             $chr=$1; $region=$2; $start=$3; $end=$4; $ori=$5; $ID=$6; $chrID="$chr-$ID";
             $line1="$chr\tGSAman\t$region\t$start\t$end\t.\t$ori\t.\tParent=$chrID;";
             push@{$Transcripts{$chrID}},$line1;
       }elsif($line=~/^(chr_\d+)\s+GSAman\s+(CDS)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+Parent=(.+?)\s/){
             $chr=$1; $region=$2; $start=$3; $end=$4; $ori=$5; $ID=$6; $chrID="$chr-$ID";
             $line1="$chr\tGSAman\t$region\t$start\t$end\t.\t$ori\t.\tParent=$chrID;";
             push@{$Transcripts{$chrID}},$line1;
       }elsif($line=~/^(chr_\d+)\s+GSAman\s+(mRNA)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?)\s/){
             $chr=$1; $region=$2; $start=$3; $end=$4; $ori=$5; $ID=$6; $chrID="$chr-$ID";
             $line1="$chr\tGSAman\t$region\t$start\t$end\t.\t$ori\t.\tID=$chrID;";
             push@{$Transcripts{$chrID}},$line1; if($ID=~/(GSAman\d+)(\D+)/){  $chrID2Info{$chrID}=$2;  }else{ }
       }else{ #print "$line";

       }
    }
}

foreach $chrID(keys%Transcripts){
   @temp=sort{ $a <=> $b }@{$Transcripts{$chrID}};
   foreach $line(@temp){
      if($line=~/ID=(.+?);/ && $1 eq $chrID){ $info=$chrID2Info{$chrID};  $line.="ManInfo=$info;"; print "$line\n";
      }elsif($line=~/Parent=(.+?);/ && $1 eq $chrID){  print "$line\n" 
      }else{ #print "$chrID\n";
      }
      #print "@temp1";
   }
}
