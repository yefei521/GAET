
$dir=$ENV{'PWD'};
open(GENE,"<","$dir/Rename_gene/TGD2024_renamed.gff3");
#chr_001 Manual  gene    256     1182    .       +       .       ID=g18410;Name=TTHERM_00161861;UserNamed=Unnamed;TGDAnno=Unknow;info=hypothetical protein;
foreach(<GENE>){
   if(/^(chr_\d+)\s+Manual\s+gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);Name=(.+?);/){
     $chr=$1; $start=$2; $end=$3; $ori=$4; $gname=$5; $Tname=$6; $Tname2chr{$Tname}=$chr;
     push@{$sorted_name{$chr}},$Tname; $Tname2gname{$Tname}=$gname;
     $mid_gene{$Tname}=int(($start+$end)/2); $start_gene{$Tname}=$start; $end_gene{$Tname}=$end;
     push@{$gff{$chr}{$Tname}{"gene"}},$_; 
   }elsif(/^(chr_\d+)\s+Manual\s+mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);Par/){
     push@{$gff{$chr}{$Tname}{"mRNA"}},$_;
   }elsif(/^(chr_\d+)\s+Manual\s+(.+?)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);Par/){
     $region=$2; push@{$gff{$chr}{$Tname}{"exon"}},$_;
   }
}
open(AS,"<","$dir/Anno_AS_NAT_NC/AS.gff3");
open(NAT,"<","$dir/Anno_AS_NAT_NC/NAT.gff3");
open(NC,"<","$dir/Anno_AS_NAT_NC/NC.gff3");
foreach(<AS>){
   if(/(chr_\d+)\tGSAman\t(mRNA)\t(\d+)\t(\d+)\t.\t(.)\t.\tID=(.+?);ManInfo=(.+?);\n/){
      $chr=$1; $start=$3; $end=$4; $ori=$5; $name=$6; $AS2chr{$Tname}=$chr;
      push@{$ASgff{$chr}{$name}{"mRNA"}},$_;  
      $mid_AS{$name}=int(($start+$end)/2); $start_AS{$name}=$start; $end_AS{$name}=$end;
   }elsif(/(chr_\d+)\tGSAman\t(exon)\t(\d+)\t(\d+)\t.\t(.)\t.\tParent=(.+?);\n/){
      push@{$ASgff{$1}{$6}{"exon"}},$_;
   }
}

foreach(<NAT>){
   if(/(chr_\d+)\tGSAman\t(mRNA)\t(\d+)\t(\d+)\t.\t(.)\t.\tID=(.+?);ManInfo=(.+?);\n/){
      $chr=$1; $start=$3; $end=$4; $ori=$5; $name=$6; $NAT2chr{$Tname}=$chr;
      push@{$NATgff{$chr}{$name}{"mRNA"}},$_;
      $mid_NAT{$name}=int(($start+$end)/2); $start_NAT{$name}=$start; $end_NAT{$name}=$end;
   }elsif(/(chr_\d+)\tGSAman\t(exon)\t(\d+)\t(\d+)\t.\t(.)\t.\tParent=(.+?);\n/){
      push@{$NATgff{$1}{$6}{"exon"}},$_;
   }
}

foreach(<NC>){
   if(/(chr_\d+)\tGSAman\t(mRNA)\t(\d+)\t(\d+)\t.\t(.)\t.\tID=(.+?);ManInfo=(.+?);\n/){
      $chr=$1; $start=$3; $end=$4; $ori=$5; $name=$6; $NAT2chr{$Tname}=$chr;
      push@{$NCgff{$chr}{$name}{"mRNA"}},$_;
      $mid_NC{$name}=int(($start+$end)/2); $start_NC{$name}=$start; $end_NC{$name}=$end;
   }elsif(/(chr_\d+)\tGSAman\t(exon)\t(\d+)\t(\d+)\t.\t(.)\t.\tParent=(.+?);\n/){
      push@{$NCgff{$1}{$6}{"exon"}},$_;
   }
}

open(IN1,"<","$dir/Anno_AS_NAT_NC/ASblastTGD2024_otf7");
open(IN2,"<","$dir/Anno_AS_NAT_NC/NATblastTGD2024_otf7");
open(IN3,"<","$dir/Anno_AS_NAT_NC/NCblastTGD2024_otf7");
#chr_020-GSAman05925     TTHERM_00704030.t1
foreach(<IN1>){
   next if $_=~/^#/;
   @temp=split(/\s/,$_);
   next if $temp[1]=~/g\d+/;
   $temp[0]=~/((chr_\d+)-.+)/; $name=$1; $name_chr=$2;
   $start_AS=$start_AS{$name}; $mid_AS=$mid_AS{$name}; $end_AS=$end_AS{$name};
    $Tname=$temp[1];
   $start_Tname=$start_gene{$Tname}; $mid_Tname=$mid_gene{$Tname}; $end_Tname=$end_gene{$Tname}; 
   $Tname_chr=$Tname2chr{$Tname};
   if(!exists$AS{$name} && $name_chr eq $Tname_chr){ $mark=0;
      if($mid_AS>$start_Tname && $mid_AS<$end_Tname){ $mark=1;  }
      if($mid_Tname>$start_AS && $mid_Tname<$end_AS){ $mark=1;  }
      if($mark==1){ push@{$blast_AST2024{$Tname}},$name; $AS{$name}++; 
      }
    }
}
$n=keys%AS;
print "AS num:$n\n";
foreach(<IN2>){
   next if $_=~/^#/;
   @temp=split(/\s/,$_);
   next if $temp[1]=~/g\d+/;
   $temp[0]=~/((chr_\d+)-.+)/; $name=$1; $name_chr=$2;
   $start_NAT=$start_NAT{$name}; $mid_NAT=$mid_NAT{$name}; $end_NAT=$end_NAT{$name};
   $Tname=$temp[1];
   $start_Tname=$start_gene{$Tname}; $mid_Tname=$mid_gene{$Tname}; $end_Tname=$end_gene{$Tname}; 
   $Tname_chr=$Tname2chr{$Tname};
   if(!exists$NAT{$name} && $name_chr eq $Tname_chr){ $mark=0;
       if($mid_NAT>$start_Tname && $mid_NAT<$end_Tname){ $mark=1;  }
       if($mid_Tname>$start_NAT && $mid_Tname<$end_NAT){ $mark=1;  }
       if(($start_NAT>$start_Tname && $start_NAT<$end_Tname) or ($end_NAT>$start_Tname && $end_NAT<$end_Tname) or ($start_Tname>$end_NAT && $start_Tname-$end_NAT<300) or ($start_NAT>$end_Tname && $start_NAT-$end_Tname<300)){ $mark=1; }
       if($mark==1){  push@{$blast_NATT2024{$Tname}},$name;  $NAT{$name}++;  }
   }  
}
$n=keys%NAT;
print "NAT num:$n\n";
foreach(<IN3>){
   next if $_=~/^#/;
   @temp=split(/\s/,$_);
   next if $temp[1]=~/g\d+/;
   $temp[0]=~/((chr_\d+)-.+)/; $name=$1; $name_chr=$2;
   $Tname=$temp[1];
   $Tname_chr=$Tname2chr{$Tname};
   if(!exists$NAT{$name} && $name_chr eq $Tname_chr){ $mark=0;
       if($mid_NAT>$start_Tname && $mid_NAT<$end_Tname){ $mark=1;  }
       if($mid_Tname>$start_NAT && $mid_Tname<$end_NAT){ $mark=1;  }
       if($mark==1){ push@{$blast_NCT2024{$Tname}},$name;  $NC{$name}++;  }
   }
}
$n=keys%NC;
print "NC num:$n\n";
open(OUT1,">","$dir/Anno_AS_NAT_NC/TGD2024_renamed_addAS-NAT-NC.gff3");
#if(exists$NCgff{"chr_066"}{"chr_066-GSAman0079670"}){@temp=@{$NCgff{"chr_066"}{"chr_066-GSAman0079670"}{"exon"}}; print "@temp"; }
foreach $chr(sort{$a cmp $b}keys%sorted_name){
   @gene=@{$sorted_name{$chr}}; $i=0;
   while($i<=$#gene){
      $Tname=$gene[$i]; $gname=$Tname2gname{$Tname};
      @temp=@{$gff{$chr}{$Tname}{"gene"}};  foreach $l(@temp){ print OUT1 "$l"; }
      @temp=@{$gff{$chr}{$Tname}{"mRNA"}};  foreach $l(@temp){ print OUT1 "$l"; }
      @temp=@{$gff{$chr}{$Tname}{"exon"}};  foreach $l(@temp){ print OUT1 "$l"; }
      @AS=@{$blast_AST2024{$Tname}};
      $j=2;
      foreach$AS_name(@AS){
         next if exists$OUT_AS{$AS_name};
         $OUT_AS{$AS_name}++;
         @temp=@{$ASgff{$chr}{$AS_name}{"mRNA"}}; 
         foreach $l(@temp){ if($l=~/^(.+?)\tID=(.+?);(ManInfo=.+?);\n/){ print OUT1 "$1\tID=$gname.t$j;Parent=$gname;$3;\n"; } }
         @temp=@{$ASgff{$chr}{$AS_name}{"exon"}};
         foreach $l(@temp){ if($l=~/^(.+?)\tParent=(.+?);\n/){ print OUT1 "$1\tID=$gname.t$j.exon;Parent=$gname.t$j;\n"; } }
         $j++;
      }
      @NAT=@{$blast_NATT2024{$Tname}};
      foreach$NAT_name(@NAT){
         next if exists$OUT_NAT{$NAT_name};
         @temp=@{$NATgff{$chr}{$NAT_name}{"mRNA"}}; 
         foreach $l(@temp){ if($l=~/^(.+?)\tID=(.+?);(ManInfo=.+?);\n/){ print OUT1 "$1\tID=$gname.t$j;Parent=$gname;$3;\n"; } }
         @temp=@{$NATgff{$chr}{$NAT_name}{"exon"}};
         foreach $l(@temp){ if($l=~/^(.+?)\tParent=(.+?);\n/){ print OUT1 "$1\tID=$gname.t$j.exon;Parent=$gname.t$j;\n"; }  }
         $j++;
         $OUT_NAT{$NAT_name}++;
      }
      @NC=@{$blast_NCT2024{$Tname}};
      foreach$NC_name(@NC){ 
         next if exists$OUT_NC{$NC_name};
         $OUT_NC{$NC_name}++;
         #@temp=@{$NCgff{$chr}{$NC_name}{"mRNA"}}; foreach $l(@temp){ print OUT1 "$l"; }
         #@temp=@{$NCgff{$chr}{$NC_name}{"exon"}}; foreach $l(@temp){ print OUT1 "$l"; }
         @temp=@{$NCgff{$chr}{$NC_name}{"mRNA"}};
         foreach $l(@temp){ if($l=~/^(.+?)\tID=(.+?);(ManInfo=.+?);\n/){ print OUT1 "$1\tID=$2.t1;Parent=$2;$3;\n"; } }
         @temp=@{$NCgff{$chr}{$NC_name}{"exon"}};
         foreach $l(@temp){ if($l=~/^(.+?)\tParent=(.+?);\n/){ print OUT1 "$1\tID=$2.t1.exon;Parent=$2.t1;\n"; }  }
      }
      $i++;
   }
   foreach $NC_name(keys%{$NCgff{$chr}}){  
         next if exists$OUT_NC{$NC_name};
         @temp=@{$NCgff{$chr}{$NC_name}{"mRNA"}}; 
         foreach $l(@temp){ if($l=~/^(.+?)\tID=(.+?);(ManInfo=.+?);\n/){ print OUT1 "$1\tID=$2.t1;Parent=$2;$3;\n"; } }
         @temp=@{$NCgff{$chr}{$NC_name}{"exon"}}; 
         foreach $l(@temp){ if($l=~/^(.+?)\tParent=(.+?);\n/){ print OUT1 "$1\tID=$2.t1.exon;Parent=$2.t1;\n"; }  }
         $OUT_NC{$NC_name}++; 
  }
}

