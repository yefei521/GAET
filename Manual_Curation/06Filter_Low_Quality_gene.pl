$dir=$ENV{'PWD'};
open(IN3,"<","$dir/Rename_gene/TGD2021blastTGD2024_otf7");
#TTHERM_01528530 chr_137-GSAman06240
foreach(<IN3>){
   if($_=~/^TTH/){
   @temp=split(/\s/,$_);
   next if $temp[0]=~/g\d+/;
   $name=$temp[0];  $Blast{$name}++;
   #print "$name\n";
   }elsif($_=~/#\sQuery:\s+(TTH.+?)\s/){
     $Tname=$1; $name{$Tname}++;
   }
}
$n=keys%Blast; print "$n\n";
$n=keys%name; print "$n\n";
open(OUT1,">","$dir/Rename_gene/LowQuality.txt");
foreach $Tname(keys%name){
   if(!exists$Blast{$Tname}){ print OUT1 "$Tname\n";  }
}
