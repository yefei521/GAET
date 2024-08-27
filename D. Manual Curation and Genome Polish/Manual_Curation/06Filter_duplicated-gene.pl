$dir=$ENV{'PWD'};
open(IN1,"<","$dir/Rename_gene/TGD2024_need_rename_gene.gff3");

@IN1=<IN1>;close IN1;
foreach(@IN1){
    if(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1;.+?=(.+?)\s/){
        $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5;    $ori2024{$name}=$ori;
        $mid_point_2024{$chr}{$name}=int(($start+$end)/2); push@{$gff2024{$chr}{$name}},$_;
        if($name=~/(chr_\d+)-GSAman(\d+)/ && length($2)==8){ $New_Tname{$name}="TTHERM_$2"; $Tname="TTHERM_$2"; $Used_Tname{$Tname}++; }
    }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s.+?ID=(.+?);Name=(.+?.t1);Note=/){
        $name=$4; push@{$gff2024{$1}{$name}},$_;
    }elsif(/(chr_\d+)\s.+?five_prime_UTR\s+(\d+)\s+(\d+)\s.+?ID=(.+?).t1.utr;Parent=(.+?).t1\s/){
        $name=$4; push@{$gff2024{$1}{$name}},$_;
    }elsif(/(chr_\d+)\s.+?CDS\s+(\d+)\s+(\d+)\s.+?ID=(.+?).t1.cds;Parent=(.+?).t1\s/){
        $name=$4; push@{$gff2024{$1}{$name}},$_; push@{$CDS2024{$name}},$2,$3;
    }elsif(/(chr_\d+)\s.+?three_prime_UTR\s+(\d+)\s+(\d+)\s.+?ID=(.+?).t1.utr;Parent=(.+?).t1\s/){
        $name=$4; push@{$gff2024{$1}{$name}},$_;
    }else{print "$_";}
}

foreach(@IN1){
   if(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1;.+?=(.+?)\s/){
      $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5;
      foreach $name2(keys%{$mid_point_2024{$chr}}){$mid=$mid_point_2024{$chr}{$name2}; if($mid>=$start && $mid<=$end){  push@{$check_2024_dup{$name2}},$name; } }
   }
}

open(OUT4,">","$dir/Rename_gene/TGD2024_duplicated_nedd_check.txt");
foreach $name(sort keys%check_2024_dup){
    @temp=@{$check_2024_dup{$name}};
    if(@temp>1){$line=join("|",@temp); print OUT4 "query:$name\t$line\n"; }
}


