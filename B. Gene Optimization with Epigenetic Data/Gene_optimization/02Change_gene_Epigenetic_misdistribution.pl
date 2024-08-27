use strict;
my$dir=$ENV{'PWD'};
my(%gff, @gff, $chr, $gene_id, %gff_6mA, $gene_line, $gene_start, $gene_end, @m6A_corrd, @m6A_corrd_RE, $number, $length, $mean_len, $mark, $n1, $n2, $n3, $n4);
#------------------------------------
open(IN1,"<$dir/06Merge_Deleted636-81_UTR2901-550_To13gff3.gff3");
print"Step1: Begin import gene information from 06Merge_Deleted636-81_UTR2901-550_To13gff3.gff3\nPlease wait .......\n";
foreach(<IN1>){
     if(/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+.+?ID=(chr_\d+);/){ $gff{$3}=$_; push@gff,$3;
       }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){  $gene_id=$5; $gff{$gene_id}.=$_; $chr=$1; push@gff,$5;
       }elsif(/$chr.+?mRNA.+?Parent=$gene_id/){  $gff{$gene_id}.=$_;
       }elsif(/$chr.+?UTR.+?Parent=$gene_id/){   $gff{$gene_id}.=$_;
       }elsif(/$chr.+?CDS.+?Parent=$gene_id/){   $gff{$gene_id}.=$_;
       }elsif(/chr.+?rDNA.+?ID=(.+?);Name/){     $gff{$1}=$_; push@gff,$1;
       }
}

#---------------------------------------
open(IN2,"<$dir/06gff3_m6A_unregular_ingene.txt");
print"STep2: Begin import unregular distribution 6mA information within gene from 06gff3_m6A_unregular_ingene.txt";
foreach(<IN2>){
   if(/^(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){ $gene_id=$5; $gff_6mA{$gene_id}.=$_; $chr=$1;
  }elsif(/^$chr.+?group/){$gff_6mA{$gene_id}.=$_;
  }elsif(/^\d+/){
  }else{ #print"$_";
  }
}
#--------------------------------------
open(OUT1,">$dir/02gff3_gene_with_6mA.gff");
print"Step3: Begin output gene corrd information and 6mA distribution corrd in 02gff3_gene_with_6mA.gff\n Please wait .......\n";
foreach(%gff_6mA){ print OUT1"$gff{$_}$gff_6mA{$_}\n"; }

#-------------------------------------
open(IN3,"<$dir/02gff3_gene_with_6mA.gff");
open(OUT2,">$dir/02New_gene_accoding_6mA_distribution.txt");
open(OUT3,">$dir/02gene_corrd_need_change.txt");
open(OUT4,">$dir/02gene_oritation_nedd_check.txt");
open(OUT5,">$dir/02gene_need_check.txt");
#chr_070 gene12948 group1: 11757|-3778    11760|-3781
$n1=$n2=$n3=$n4=0;
print"Step4: Begin seareach gene by relative corrd of 6mA group\n Please wait ........\n";
foreach(<IN3>){
    if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){ $gene_id=$5;  $chr=$1; $gene_start=$2; $gene_end=$3; $gene_line=$_;
   }elsif(/$chr.+?group/){ #print "1\t"; 
       @m6A_corrd=split(/\s+/,$_); my$n=@m6A_corrd; undef@m6A_corrd_RE;
       foreach my$i(@m6A_corrd){ if($i=~/(\d+)\|(.+)/){ if($2!=$mark && abs($2-$m6A_corrd_RE[-1])<1200 ){push@m6A_corrd_RE,$2; $mark=$2; } $mark=$m6A_corrd_RE[-1]; }  }
       $number=0; $length=0; $mean_len=0; 
       if(@m6A_corrd_RE>0){
         foreach my$j(@m6A_corrd_RE){ $length+=$j; $number++; }  $mean_len=$length/$number;
         if($mean_len<-1500){                         print OUT2"$gff{$gene_id}$_\n";    $n1++;
           }elsif($mean_len<1500 && $mean_len>-1500){ print OUT3"$_$gff{$gene_id}\n"; $n2++;
           }elsif($mean_len>2500){                    print OUT4"$gff{$gene_id}$_\n";     $n3++;
           }else{                                     print OUT5"$gff{$gene_id}$_\n";     $n4++;
           }
       }
   }elsif(/^\d+/){
   }else{ #print"$_";
   }
}
print"New_gene_accoding_6mA_distribution: $n1\n";
print"gene_corrd_need_change: $n2\n";
print"gene_oritation_nedd_check: $n3\n";
print"gene_need_check: $n4\n";
