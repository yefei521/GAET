# perl ./5.pl gffcmp.combined.gtf_rep_more2_GT-AG_UTR.gff3 
use strict;
my(%gff, $name, $mark, $gene_id, %Error_gene, @temp, $len5, $len3, $len_cds, $n, @CDS, );
my$dir=$ENV{'PWD'};
open(IN1,"<$ARGV[0]"); 
open(OUT1,">$ARGV[0]_Filtered");
#open(OUT2,">$dir/");

my@GFF=<IN1>; close IN1;
#chr_001 AUGUSTUS        gene    258     1283    .       +       .       ID=STRG.3.1_ORF2;Name=STRG.3.1_ORF2;Note=
#chr_001 AUGUSTUS        mRNA    258     1283    .       +       .       ID=STRG.3.1_ORF2;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        five_prime_UTR  258     333     .       +       .       ID=STRG.3.1_ORF2.t1.utr5;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        CDS     334     342     .       +       .       ID=STRG.3.1_ORF2.t1.cds;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        CDS     393     452     .       +       .       ID=STRG.3.1_ORF2.t1.cds;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        CDS     504     854     .       +       .       ID=STRG.3.1_ORF2.t1.cds;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        three_prime_UTR 855     1283    .       +       .       ID=STRG.3.1_ORF2.t1.utr3;Parent=STRG.3.1_ORF2

foreach(@GFF){ 
       if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5;  $gff{$gene_id}.=$_; 
       }elsif(/chr_\d+.+?CDS.+?Parent=$gene_id/){   $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?five_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){   $gff{$gene_id}.=$_; 
       }elsif(/chr_\d+.+?three_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){   $gff{$gene_id}.=$_; 
       }elsif(/chr_\d+.+?mRNA.+?Parent=$gene_id/){  $gff{$gene_id}.=$_;
       }else{print "$_";}
}
my$n=keys%gff; print "$n\n";
my(%hash, @num,$num, $id, $UTR5_num, $UTR3_num, $UTR5_len, $UTR3_len);
foreach my$name(keys%gff){ my@temp=split(/\n/,$gff{$name}); undef@num;
   foreach my$line(@temp){
       if($line=~/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5; push@num,$line; $num=@num+1; $id="t$num";  $hash{$gene_id}{$id}="$line\n"; 
       }elsif($line=~/chr_\d+.+?CDS.+?Parent=$gene_id/){   $hash{$gene_id}{$id}.="$line\n";
       }elsif($line=~/chr_\d+.+?five_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){   $hash{$gene_id}{$id}.="$line\n"; 
       }elsif($line=~/chr_\d+.+?three_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){   $hash{$gene_id}{$id}.="$line\n"; 
       }elsif($line=~/chr_\d+.+?mRNA.+?Parent=$gene_id/){  $hash{$gene_id}{$id}.="$line\n"; }        
  }
}

foreach my$name(keys%hash){ 
    foreach my$id(keys%{$hash{$name}}){ my@temp=split(/\n/,$hash{$name}{$id}); 
       foreach my$line(@temp){
          if($line=~/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5;  $UTR5_num=0; $UTR3_num=0; $UTR5_len=0; $UTR3_len=0; 
           }elsif($line=~/chr_\d+.+?CDS.+?Parent=$gene_id/){   
           }elsif($line=~/chr_\d+.+?five_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){    $UTR5_num++; $UTR5_len+=$2-$1+1;
           }elsif($line=~/chr_\d+.+?three_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){    $UTR3_num++; $UTR3_len+=$2-$1+1;
           }elsif($line=~/chr_\d+.+?mRNA.+?Parent=$gene_id/){   }
       } #print "$UTR3_num\t$UTR3_len\n$UTR5_num\t$UTR5_len\n\n\n"; 
       if($UTR3_num>2 || $UTR5_num>2 || $UTR3_len>1000 || $UTR5_len>1000){ delete$hash{$name}{$id}; print "$UTR3_num\t$UTR3_len\n$UTR5_num\t$UTR5_len\n\n\n";  
       }else{  $hash{$name}{$id}=~s/t1/$id/g;  }
    }
}


foreach my$name(keys%hash){
    foreach my$id(sort{$a<=>$b}keys%{$hash{$name}}){
         print OUT1 "$hash{$name}{$id}";
    }
}
