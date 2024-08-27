use strict;
my$dir=$ENV{'PWD'};
my(%GTF_name, %GFF_line, $ID, $name);
open(IN1,"<$ARGV[0]") or die; #GTF
open(IN2,"<$ARGV[1]") or die; #GFF

my@GTF=<IN1>; my@GFF=<IN2>; close IN1; close IN2;
foreach my$m(@GTF){ if($m=~/(chr_\d+)\s.+?transcript\s+(\d+)\s+(\d+)\s+1000\s+(.)\s.+?gene_id\s+"(.+?)";\s+transcript_id\s+"(.+?)";.+?FPKM "(.+?)";/){   $GTF_name{$6}=$m; }  }

#chr_001 AUGUSTUS        gene    279     954     .       +       .       ID=STRG.1.1_ORF1;Name=STRG.1.1_ORF1;Note=
#chr_001 AUGUSTUS        mRNA    279     954     .       +       .       ID=STRG.1.1_ORF1;Parent=STRG.1.1_ORF1
#chr_001 AUGUSTUS        five_prime_UTR  279     333     .       +       .       ID=STRG.1.1_ORF1.t1.utr5;Parent=STRG.1.1_ORF1
#chr_001 AUGUSTUS        CDS     334     342     .       +       .       ID=STRG.1.1_ORF1.t1.cds;Parent=STRG.1.1_ORF1
#chr_001 AUGUSTUS        CDS     393     452     .       +       .       ID=STRG.1.1_ORF1.t1.cds;Parent=STRG.1.1_ORF1
#chr_001 AUGUSTUS        CDS     504     854     .       +       .       ID=STRG.1.1_ORF1.t1.cds;Parent=STRG.1.1_ORF1
#chr_001 AUGUSTUS        three_prime_UTR 855     954     .       +       .       ID=STRG.1.1_ORF1.t1.utr3;Parent=STRG.1.1_ORF1i
foreach my$m(@GFF){ if($m=~/gene.+?ID=((.+?)_ORF\d+);Name=(.+?);/){ $GFF_line{$2}{$1}=$m; $ID=$2; $name=$1; 
                   }elsif($m=~/mRNA.+?Parent=((.+?)_ORF\d+)/){ $GFF_line{$2}{$1}.=$m;  
                   }elsif($m=~/CDS.+?Parent=((.+?)_ORF\d+)/){ $GFF_line{$2}{$1}.=$m;
                   }elsif($m=~/UTR.+?Parent=((.+?)_ORF\d+)/){ $GFF_line{$2}{$1}.=$m;     
                   }
}

open(OUT1,">$ARGV[0]_GapNew_gene.gff3");
foreach my$m(keys%GFF_line){
    if(exists$GTF_name{$m}){
        foreach my$n(keys%{$GFF_line{$m}}){
           print OUT1"$GFF_line{$m}{$n}";
        }
     }
}
