use strict;
my(@GFF, @supercontig, @chr, %corrd, $name, %UTR_num);
my$dir=$ENV{'PWD'};
open(IN1,"<$ARGV[0]"); #SB210_DRS_DNA_Q7_maped_sorted.gtf_UTR.gff3_Filtered_sorted_RmAS
open(IN2,"<$dir/supercontig.gff3");
open(OUT1,">$ARGV[0]_sorted");
@GFF=<IN1>; close IN1;
@supercontig=<IN2>; close IN2;
foreach(@supercontig){ if(/^(chr_\d+)\s/){push@chr,$1; } }

#chr_114 AUGUSTUS        gene    1237989 1242412 .       +       .       ID=STRG.12599.1_ORF8;Name=STRG.12599.1_ORF8;Note=
#chr_001 AUGUSTUS        supercontig     1       1459056 .       .       .       ID=chr_001;Name=chr_001
foreach my$chr(@chr){
    undef%corrd; 
    foreach my$supercontig(@supercontig){ if($supercontig=~/$chr.+?supercontig\t(\d+)\t(\d+)\t.+?ID=(.+?);Name/){$corrd{"$1\t$2\t$3"}=$supercontig; } }
    foreach my$line(@GFF){ if($line=~/^$chr.+?gene\t(\d+)\t(\d+)\t.+?ID=(.+?);/){        $name="$1\t$2\t$3"; $corrd{$name}=$line; 
                           }elsif($line=~/$chr.+?mRNA\t(\d+)\t(\d+)\t.+?Parent=(.+)\n/){  $corrd{$name}.=$line; 
                           }elsif($line=~/$chr.+?CDS\t(\d+)\t(\d+)\t.+?Parent=(.+)\n/){    $corrd{$name}.=$line;
                           }elsif($line=~/$chr.+?UTR\t(\d+)\t(\d+)\t.+?Parent=(.+)\n/){     $corrd{$name}.=$line; $UTR_num{$name}++;
                           } 
    } 
    foreach(keys%UTR_num){ if($UTR_num{$_}>=4){delete$corrd{$_};   } }
    foreach(sort{$a<=>$b}keys%corrd){print OUT1"$corrd{$_}"; }   
}
