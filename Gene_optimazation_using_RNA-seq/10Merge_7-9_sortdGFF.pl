use strict;
my(@GFF, @supercontig, @chr, %corrd, $name, );
my$dir=$ENV{'PWD'};
open(IN1,"<$ARGV[0]"); #7Merged_RMoverlap_gene.gff3_sorted
open(IN2,"<$dir/supercontig.gff3");
open(IN3,"<$ARGV[1]"); #9Merged_RMoverlap_newgene.gff3_sorted
open(OUT1,">$ARGV[0]_$ARGV[1]_sorted");
my@GFF7=<IN1>; close IN1;
my@GFF9=<IN3>; close IN3;
push@GFF,@GFF7,@GFF9;
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
                           }elsif($line=~/$chr.+?UTR\t(\d+)\t(\d+)\t.+?Parent=(.+)\n/){     $corrd{$name}.=$line;
                           } 
    } 
    foreach(sort{$a<=>$b}keys%corrd){print OUT1"$corrd{$_}"; }   
}
