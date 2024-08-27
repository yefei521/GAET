# perl ./5.pl SB210_DRS_DNA_Q7_maped_sorted.gtf_UTR.gff3 SB210_DRS_DNA_Q7_maped_sorted.gtf
use strict;
my($name, $mark, %GFF, %Error_gene, @temp, $len5, $len3, $len_cds, $n, @CDS, );
my$dir=$ENV{'PWD'};
open(IN1,"<$ARGV[0]"); 
open(IN2,"<$ARGV[1]");
open(IN3,"<$dir/02ATG2TGA_gene_result.txt");
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
foreach my$line(@GFF){
   if($line=~/gene.+?ID=((.+?)_ORF\d+);Name/){ $name=$2; $mark=$1; $GFF{$name}.=$line;
   }elsif($line=~/Parent=$mark/){$GFF{$name}.=$line; }
}
my$m=keys%GFF; print"$m\n";

##STRG.14656.1_ORF3       ATT-GAT no
my@Error_gene=<IN3>; close IN3;
foreach my$line(@Error_gene){ if($line=~/^((.+?)_ORF\d+)\t.+?\t(no)/){$Error_gene{$2}=$1; } }

my$m=keys%Error_gene; print"$m\n";
#chr_001 StringTie       transcript      1457907 1458914 1000    +       .       gene_id "STRG.2"; transcript_id "STRG.2.1"; cov "11.518440"; FPKM "1.889423"; TPM "0.310581";
#chr_001 StringTie       exon    1457907 1458352 1000    +       .       gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; cov "2.589686";
#chr_001 StringTie       exon    1458412 1458914 1000    +       .       gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "2"; cov "19.435389";
my@GTF=<IN2>; close IN2;
foreach my$line(@GTF){ if($line=~/transcript_id "(.+?)"; cov "(.+?)"; FPKM "(.+?)";/){ if($2<10 || exists$Error_gene{$1}){ delete$GFF{$1}; } } } ##merge.gtf cov set=10
my$m=keys%GFF; print"$m\n";

foreach my$line(keys%GFF){
    @temp=split(/\n/,$GFF{$line}); $len5=0; $len3=0; undef@CDS; $mark=0;
    foreach my$line2(@temp){
       if($line2=~/gene.+?ID=((.+?)_ORF\d+);Name/){  
       }elsif($line2=~/five_prime_UTR\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+)/){$len5+=$2-$1+1; 
       }elsif($line2=~/CDS\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+)/){ push@CDS,$1,$2;  
       }elsif( $line2=~/three_prime_UTR\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+)/){$len3+=$2-$1+1; }
    }
    pop@CDS; shift@CDS; $n=0; while($n<$#CDS){ $len_cds=$CDS[$n+1]-$CDS[$n]-1; if($len_cds>4000){$mark=1; last; }  $n+=2;  }
    if( $mark==1 || $len5>3000 || $len3>3000 ){next; }
    print OUT1 "$GFF{$line}";
}

