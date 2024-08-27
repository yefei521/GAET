use strict;
my$dir=$ENV{'PWD'};
my(@GFF, @GTF, %Name, $name, %corrd, $chr, $ID, $min_sdv, $mark, %NEW, $sum, $ave, $m, %GTF, %sdv, );
open(IN1,"<$ARGV[0]"); @GFF=<IN1>; close IN1; # SB210_DRS_DNA_Q7_maped_sorted.gtf_UTR.gff3_Filtered_sorted
open(IN2,"<$ARGV[1]"); @GTF=<IN2>; close IN2; #SB210_DRS_DNA_Q7_maped_sorted.gtf

#chr_001 AUGUSTUS        gene    258     1283    .       +       .       ID=STRG.3.1_ORF2;Name=STRG.3.1_ORF2;Note=
#chr_001 AUGUSTUS        mRNA    258     1283    .       +       .       ID=STRG.3.1_ORF2;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        five_prime_UTR  258     333     .       +       .       ID=STRG.3.1_ORF2.t1.utr5;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        CDS     334     342     .       +       .       ID=STRG.3.1_ORF2.t1.cds;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        CDS     393     452     .       +       .       ID=STRG.3.1_ORF2.t1.cds;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        CDS     504     854     .       +       .       ID=STRG.3.1_ORF2.t1.cds;Parent=STRG.3.1_ORF2
#chr_001 AUGUSTUS        three_prime_UTR 855     1283    .       +       .       ID=STRG.3.1_ORF2.t1.utr3;Parent=STRG.3.1_ORF2
foreach my$line(@GFF){ if($line=~/gene\t(\d+)\t(\d+)\t.+?ID=(((STRG\.\d+)\.\d+)_ORF\d+);/){ $Name{$4}=$5; }  }
my$x=keys%Name; print"$x\n";
foreach my$line(@GFF){ if($line=~/^(chr_\d+)\s.+?gene\t(\d+)\t(\d+)\t.+?ID=(.+?)_ORF\d+;/){        $name=$4; $corrd{$name}=$line; $chr=$1;
                           }elsif($line=~/$chr.+?mRNA\t(\d+)\t(\d+)\t.+?Parent=(.+)\n/){  $corrd{$name}.=$line;
                           }elsif($line=~/$chr.+?CDS\t(\d+)\t(\d+)\t.+?Parent=(.+)\n/){    $corrd{$name}.=$line;
                           }elsif($line=~/$chr.+?UTR\t(\d+)\t(\d+)\t.+?Parent=(.+)\n/){     $corrd{$name}.=$line; }
}
my$x=keys%corrd; print"$x\n";
foreach my$line(@GTF){ if($line=~/transcript\t.+?transcript_id "(.+?)";/){$name=$1; }elsif($line=~/exon\t.+?transcript_id "($name)";\s.+?cov\s+"(.+?)";/){$GTF{$name}.="$2\n";  }  }
my$x=keys%GTF; print"$x\n";
foreach my$line(keys%GTF){ my@temp=split(/\n/,$GTF{$line});
                           foreach(@temp){$sum+=$_; } $ave=$sum/@temp;
                           foreach(@temp){$m+=($_-$ave)*($_-$ave) }
                           $sdv{$line}=$m/@temp;
}
my$x=keys%sdv; print"$x\n";
foreach my$line(sort{$a<=>$b}keys%corrd){ 
        $line=~/^((STRG\.\d+)\.\d+)/; $ID=$2; $min_sdv="none";
        foreach my$line2(sort{$a<=>$b}keys%Name){
                if($Name{$line2} eq $ID){ 
                        if($min_sdv eq "none"){ $min_sdv=$sdv{$line2}; $mark=$line2; delete$Name{$line2};
                        }else{ if($sdv{$line2}<$min_sdv){ $min_sdv=$sdv{$line2}; $mark=$line2; }  delete$Name{$line2};    }
                } 
        }
        $NEW{$ID}=$mark;  
}
my$x=keys%NEW; print"$x\n";
open(OUT1,">$ARGV[0]_RmAS");
open(OUT2,">$ARGV[0]_RmAS.name");
foreach my$line(sort{$a<=>$b}keys%NEW){  print OUT1"$corrd{$NEW{$line}}"; print OUT2"$NEW{$line}\n"; }


#chr_001 StringTie       transcript      1457545 1458860 1000    -       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "16.000162"; FPKM "2.624581"; TPM "0.431425";
#chr_001 StringTie       exon    1457545 1458673 1000    -       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "17.149868";
#chr_001 StringTie       exon    1458754 1458860 1000    -       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2"; cov "3.869159";




