##This scripts used for analysis the 6mA sites within penetrationaround TSS.
###Usage: perl perl.pl gene.bed m6A.bed well-gene.txt
#perl 01Count_H3K4me3-sigmaP_gene1kb.pl 2-upd-Genome-GFF3-latest-2.gff3_gene H3K4me3-ChIP_140-160.srt.rmdup.bed 2kgene_well-gene-top8500.txt &
my (@ar,@GFF,@bed,$sequence,$chr,$m,$start,$anotation,$oritation,$TSS,$bin_size,%bin_number,%bin_number_sigmap,@keys,$n,%m6A);
open (IN1, "<", "$ARGV[0]") or die "Can't open IN1 : $!";
open (IN2, "<", "$ARGV[1]") or die "Can't open IN2 : $!";
open (IN3, "<", "$ARGV[2]") or die "Can't open IN3 : $!";
open (OUT1, ">", "3end-1kb_H3K4me3-sigmaP");

#chr_001 326     475     A00164:299:H7K5KDSXX:2:1336:14877:1266/2        60      +
foreach $sequence(<IN2>){ if($sequence=~/^(chr_\d+)\s+(\d+)\s+(\d+)\s+.+?\s+\d+\s+(.)\s/){$mid=int(($2+$3)/2); $m6A{"$1\t$mid"}++; $num++; } }
foreach(<IN3>){ if(/^(.+?)\s/){$well_gene{$1}++;} }

#chr_001 AUGUSTUS        gene    226     854     .       +       .       ID=g18410;Name=TTHERM_00161861;Note="hypothetical protein"
foreach $anotation(<IN1>){
     if ($anotation=~/^(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s.+Name=(.+?);/){ 
           next if !exists$well_gene{$5};
           my$len=$3-$2; my$start=$2; my$end=$3; my$name=$5;
           if($4 eq "+"){ $m=$3-1000; if($len>1000){$n=$3;}else{$n=$3;}
                          while($m<=$n){if(exists$m6A{"$1\t$m"}){ $bin_number_sigmap{$name}+=$m6A{"$1\t$m"};    } $m++; }  }
           elsif($4 eq "-"){ $n=$2; if($len>1000){$m=$2+1000;}else{ $m=$2; }
                          while($n<=$m){if(exists$m6A{"$1\t$n"}){ $bin_number_sigmap{$name}+=$m6A{"$1\t$n"};    } $n++; }  }
     }
}
foreach(keys%bin_number_sigmap){
   $ratio=$bin_number_sigmap{$_}/$num;
   print OUT1"$_\t$bin_number_sigmap{$_}\t$ratio\n";

}
