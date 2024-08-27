#perl 01Select_ATG2TGA_cds.pl total_transcripts_exon.fasta_orf-1
use strict;
my(@cds,@list,$n,$name,%list);
open (IN1, "<", "$ARGV[0]")
                or die "Can't open IN1 : $!";
open (OUT1, ">", "$ARGV[0]_ATG2TGA");
open (OUT2, ">", "$ARGV[0]_ATG2TGA_name");
{@cds=<IN1>;close IN1;}
$n=0;
#>lcl|STRG.3.1:760-924 ORF1_STRG.3.1:759:923
my(%fasta,%cds_ATG2TGA);
foreach(@cds){
   if(/^>.+?(ORF\d+)_(STRG.\d+\.\d+):(\d+):(\d+)\n/){
     my$len=abs($4-$3)+1;
     $name=">$2\t$1\t$3\t$4\t$len"; $n=1; 
#   }elsif(/^>MSTRG.\d+.\d+:c\d+-\d+\s+ORF\d+_(MSTRG.\d+.\d+):(\d+):(\d+)\n/){
#     $n=0;
   
   }elsif($n==1){
     chomp;
     $fasta{$name}.=$_;  }
}
#$n=keys %fasta; print "$n\n";
#print "$fasta{$name}\n";
foreach(keys %fasta){
     if($fasta{$_}=~/^ATG.+?TGA$/){
       print OUT1 "$_\n$fasta{$_}\n";
       my@a=split />/;
       print OUT2 "$a[1]\n";}
}
close OUT1;
close OUT2;
my(%max, %longest);
my$dir=$ENV{'PWD'};
open(IN2,"<$dir/$ARGV[0]_ATG2TGA_name");
#open(IN3,"<$dir/11$ARGV[0]_ATG2TGA");
open(OUT3,">$dir/$ARGV[0]_ATG2TGA_longest");
foreach my$i(<IN2>){chomp$i;
    if($i=~/(.+?)\t(ORF\d+)\t(\d+)\t(\d+)\t(\d+)/){
        if($5>$max{$1}){ $longest{$1}=$i; $max{$1}=$5; }
    }
}
foreach(keys%longest){chomp;
    if($longest{$_}=~/(.+?)\t(ORF\d+)\t(\d+)\t(\d+)\t(\d+)/){
        my$name=">$longest{$_}";
        print OUT3"$name\n$fasta{$name}\n";
    }
}
