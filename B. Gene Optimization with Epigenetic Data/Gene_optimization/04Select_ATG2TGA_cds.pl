#perl 01Select_ATG2TGA_cds.pl total_transcripts_exon.fasta_orf-1
use strict;
my(@cds1,@cds2,@list,$n,$name,%list);
open (IN1, "<", "03Excavate_new_transcript.fasta_orffinder_otf1") or die "Can't open IN1 : $!";
open(IN2,"<","03Excavate_new_transcript.fasta_orffinder_minus_otf1") or die "Can't open IN1 : $!" ;
open (OUT1, ">", "04temp_plus.fa");
open (OUT2, ">", "04temp_plus_name");
 
   {@cds1=<IN1>;close IN1;} #open plus, select ATG2TGA fasta
   $n=0;
   #>lcl|WT_C10_2_rmdup_STRG.161.1:244-582 ORF1_WT_C10_2_rmdup_STRG.161.1:243:581
   my(%fasta,%cds_ATG2TGA);
   foreach(@cds1){
       if(/^>.+?(ORF\d+)_(.+?_STRG.\d+\.\d+):(\d+):(\d+)\n/){
          my$len=$4-$3+1;
          $name=">$2\t$1\t$3\t$4\t$len\tplus"; $n=1; 
       }elsif($n==1){
          chomp;
          $fasta{$name}.=$_;  }
       }
   foreach(keys %fasta){
      if($fasta{$_}=~/^ATG.+?TGA$/){
          print OUT1 "$_\n$fasta{$_}\n";
          my@a=split />/;
          print OUT2 "$a[1]\n";}
   }
   undef %fasta; undef %cds_ATG2TGA; 
   close OUT1; close OUT2;
 
#open minus , select ATG2TGA fasta
  
   open (OUT3, ">", "04temp_minus.fa");
   open (OUT4, ">", "04temp_minus_name");
   {@cds2=<IN2>;close IN2;}
   $n=0;
   my(%fasta,%cds_ATG2TGA);
   foreach(@cds2){
       if(/^>.+?(ORF\d+)_(.+?_STRG.\d+\.\d+):(\d+):(\d+)\n/){
          my$len=$3-$4+1;
          $name=">$2\t$1\t$3\t$4\t$len\tminus"; $n=1;
       }elsif($n==1){
          chomp;
          $fasta{$name}.=$_;  }
       }
   foreach(keys %fasta){
      if($fasta{$_}=~/^ATG.+?TGA$/){
          print OUT3 "$_\n$fasta{$_}\n";
          my@a=split />/;
          print OUT4 "$a[1]\n";}
   }
   close OUT3; close OUT4;
   undef %fasta; undef %cds_ATG2TGA;

my($mark,%fasta,%max, %longest,$name,$name_to,$length);
my$dir=$ENV{'PWD'};
open(PLUS,"<$dir/04temp_plus.fa");
open(MINUS,"<$dir/04temp_minus.fa");
open(LONG,">$dir/04orffinder_ATG2TGA_longest.fa");
#>WT_C10_2_rmdup_STRG.4448.1     ORF2    894     1190    297     plus
foreach my$i(<PLUS>){chomp$i;
    if($i=~/>(.+?)\t(ORF\d+)\t(\d+)\t(\d+)\t(\d+)\tplus/){
        if($5>$max{$1}){ $longest{$1}=$i; $max{$1}=$5; $mark=1; $name=$i;}else{$mark=0;}
    }elsif($mark==1){$fasta{$name}=$i; $mark=0;}
}

foreach my$i(<MINUS>){ chomp$i;
    if($i=~/>(.+?)\t(ORF\d+)\t(\d+)\t(\d+)\t(\d+)\tminus/){
       $name=$1; $length=$5; $name_to=$i;
       if($length>$max{$name}){$longest{$name}=$i; $max{$name}=$length; $mark=1;}else{$mark=0;}
    }elsif($mark==1){$fasta{$name_to}=$i; $mark=0;}
}
foreach(keys%longest){chomp;
    if($longest{$_}=~/>(.+?)\t(ORF\d+)\t(\d+)\t(\d+)\t(\d+)/){
        #my$name=">$longest{$_}";
        print LONG"$longest{$_}\n$fasta{$longest{$_}}\n";
    }
}
