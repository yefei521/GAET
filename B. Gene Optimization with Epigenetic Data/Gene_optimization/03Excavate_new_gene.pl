#This script used to excavate potential gene acorrding result 02New_gene_accoding_6mA_distribution.txt 
use strict;
my$dir=$ENV{'PWD'};
my(@gtf, @m6A, @m6A_corrd, $chr, $mark, %transcript, %transcript_FPKM, $transcript_id, $line, %transcript_id, %transcript_fasta);
#-----------------------------
#chr_070 gene12948 group1: 11757|-3778    11760|-3781  
print"STEP1: STrating import data from file 02New_gene_accoding_6mA_distribution.txt\nPlease wait .............\n";
#open(IN1,"<$dir/02New_gene_accoding_6mA_distribution.txt");
open(OUT1,">$dir/03Excavate_new_transcript.gtf");
open(OUT2,">$dir/03gene_with_6mA_need_check.gtf");
@gtf=glob"$dir/*rename.gtf";
my@txt=glob"$dir/02*txt";
foreach my$txt(@txt){ chomp $txt;
  open(TXT,"<$txt") or die;
  foreach my$m(<TXT>){ undef @m6A; undef@m6A_corrd;
   if($m=~/^(chr_\d+)\s+.+?group\d+/){ 
      @m6A=split(/\s+/,$m); $chr=$1;
      foreach my$i(@m6A){ push@m6A_corrd,$i if($i=~/^(\d+)\|.+/); }
      open(GENE,"$dir/chr_gene/$chr.txt") or die;
      foreach(<GENE>){
          if(/$chr.+?gene\t(\d+)\t(\d+)\t/){if($1<$m6A_corrd[0] && $2>$m6A_corrd[-1]){$mark=1; $line=$_; }else{$mark=0;}  }
          last if $mark==1;
      }
      print OUT2"$m$line" if $mark==1;
      if($mark==0){
         foreach(@gtf){ chomp; open(GTF,"<$_")or die;
            foreach my$j(<GTF>){ 
               if($j=~/$chr\tStringTie\ttranscript\t(\d+)\t(\d+)\t.+?transcript_id "(.+?)";.+?FPKM "(.+?)"/){
                  if($1<$m6A_corrd[0] && $2>$m6A_corrd[-1]){ $transcript{$3}=$j; $transcript_FPKM{$3}=$4; $transcript_id=$3; print OUT1"\n$m$j"; $transcript_id{$3}=$3; }
              }elsif($j=~/$chr.+?exon.+?"$transcript_id";/){ $transcript{$3}.=$j;  print OUT1"$j"; }             
            }
         }
      } 
   }
 }
}
##
#Strategy: extract all of trancripts(including multi-transcripts for one gene) to predict ATG-TGA, because maybe not the higest FPKM transcript will be the RIGHT gene model.
#---------------------------------------
print"STEP2: Starting import fasta from file \nPlease wait .............\n";
open(OUT3,">$dir/03Excavate_new_transcript.fasta");
my@fasta_file=glob"$dir/*rename.fasta";
foreach my$i(@fasta_file){chomp $i;
    open(FASTA,"<$i");
    foreach my$j(<FASTA>){
        if($j=~/^>(.+?)\n/){ $transcript_id=$1;
       }else{ $transcript_fasta{$transcript_id}.=$j }
    }
}
my$n=keys%transcript_id;
foreach(keys%transcript_id){
    print OUT3">$_\n$transcript_fasta{$_}";
}
print"$n transcripts fasta exported to be predicted.\n FINiSHED thank you!\n";
#------------------------------
