use strict;
my$dir=$ENV{'PWD'};
my(%intron, $chr, $gene_id, %gene, @temp, @CDS, $gene_id, $length, %aver_intron_length, $n);
open(IN1,"<$dir/04T-t-manualcheck_220224_rmErrorGene_rmOverlapGene.gff3");
open(IN2,"<$dir/05Nano_new_gene.gtf");
#open(IN3,"<$dir/");
open(OUT2,">$dir/05intron_len2long_name.txt");
foreach(<IN1>){
       if(/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+(\d+)\s.+?ID=(chr_\d+);/){   $chr=$1;
       }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5;  
       }elsif(/chr_\d+.+?mRNA.+?Parent=$gene_id/){  
       }elsif(/chr_\d+.+?UTR.+?Parent=$gene_id/){   
       }elsif(/chr_\d+.+?CDS.+?Parent=$gene_id/){  $gene{$gene_id}.=$_;  
       }elsif(/chr_\d+.+?MANUAL_ENTRY.+?ID=(.+?);Name/){  
       }elsif(/chr_\d+.+?rDNA.+?ID=(.+?);Name/){ 
       }
    }
my$n=keys%gene; print "$n\n"; my$m=0;
foreach(keys%gene){
      @temp=split(/\n/,$gene{$_}); undef @CDS;
      foreach my$line(@temp){  if($line=~/CDS\s+(\d+)\s+(\d+)\s.+?Parent=(.+)/){  push@CDS,$1,$2; $gene_id=$3; }  }
      shift@CDS; pop@CDS; $n=0; $length=0; 
      while($n<$#CDS){ $length+=$CDS[$n+1]-$CDS[$n]; $n+=2;  $intron{$m}=$CDS[$n+1]-$CDS[$n]; $m++; if( ($CDS[$n+1]-$CDS[$n])>2000 ){print OUT2"$gene_id\t$CDS[$n]\t$CDS[$n+1]\n"; } }
      if(@CDS>0){ $aver_intron_length{$gene_id}=(2*$length)/@CDS; }      
}
open(OUT1,">$dir/05intron.txt");
foreach(keys%aver_intron_length){   print OUT1"$_\t$aver_intron_length{$_}\n"; }
