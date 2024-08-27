use strict;
my(%Error, $supercontig, $gene_id, %mRNA, %gff, );
my$dir=$ENV{'PWD'};
open(IN1,"<$dir/01Error_gene.txt");
open(IN2,"<$dir/02ATG2TGA_gene_result.txt");
open(IN3,"<$dir/Thermophila_tetrahymena_2021_UTR_manualcheck_220224.gff3");
open(IN4,"<$dir/05intron_len2long_name.txt");
open(OUT,">$dir/03T-t-manualcheck_220224_rmErrorGene.gff3");
#-       g13661  8
#+       g24920  8
#+       g9331   8
while(<IN1>){ if(/^.\s+(.+?)\s+/){ $Error{$1}++;  } }
while(<IN4>){ if(/^(.+?)\.t/){$Error{$1}++; } }

#g21752.t1       ATG-TGA yes
#g4502.t1        TTT-TAC no
#g16296.t1       AGA-GTT no
while(<IN2>){ if(/^(.+)\..+?(no)/){ $Error{$1}++;  }  }

while(<IN3>){
       if(/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+.+?ID=(chr_\d+);/){ $supercontig++;   print OUT"$_";
       }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5; $mRNA{$gene_id}=$6; $gff{$gene_id}.=$_; if(exists$Error{$gene_id}){}else{print OUT"$_"; }
       }elsif(/chr_\d+.+?mRNA.+?Parent=$gene_id/){  $gff{$gene_id}.=$_; if(exists$Error{$gene_id}){}else{ print OUT"$_"; }
       }elsif(/chr_\d+.+?UTR.+?Parent=$gene_id/){   $gff{$gene_id}.=$_; if(exists$Error{$gene_id}){}else{ print OUT"$_"; }
       }elsif(/chr_\d+.+?CDS.+?Parent=$gene_id/){   $gff{$gene_id}.=$_; if(exists$Error{$gene_id}){}else{ print OUT"$_"; }
       }elsif(/chr_\d+.+?MANUAL_ENTRY.+?ID=(.+?);Name/){ print OUT"$_";
       }elsif(/chr_\d+.+?rDNA.+?ID=(.+?);Name/){  print OUT"$_"; # push@gff,$1; $gff{$1}=$_;
       }
    }
