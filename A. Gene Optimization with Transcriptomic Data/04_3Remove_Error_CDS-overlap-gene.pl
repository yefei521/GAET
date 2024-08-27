use strict;
my$dir=$ENV{'PWD'};
my(%Error_gene, $gene_id, );
open(IN1,"<$dir/03T-t-manualcheck_220224_rmErrorGene.gff3") or die;
open(IN2,"<$dir/04Overlap_gene_CDS.txt") or die; #Repeat gene which overlap CDS
open(OUT1,">$dir/04T-t-manualcheck_220224_rmErrorGene_rmOverlapGene.gff3") or die;

#Error
#chr_001 AUGUSTUS        gene    873443  875279  .       +       .       "ID=g18644;Name=TTHERM_00155430;Note=""NTP9 equilibrative nucleoside transporter family protein"""
#chr_001 AUGUSTUS        gene    873481  873714  .       -       .       ID=WT_C5_2_rmdup_STRG.251;Name=WT_C5_2_rmdup_STRG.251.1;Note=
foreach(<IN2>){
   if(/ID=(.+?);Name=(.+?);Note/){ $Error_gene{$1}++;  }
}
my@GFF=<IN1>; close IN1;
foreach(@GFF){
       if(/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+.+?ID=(chr_\d+);/){  print OUT1"$_";
       }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5; if(exists$Error_gene{$gene_id}){}else{ print OUT1"$_"; }
       }elsif(/chr_\d+.+?mRNA.+?Parent=$gene_id/){  if(exists$Error_gene{$gene_id}){}else{ print OUT1"$_"; }
       }elsif(/chr_\d+.+?UTR.+?Parent=$gene_id/){   if(exists$Error_gene{$gene_id}){}else{ print OUT1"$_"; }
       }elsif(/chr_\d+.+?CDS.+?Parent=$gene_id/){   if(exists$Error_gene{$gene_id}){}else{ print OUT1"$_"; }
       }elsif(/chr_\d+.+?MANUAL_ENTRY.+?ID=(.+?);Name/){ print OUT1"$_";
       }elsif(/chr_\d+.+?rDNA.+?ID=(.+?);Name/){ print OUT1"$_";
       }
    }


