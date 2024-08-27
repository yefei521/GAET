use strict;
my$dir=$ENV{'PWD'};
my();
open(IN1,"<","$dir/gffcmp.combined.gtf_rep_more2_GT-AG_UTR.gff3_Filtered");
open(IN2,"<","$dir/Upset-plot/7Merged_RMoverlap_gene.gff3");
my@IN1=<IN1>;
my@IN2=<IN2>;

foreach(@IN1){
       if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5;  $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?CDS.+?Parent=$gene_id/){   $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?five_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){   $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?three_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){   $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?mRNA.+?Parent=$gene_id/){  $gff{$gene_id}.=$_;
       }else{print "$_";}
}


foreach(@IN1){
   
}


