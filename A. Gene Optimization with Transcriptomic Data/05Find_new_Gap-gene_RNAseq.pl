use strict;
my($n,@chr, $chr, @gene_corrd, %gene, $contig_end, $gene_id, $gene_start, $gene_end, $former, $next, %new_nano_gene, $min);
my$dir=$ENV{'PWD'};
open(IN1,"<$dir/04T-t-manualcheck_220224_rmErrorGene_rmOverlapGene.gff3");
# open(IN1,"<$dir/04chr_066_temp.gff3") or die;
open(IN2,"<$dir/SB210_DRS_DNA_Q7_maped_sorted.gtf"); #nanopore RNA for define the region 
#open(IN2,"<$dir/temp_nano,gtf") or die;
open(IN3,"<$dir/veg_assembly_accepted_hits_unique_sorted_rmdup.gtf") or die;
open(IN4,"<$dir/starvation_assembly_accepted_hits_unique_sorted_rmdup.gtf") or die;
open(IN5,"<$dir/conjugation_assembly_accepted_hits_unique_sorted_rmdup.gtf") or die;

my@GFF=<IN1>; close IN1; foreach(@GFF){ if(/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+(\d+)\s.+?ID=(chr_\d+);/){ push@chr,$1;  }  }
my@NAnopore_GTF=<IN2>;
foreach $chr(@chr){
   print "Processing:$chr\n"; undef@gene_corrd; undef%gene;
   foreach(@GFF){
       if(/($chr)\s+.+?supercontig\s+(\d+)\s+(\d+)\s.+?ID=(chr_\d+);/){  push@gene_corrd,$2; $contig_end=$3; 
       }elsif(/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5; push@gene_corrd,$2,$3; $gene{$_}++;
       }elsif(/$chr.+?mRNA.+?Parent=$gene_id/){  
       }elsif(/$chr.+?UTR.+?Parent=$gene_id/){   #if(exists$Error_gene{$gene_id}){}else{ print OUT1"$_"; }
       }elsif(/$chr.+?CDS.+?Parent=$gene_id/){   #if(exists$Error_gene{$gene_id}){}else{ print OUT1"$_"; }
       }elsif(/$chr.+?MANUAL_ENTRY.+?ID=(.+?);Name/){  #print OUT1"$_";
       }elsif(/$chr.+?rDNA.+?ID=(.+?);Name/){ #print OUT1"$_";
       }
    }
    push@gene_corrd,$contig_end; 
    foreach(keys%gene){
        if(/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){
           $gene_start=$2; $gene_end=$3; 
           $former=&MAX_former;  $next=&MIN_next; $gene_id="";
           foreach my$line(@NAnopore_GTF){ 
               if($line=~/$chr\tStringTie\ttranscript\t(\d+)\t(\d+)\t.+?transcript_id\s+"(.+?)";.+?cov\s+"(.+?)";/){  
                     if($1>$former && $2<$gene_start && $4>5){ $gene_id=$3; $new_nano_gene{$gene_id}=$line;  }
                     if($1>$gene_end && $2<$next && $4>5){ $gene_id=$3; $new_nano_gene{$gene_id}=$line;  }
               }elsif($line=~/$chr\tStringTie\texon\t(\d+)\t(\d+)\t.+?transcript_id\s+"($gene_id)";/){
                     $new_nano_gene{$gene_id}.=$line; 
               }
           }$gene_id="";
        }
    }
}

sub MAX_former{ 
    $min=$gene_corrd[-1]-$gene_corrd[0]; my$mark=$gene_start;
    foreach(sort{$a<=>$b}@gene_corrd){ if($_<$gene_start && ($gene_start-$_)<$min ){ $min=$gene_start-$_; $mark=$_; }  }
     $mark;
}
sub MIN_next{ 
    $min=$gene_corrd[-1]-$gene_corrd[0]; my$mark=$gene_end;
    foreach(sort{$a<=>$b}@gene_corrd){ if($_>$gene_end && ($_-$gene_end)<$min ){ $min=$_-$gene_end; $mark=$_; }  }
     $mark;
}
open(OUT1,">$dir/05Nano_new_gene.gtf");
foreach(keys%new_nano_gene){  print OUT1"$new_nano_gene{$_}"; }

