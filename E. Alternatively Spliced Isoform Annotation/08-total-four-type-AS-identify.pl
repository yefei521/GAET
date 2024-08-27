use strict;
my(%AS, %AS_other, $chr, $chr_len, %TAS, %TAS_other, @lines);
my$dir=$ENV{'PWD'};
open(IN1,"<","$dir/gffcmp.combined.gtf_rep_more2_GT-AG_UTR.gff3_Filtered_intron.bed");
open(IN2,"<","$dir/7Merged_RMoverlap_gene.gff3_intron.bed");
open(IN3,"<","$dir/chr-list.txt");
open(IN4,"<","$dir/7Merged_RMoverlap_gene.gff3");
open(IN5,"<","$dir/gffcmp.combined.gtf_rep_more2_GT-AG_UTR.gff3");

open(OUT1,">","$dir/08-TA-cassette-exon-skip.bed");
open(OUT2,">","$dir/08-TA-Alter-5-splice-site.bed"); #5'AS
open(OUT3,">","$dir/08-TA-Alter-3-splice-site.bed"); #3'AS
open(OUT4,">","$dir/08-TA-Mutually-exclulsive-exon.bed");
open(OUT5,">","$dir/08-TA-Intron-retetion.bed");

my@IN1=<IN1>;
my@IN2=<IN2>;
my@IN3=<IN3>;
my@GFF=<IN4>; close IN4;
my@GFF_AS=<IN5>; close IN5;
foreach(@IN3){ 
   if(/SN:(chr_\d+)\s+LN:(\d+)\s/){ $chr=$1; $chr_len=$2; undef%AS;undef%AS_other; # print "$chr\n";
      foreach(@IN2){if(/^$chr\s+(\d+)\s+(\d+)\s/){ 
            if(exists$AS{$1} or $AS{$1}==$2){ $AS{$1}=$2;}else{$AS{$1}=$2; }  
            if(exists$TAS{$2} or $TAS{$2}==$1){ $TAS{$2}=$1;}else{$TAS{$2}=$1; } }  }
      foreach my$line(@IN1){if($line=~/^$chr\s+(\d+)\s+(\d+)\s/){ my$start=$1; my$end=$2; chomp$line; my$mark=0; my$mark2=0;
            if(exists$AS{$start} and $AS{$start}!=$end ){  #print "$line\n";
                 foreach(keys%AS){ 
                      if($end==$AS{$_} && $TAS{$end}>$start){print OUT1"$line\t$start,$AS{$start}\t$_,$AS{$_}\n"; $mark++;   
                      }elsif($end==$AS{$_} && $TAS{$end}>=$start){ print OUT5"$line\t$start,$AS{$start}\t$_,$AS{$_}\n"; $mark++; }  }
		 if($mark==0 && exists$AS{$start} and $AS{$start}<$end){ print OUT3"$line\t$start,$AS{$start}\t$start,$AS{$start}\n";   }   
                 if($mark==0 && exists$AS{$start} and $AS{$start}>$end){ print OUT3"$line\t$start,$AS{$start}\t$start,$AS{$start}\n";   }
            }elsif(!exists$AS{$start} && exists$TAS{$end}){ 
                 if(exists$TAS{$end} and $TAS{$end}<$start){ print OUT2"$line\t$start,$AS{$start}\t$start,$AS{$start}\n";   }
                 if(exists$TAS{$end} and $TAS{$end}>$start){ print OUT2"$line\t$start,$AS{$start}\t$start,$AS{$start}\n";   }
            }elsif(exists$AS{$start} and $AS{$start}==$end){  
            }else{push@lines,$line;  }   
      }  
      }
   }
}

my(%gff, $gene_id, @lines2);
foreach(@GFF){ 
       if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5;  $gff{$gene_id}.=$_; 
       }elsif(/chr_\d+.+?CDS.+?Parent=$gene_id/){   $gff{$gene_id}.=$_;
       }elsif(/chr_\d+.+?five_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){   $gff{$gene_id}.=$_; 
       }elsif(/chr_\d+.+?three_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){   $gff{$gene_id}.=$_; 
       }elsif(/chr_\d+.+?mRNA.+?Parent=$gene_id/){  $gff{$gene_id}.=$_;
       } #else{print "$_";}
}
my$n=keys%gff; print "$n\n";
my(%hash, @num,$num, $id, $UTR5_num, $UTR3_num, $UTR5_len, $UTR3_len, $mark, $chr);
foreach my$name(keys%gff){ my@temp=split(/\n/,$gff{$name}); undef@num;
   foreach my$line(@temp){
       if($line=~/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ $gene_id=$5; $chr=$1;
       }elsif($line=~/chr_\d+.+?CDS\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){    push@{$hash{$chr}{$gene_id}},$1,$2;
       }elsif($line=~/chr_\d+.+?five_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){    push@{$hash{$chr}{$gene_id}},$1,$2;
       }elsif($line=~/chr_\d+.+?three_prime_UTR\s+(\d+)\s+(\d+)\s.+?Parent=$gene_id/){   push@{$hash{$chr}{$gene_id}},$1,$2;
       }elsif($line=~/chr_\d+.+?mRNA.+?Parent=$gene_id/){  }        
  }
}
my$mark=0;
foreach my$line(@lines){
    if($line=~/^(chr_\d+)\s+(\d+)\s+(\d+)\s/){$chr=$1; my$start=$2; my$end=$3; $mark=0; 
      foreach my$i(keys%{$hash{$chr}}){ 
          my$j=0; my$n=@{$hash{$chr}{$i}}; my$k=0;
          while($j<$n){ 
               if($hash{$chr}{$i}[$j+2]-$hash{$chr}{$i}[$j+1]==1){$k=$j+3;  }else{$k=$j+1; }
               if($hash{$chr}{$i}[$j]<=$start && $hash{$chr}{$i}[$k]>=$end ){ print OUT4"$line\t$hash{$chr}{$i}[$j],$hash{$chr}{$i}[$k]\t$start,$end\n"; $mark++; last; } 
               $j=$k+1;  }
          last if $mark>0;
      }
      print OUT5"$line\t$start,$end\n";
   }
}

