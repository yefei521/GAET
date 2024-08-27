use strict;
my(%AS, %AS_other, $chr, $chr_len, );
my$dir=$ENV{'PWD'};
open(IN1,"<","$dir/gffcmp.combined.gtf_rep_more2_GT-AG_UTR.gff3_Filtered_intron.bed");
open(IN2,"<","$dir/7Merged_RMoverlap_gene.gff3_intron.bed");
open(IN3,"<","$dir/chr-list.txt");
open(OUT1,">","$dir/08-01cassette-exon-skip.bed");
my@IN1=<IN1>;
my@IN2=<IN2>;
my@IN3=<IN3>;
foreach(<IN3>){ 
   if(/SN:(chr_\d+)\s+LN:(\d+)\s/){ $chr=$1; $chr_len=$2; undef%AS;undef%AS_other;  
      foreach(@IN2){if(/^$chr\s+(\d+)\s+(\d+)\s/){ if(!exists$AS{$1} or $AS{$1}==$2){ $AS{$1}=$2;}else{$AS_other{$1}=$2; }  }  }
      foreach my$line(@IN1){if($line=~/^$chr\s+(\d+)\s+(\d+)\s/){ my$start=$1; my$end=$2; chomp$line; 
            if(exists$AS{$1} and $AS{$1}!=$2 ){ foreach(keys%AS){if($end==$AS{$_}){print OUT1"$line\t$start,$AS{$start}\t$_,$AS{$_}\n";  }  }   }   }  }
   }
}

