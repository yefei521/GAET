use strict;
my(%AS, %AS_other, $chr, $chr_len, @IR, %gff);
my$dir=$ENV{'PWD'};
open(IN1,"<","$dir/gffcmp.combined.gtf_rep_more2_GT-AG_UTR.gff3_Filtered_intron.bed");
open(IN2,"<","$dir/7Merged_RMoverlap_gene.gff3_intron.bed");
open(IN3,"<","$dir/chr-list.txt");
open(IN4,"<","$dir/Upset-plot/7Merged_RMoverlap_gene.gff3");
my@IN4=<IN4>; close IN4;
my@IN1=<IN1>; close IN1;
my@IN2=<IN2>; close IN2;
my@IN3=<IN3>; close IN3;
open(OUT1,">","$dir/08-03Alter-5-splice-site.bed");
foreach(@IN3){
   if(/SN:(chr_\d+)\s+LN:(\d+)\s/){ $chr=$1; $chr_len=$2; undef%AS;undef%AS_other;
      foreach(@IN2){if(/^$chr\s+(\d+)\s+(\d+)\s/){ if(!exists$AS{$1} or $AS{$1}==$2){ $AS{$1}=$2;}else{$AS_other{$1}=$2; }  }  }
      foreach my$line(@IN1){if($line=~/^$chr\s+(\d+)\s+(\d+)\s+\+/){ my$start=$1; my$end=$2; chomp$line; my$MARK3=0; my$MARK5=0; my$a=0; my$b=0;
            if(!exists$AS{$1} and !exists$AS_other{$1}){ $MARK5++; }
            foreach my$keys(keys%AS){if($end==$AS{$keys}){$a=1;} } if($a>0){$MARK3++; }
            foreach my$keys(keys%AS_other){if($end==$AS_other{$keys}){$b=1;} } if($b>0){$MARK3++; }
            if($MARK5==0 and $MARK3==0){  push@IR,$line;
                foreach(@IN4){if(/$chr.+?mRNA\s+(\d+)\s+(\d+).+?ID=(.+?);/){ if($1<$start && $2>$end){print OUT1"$line\n"; }  }   }
            }
         }
      }
      foreach my$line(@IN1){if($line=~/^$chr\s+(\d+)\s+(\d+)\s+\-/){ my$start=$1; my$end=$2; chomp$line; my$MARK3=0; my$MARK5=0; my$a=0; my$b=0;
            if(!exists$AS{$1} and !exists$AS_other{$1}){ $MARK3++; }
            foreach my$keys(keys%AS){if($end==$AS{$keys}){$a=1;} } if($a>0){$MARK5++; }
            foreach my$keys(keys%AS_other){if($end==$AS_other{$keys}){$b=1;} } if($b>0){$MARK5++; }
            if($MARK5==0 and $MARK3==0){  push@IR,$line;
                foreach(@IN4){if(/$chr.+?mRNA\s+(\d+)\s+(\d+).+?ID=(.+?);/){ if($1<$start && $2>$end){print OUT1"$line\n"; }  }   }
            }
         }
      }
   }
}
