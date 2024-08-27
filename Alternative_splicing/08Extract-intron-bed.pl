use strict;
my(%CDS, @corrd, $name, @CDS, $chr, $ori);
my$dir=$ENV{'PWD'};
open(IN1,"<","$dir/Upset-plot/7Merged_RMoverlap_gene.gff3");
open(OUT1,">","$dir/7Merged_RMoverlap_gene.gff3_intron.bed");
#open(IN1,"<","$dir/gffcmp.combined.gtf_rep_more2_GT-AG_UTR.gff3_Filtered");
#open(OUT1,">","$dir/gffcmp.combined.gtf_rep_more2_GT-AG_UTR.gff3_Filtered_intron.bed");
#chr_001 AUGUSTUS        mRNA    1929    7152    .       +       .       ID=g18411.t1;Parent=g18411
#chr_001 AUGUSTUS        CDS     1929    2780    .       +       0       ID=g18411.t1.cds;Parent=g18411.t1
#chr_001 AUGUSTUS        CDS     3123    3573    .       +       0       ID=g18411.t1.cds;Parent=g18411.t1
#chr_001 AUGUSTUS        CDS     3719    4573    .       +       2       ID=g18411.t1.cds;Parent=g18411.t1
#chr_001 AUGUSTUS        CDS     4636    5283    .       +       2       ID=g18411.t1.cds;Parent=g18411.t1
#chr_001 AUGUSTUS        CDS     5353    5563    .       +       2       ID=g18411.t1.cds;Parent=g18411.t1
#chr_001 AUGUSTUS        CDS     5614    6218    .       +       1       ID=g18411.t1.cds;Parent=g18411.t1
#chr_001 AUGUSTUS        CDS     6299    7152    .       +       2       ID=g18411.t1.cds;Parent=g18411.t1
foreach(<IN1>){
   if(/mRNA.+?ID=(.+?);Parent/){$name=$1; 
  }elsif(/^(chr_\d+)\s+.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=$name/){ $CDS{$name}.=$_; 
  }
}

foreach(keys%CDS){ @CDS=split(/\n/,$CDS{$_}); undef@corrd;
  foreach(@CDS){ 
     if(/^(chr_\d+)\s+.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?Parent=(.+)/){ $chr=$1; $ori=$4; $name=$5; push@corrd,$2,$3; } 
  }
  pop@corrd; shift@corrd;
  my$i=0;while($i<@corrd){ my$start=$corrd[$i]+1; my$end=$corrd[$i+1]-1; print OUT1"$chr\t$start\t$end\t$ori\t$name\n"; $i+=2; }
  
}

