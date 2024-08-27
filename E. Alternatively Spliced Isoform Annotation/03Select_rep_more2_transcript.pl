use strict;
my$dir=$ENV{'PWD'};
my(%gff, $name);
open(IN1,"<$dir/gffcmp.combined.gtf");
open(OUT1,">$dir/gffcmp.combined.gtf_rep_more2");

foreach(<IN1>){
   if(/transcript_id "(.+?)";.+? class_code "j";.+? num_samples "(.)";/){  if($2>1){ $name=$1; $gff{$name}=$_;}
   }elsif(/transcript_id "($name)";/){$gff{$name}.=$_;  }
}

foreach(keys%gff){   print OUT1"$gff{$_}";  }


#chr_001 StringTie       exon    109813  109852  .       +       .       transcript_id "total_00000021"; gene_id "XLOC_000002"; exon_number "1";
#foreach(keys%gff){
#   @temp=split(/\n/,$gff{$_}); undef@exon;
#   foreach my$a(@temp){
#      if(/transcript_id "(.+?)";.+? class_code "j";.+? num_samples "(.)";/){ 
#      }elsif(/^(chr_\d+)\tStringTie\texon\t(\d+)\t(\d+)\t.\t(.)\t.+?/){ push@exon,$2,$3;  $chr=$1; $ori=$4; }
#   }shift@exon; pop@exon;
#   my$b=0;while($b<@exon){my$start=$exon[$i]+1; my$end=$exon[$i]+2; my$cmd=`samtools faidx ~/learning/TGD_Reference/TGD2020_latest/1-Genome_assembly_new.fasta $chr:$start-$end`; }
#}
