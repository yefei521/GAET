use strict;
my$dir=$ENV{'PWD'};
my(%GFF, %CDS, @temp, @chr, %fa, $Par, $name);
open(IN1,"<$dir/gffcmp.combined.gtf_rep_more2");
open(OUT1,">$dir/gffcmp.combined.gtf_rep_more2_intron.gff");

#hap_tig00000120_pilon_pilon     AUGUSTUS        gene    1828    3612    1       +       .       ID=hap_tig00000120_pilon_pilon.g1
#hap_tig00000120_pilon_pilon     AUGUSTUS        transcript      1828    3612    1       +       .       ID=hap_tig00000120_pilon_pilon.g1.t1;Parent=hap_tig00000120_pilon_pilon.g1
#hap_tig00000120_pilon_pilon     AUGUSTUS        start_codon     1828    1830    .       +       0       Parent=hap_tig00000120_pilon_pilon.g1.t1
#hap_tig00000120_pilon_pilon     AUGUSTUS        CDS     1828    3612    1       +       0       ID=hap_tig00000120_pilon_pilon.g1.t1.cds;Parent=hap_tig00000120_pilon_pilon.g1.t1
#hap_tig00000120_pilon_pilon     AUGUSTUS        stop_codon      3610    3612    .       +       0       Parent=hap_tig00000120_pilon_pilon.g1.t1
foreach(<IN1>){   if(/^(chr_\d+)\tStringTie\texon\t(\d+)\t(\d+)\t.\t(.)\s.+?transcript_id "(.+?)";/){ my$start=$2-1; my$end=$3+1;  $CDS{"$1\t$4\t$5"}.="$start\t$end\t";   }   }

foreach my$name(keys%CDS){ @temp=split(/\t/,$CDS{$name}); pop@temp; shift@temp;
       next if @temp==0; 
       @chr=split(/\t/,$name);
       my$i=0; my$j=1;  
       while($i<@temp){ 
           print OUT1"$chr[0]\tAUGUSTUS\tgene\t$temp[$i]\t$temp[$i+1]\t.\t$chr[1]\t.\tID=$chr[2].$j\n$chr[0]\tAUGUSTUS\tCDS\t$temp[$i]\t$temp[$i+1]\t.\t$chr[1]\t.\tID=$chr[2].$j.cds;Parent=$chr[2].$j\n";  
           $i+=2;  $j++;}
}
close IN1; close OUT1;
my$cwd=`gffread gffcmp.combined.gtf_rep_more2_intron.gff -g ~/learning/TGD_Reference/TGD2020_latest/1-Genome_assembly_new.fasta -x gffcmp.combined.gtf_rep_more2_intron.fasta`;
print "$cwd";


open(IN2,"<$dir/gffcmp.combined.gtf_rep_more2_intron.fasta");
open(IN3,"<$dir/gffcmp.combined.gtf_rep_more2");
open(OUT2,">$dir/gffcmp.combined.gtf_rep_more2_GT-AG");
foreach(<IN2>){
   if(/^>((.+?)\.\d+)\n/){ $Par=$1; $name=$2; }else{chomp; $fa{$Par}.=$_;}
}
foreach(<IN3>){
   if(/transcript_id "(.+?)";.+? class_code "j";.+? num_samples "(.)";/){  if($2>1){ $name=$1; $GFF{$name}=$_;}
   }elsif(/transcript_id "($name)";/){$GFF{$name}.=$_;  }
}
foreach(keys%fa){ $_=~/^(.+?)\.\d+/; $name=$1;
   if($fa{$_}=~/^GT.+AG$/){  
   }elsif($fa{$_}=~/^CT.+AC$/){ 
   }elsif($fa{$_}=~/^GC.+AG$/){ 
   }elsif($fa{$_}=~/^CT.+GC$/){ #print OUT2"$GFF{$name}";
   }else{ delete$GFF{$name};}
}
foreach(keys%GFF){print OUT2"$GFF{$_}"; }
my$cmd=`gffread gffcmp.combined.gtf_rep_more2_GT-AG  -g ~/learning/TGD_Reference/TGD2020_latest/1-Genome_assembly_new.fasta -w gffcmp.combined.gtf_rep_more2_GT-AG_exon.fasta`;
print "$cmd";
