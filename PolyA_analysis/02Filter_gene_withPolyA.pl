$dir=$ENV{'PWD'};
open(IN1,"<","$dir/polya_seq_2ManUTR+lowQ_exon_otf7");
while(<IN1>){
    next if $_=~/^#/;
    @a=split(/\t/,$_);
    next if exists$seq{$a[0]};
    $seq{$a[0]}++;
    push@{$hash{$a[1]}},$a[0];
}

open(IN2,"<","$dir/WT_mix_seqkit2DNA_MaxIntron2k.sam_polya.txt");
#chr_001 235664  236255  3dc132ab-8754-42da-9d29-c2b542cffdaa    14      -       TCTAAAAAAAAAAAAGA       polya:AAAAAAAAAAAAGA
while(<IN2>){
    if(/^(chr_\d+)\s+(\d+)\s+(\d+)\s+(.+?)\s+\d+\s+(.)\s/){
       if($5 eq "+"){   $pos_TES{$4}=$3; 
       }elsif($5 eq "-"){ $pos_TES{$4}=$2;  }
    }
}


#$cmd=`rm -r $dir/gene_TES`; print $cmd;
$cmd=`mkdir -p $dir/gene_TES`; print $cmd;
open(OUT1,">","$dir/polya_seq_2ManUTR+lowQ_exon_otf7_filter");
open(OUT2,">","$dir/polya_seq_2ManUTR+lowQ_exon_otf7_pos_TES");
foreach $name(keys%hash){
   print OUT1 "$name\t";
   print OUT2 "$name\t";
   $out=$name;
   open($out,">","$dir/gene_TES/$name");
   print $out "name\tpos\tnum\n";
   foreach $seqid(@{$hash{$name}}){
       print OUT1 "$seqid\t"; 
       print OUT2 "$pos_TES{$seqid}\t";
       print $out "$seqid\t$pos_TES{$seqid}\t1\n";
    }
   print OUT1 "\n"; print OUT2 "\n";
   close $out;
}

