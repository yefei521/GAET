my$dir=$ENV{'PWD'};
open(IN1,"<","$dir/WT_mix_seqkit2DNA_MaxIntron2k.sam_polya.txt") or die;
open(IN2,"<","$dir/polya_seq_2ManUTR+lowQ_exon_otf7_filter") or die;

foreach(<IN1>){
   if(/^(chr_\d+)\s+(\d+)\s+(\d+)\s+(.+?-.+?)\s+(\d+)\s+(.)\s/){ #chr_001 235664  236255  3dc132ab-8754-42da-9d29-c2b542cffdaa    14      -       TCTAAAAAAAAAAAAGA       polya:AAAAAAAAAAAAGA
        $polya_len{$4}=$5;
    }
}

open(OUT1,">","$dir/ManUTR+lowQ_gene_polyA_length.txt");
print OUT1"name\tpolyA_num\tmedian-length\tlength1\tlength2\t...\n";
foreach(<IN2>){
    @a=split(/\t/); #YF00007972.t1   38149476-2524-4f9b-a8b5-bffc6c75c61b    b741d5b0-f358-48a5-a3d8-eb0c8c034cde
    $name=shift@a; $num_polya=0; undef@len;
    foreach $id(@a){ if($id=~/.+?-.+?/){ push@len,$polya_len{$id}; $num_polya++; } }
    @sort_len=sort{$b<=>$a}@len;
    $median_id=int(@sort_len/2)-1;
    if($median_id<0){$median_id=0;}
    print OUT1"$name\t$num_polya\t$sort_len[$median_id]\t";
    foreach $len(@sort_len){print OUT1"$len\t"; }
    print OUT1"\n";
}
