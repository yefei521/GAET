use strict;
my$dir=$ENV{'PWD'}; 
my(@GFF, $chr, $gene_id, %gene, @Chr_gff, %temp, $gene_start, $gene_end, $n, %gene2line, $name, $mark, %GFF, );
my(@CDS, %CDS_ratio, %repeat_line, %repeat_gene, @repeat_GENE, @temp, $gene_name);
open(IN1,"<$dir/8Merge_raw_GapNew_gene.gff3") or die;
#open(IN1,"<$dir/temp.gff");
@GFF=<IN1>; close IN1;
#chr_001 AUGUSTUS        supercontig     1       1459056 .       .       .       ID=chr_001;Name=chr_001
my$cmd=`rm chr/chr*`;
print "$cmd";
foreach(@GFF){
       if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ 
           $gene_id=$5; $chr=$1;
           open(TXT,">>$dir/chr/$chr.gff3") or die; 
           print TXT"$_"; close TXT; 
           $gene{$_}=$_; $gene2line{$gene_id}=$_;
       }
    }
print"0\n";
foreach my$line(@GFF){
   if($line=~/gene.+?ID=((.+?)_ORF\d+);Name/){ $name=$1; $mark=$1; $GFF{$name}.=$line;
   }elsif($line=~/Parent=$mark/){$GFF{$name}.=$line; }
}
print"1\n";
foreach my$line(keys%GFF){
    @temp=split(/\n/,$GFF{$line}); my$len5=0; my$len3=0; my $mark=0; my$gene_len=0;
    foreach my$line2(@temp){
       if($line2=~/gene\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+);Name/){ $gene_len=$2-$1+1; $gene_name=$3;
       }elsif($line2=~/five_prime_UTR\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+)/){$len5+=$2-$1+1;
       }elsif($line2=~/CDS\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+)/){ 
       }elsif( $line2=~/three_prime_UTR\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+)/){$len3+=$2-$1+1; }
    }
    my$CDS_ratio=($gene_len-$len5-$len3)/$gene_len; $CDS_ratio{$gene_name}=$CDS_ratio;
}
print"2\n";
open(IN2,"<$dir/chrname") or die;
open(OUT1,">$dir/9Merged_RMoverlap_newgene.gff3");
%temp=%gene;
foreach(<IN2>){ 
    if(/^(chr_\d+)\s.+?supercontig\s+(\d+)\t(\d+)\s+/){ 
        my$chr=$1; my$contig_start=$2; my$contig_end=$3;  
        open(TXT,"<$dir/chr/$chr.gff3") or die;
        @Chr_gff=<TXT>; close TXT; 
        foreach my$line(@Chr_gff){ undef@repeat_GENE; 
             if($line=~/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/ && !exists$repeat_line{$line}){ 
               $gene_start=$2; $gene_end=$3; my$name=$5;
               delete$temp{$line}; push@repeat_GENE,$5; $repeat_line{$line}++; undef%repeat_gene;
               foreach my$m(keys%temp){
                 if($m=~/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name=(.+?);/){ my$x=$2; my$y=$3;
                   if( ($gene_start<=$x) and ($gene_end>=$y)  ){ push@repeat_GENE,$5; $repeat_line{$m}++; delete$temp{$m}; 
                   }elsif( ($gene_start<=$x) and ($gene_end<=$y) and ($x<=$gene_end)   ){ push@repeat_GENE,$5; $repeat_line{$m}++; delete$temp{$m}; 
                   }elsif( ($gene_start>=$x) and ($gene_end>=$y) and ($y>=$gene_start) ){ push@repeat_GENE,$5; $repeat_line{$m}++; delete$temp{$m};
                   }elsif( ($gene_start>=$x) and ($gene_end<=$y) ){ push@repeat_GENE,$5; $repeat_line{$m}++;  delete$temp{$m};  }
                 }
               }
               foreach(@repeat_GENE){ $repeat_gene{$_}=$CDS_ratio{$_}; }
               foreach(reverse sort{$repeat_gene{$a}<=>$repeat_gene{$b}}keys%repeat_gene){print OUT1"$GFF{$_}";last; }
            }
        }
    } 
}

