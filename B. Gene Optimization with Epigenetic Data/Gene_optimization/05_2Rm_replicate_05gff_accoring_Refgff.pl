use strict;
my$dir=$ENV{'PWD'};
my(%gff, @gff, $gene_id, $gene_line, $chr, %gene_id, %min_UTR, %gff_ref);
open(IN1,"<$dir/05higest_STRG_loci.gtf_UTR.gff3");
open(IN2,"<$dir/06Merge_Deleted636-81_UTR2901-550_To13gff3.gff3");
#-----------------------------------------------------------------
print"STEP1: STart import gene imformation from 05higest_STRG_loci.gtf_UTR.gff3\n";
foreach(<IN1>){
     if(/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+.+?ID=(chr_\d+);/){ $gff{$3}=$_; push@gff,$_;
       }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){  $gene_id=$5; $gene_line=$_; $gff{$gene_line}=$_; $chr=$1; push@gff,$_;
       }elsif(/$chr.+?mRNA.+?Parent=$gene_id/){  $gff{$gene_line}.=$_;
       }elsif(/$chr.+?UTR.+?Parent=$gene_id/){   $gff{$gene_line}.=$_;
       }elsif(/$chr.+?CDS.+?Parent=$gene_id/){   $gff{$gene_line}.=$_;
       }elsif(/chr.+?rDNA.+?ID=(.+?);Name/){     $gff{$_}=$_; push@gff,$_;
       }
}
my $n=keys%gff;
print"$n gene with duplicate were importted\n";
#----------------------------------------------------------------
my($gene_start, $gene_end, @temp, $number, %UTR_min, %min, %min_genelen, %gff_rm, %gff_id);
print"STEP2: Start caculate best gene model by their UTR number and gene length\n";
open(OUT1,">$dir/05-2Rm_duplicate_Merge.gff3");
open(OUT2,">$dir/05-2Rm_duplicate.gff3");
foreach(@gff){
    if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){
       $gene_id=$5; $gene_line=$_; $gene_start=$2; $gene_end=$3; $chr=$1; 
       if(!exists$gene_id{$5}){ undef@temp; $number=0; undef %min_UTR; undef %min; undef %min_genelen;
         $min_UTR{$gene_id}=100;
         foreach my$i(keys%gff){
            if($i=~/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){
               if($gene_start-$3<0){$gene_id{$5}=$5; push@temp,$i;}                       
            }
         }
         foreach my$j(@temp){
            if($j=~/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){
                  $number=$gff{$j}=~s/UTR/UTR/g; 
                  if($number<$min_UTR{$gene_id}){$min_UTR{$gene_id}=$number; $min_genelen{$gene_id}=$3-$2; $min{$gene_id}=$j; #print "$gene_id\n"; 
                 }elsif($number==$min_UTR{$gene_id} && $min_genelen{$gene_id}<($3-$2)){$min_UTR{$gene_id}=$number; $min_genelen{$gene_id}=$3-$2; $min{$gene_id}=$j; } 
            }   
         }
         $gff_rm{$min{$gene_id}}=0;
       }
   }
}

foreach(keys%gff_rm){
    print OUT2"$gff{$_}";
}

my$n=keys%gff_rm;
print"$n gene were remove deplicate.\n";
#----------------------------------------------------------------------
undef@gff; my$n=0;
foreach(<IN2>){
     if(/(chr_\d+)\s+.+?supercontig\s+(\d+)\s+.+?ID=(chr_\d+);/){ $gff_ref{$3}=$_; push@gff,$_;
       }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){  $gene_id=$5; $gene_line=$_; $gff_ref{$gene_line}=$_; $chr=$1; push@gff,$_; $n++;
       }elsif(/$chr.+?mRNA.+?Parent=$gene_id/){  $gff_ref{$gene_line}.=$_;
       }elsif(/$chr.+?UTR.+?Parent=$gene_id/){   $gff_ref{$gene_line}.=$_;
       }elsif(/$chr.+?CDS.+?Parent=$gene_id/){   $gff_ref{$gene_line}.=$_;
       }elsif(/chr.+?rDNA.+?ID=(.+?);Name/){     $gff_ref{$_}=$_; push@gff,$_;
       }
}
print "$n gene imported\n";
$n=0; my$m=0; my$mark=3;
foreach my$i(@gff){
     if($i=~/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){
          $gene_start=$2; $chr=$1; $m++; $mark=0;
        #  print OUT1 "$gff_ref{$i}" if keys%gff_rm==0;
          foreach my$j(keys%gff_rm){ 
             if($j=~/($chr)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);Name/){
                if($2<$gene_start && $gff_rm{$j}==0){ 
                      print OUT1"$gff{$j}$gff_ref{$i}"; $gff_rm{$j}=1; $n++;  $mark=1; 
                } #else{print OUT1"$gff_ref{$i}"; $n++; }
             } 
          }
          print OUT1"$gff_ref{$i}" if$mark==0;
     }else{ print OUT1"$gff_ref{$i}"; $m++; $n++; }
}
print "n: $n gene \n m: $m\n";

