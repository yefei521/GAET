#!/usr/bin/perl -w
use threads;
use File::Spec;
use threads::shared;
use Fcntl qw(:flock);
$dir=$ENV{'PWD'};
open(IN2,"<$dir/chr.txt"); #chr.txt
open(IN1,"<","$dir/Manual_check-total-gene.gff3_Right_UTR.gff3");
@IN2=<IN2>; close IN2;
@IN1=<IN1>; close IN1;
#chr_001 226     854    g18410  .   +
foreach(@IN2){
   if(/^(chr_\d+)\s+/){ $chr{$1}=0; }
}
foreach(@IN1){
   if(/^(chr_\d+)\s+.+?gene\s+(\d+)\s+(\d+)\s+.+?\s+(.)\s+.+?\s+ID=(.+?);/){
        $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5;
        $gff{$chr}{$name}="$start\t$end";
        $ori{$chr}{$name}=$ori;
        $Gene_len{$name}=$end-$start+1;
   }
}

#$C=200737460; $I=174763830;
#---------------------------
print "START.....\n";
@chr=keys%chr;
foreach my $i ( 1 .. 180){
    $p=shift@chr;  $Inputp=$p;
    if( fork() ){  wait if($i+1 >= 20);
    }else{
       $hash2=&MATCH_INPUT($Inputp,$dir); #%sigmaC=%{$hash2};
       print "END:threads$p\n";
       exit();
    }
    sleep 0.1;
}

until(wait() == -1){
 $cmd=`cat $dir/ReadsCount/*txt > ReadsCount.txt`;
 print $cmd;
 $cmd=`echo "Name\tpoly-A_len\tName\tReadsCount\tGene-len\n" > ManUTR+lowQ_gene_polyA_length.txt_exp`;
 print $cmd;
 $cmd=`awk 'NR==FNR{a[$1]=$0;} NR>FNR{match($1,/(.+)\.t/,b); if(b[1] in a ){print $1,$4,a[b[1]];}}' ReadsCount.txt ManUTR+lowQ_gene_polyA_length.txt > ManUTR+lowQ_gene_polyA_length.txt_exp`;
 print $cmd;
}

#------------------------
sub MATCH_INPUT{
   my($zmwp,$dir)=@_; print "START: threads:$zmwp\n";
   my$zmwp_rep=$zmwp; $zmwp_in=$zmwp."_in"; $zmwp_out=$zmwp."_out"; $chr=$zmwp;
   open($zmwp_in,"<$dir/WT_mix_seqkit2DNA_MaxIntron2k/$zmwp") or die;
   open($zmwp_out,">$dir/ReadsCount/$zmwp.txt") or die;
   while(<$zmwp_in>){ 
       if(/^(.+?)\s+(0)\s+$chr\s+(\d+)\s+.+?\s+(.+?)\s+.+?\s+.+?\s+.+?\s+(.+?)\s/){ #6167d4cf-9149-4179-8de5-258a6d61c8b5    0       chr_133 432544  60      27M2I4M1I29M5D5M1D2M1D24M1I48M2D53M2I43M1D21M20S        *       0       0       TGCCCTTATGCCCACAATTAAGTATGCATATTCAAAAG
           $id=$1; $flag=$2; $start=$3; $match=$4; $sequence=$5;
           $temp=$match; $match_len=0; undef $poly_seq_filter;
           while($temp=~/^(\d+)([MDSIHN])(.+)/){  $temp=$3; $MM=$2; $l=$1;  if($MM=~/M|D|N/){$match_len+=$l;}  }
           $mid=$start+int($match_len/2);
           if($flag==0){$ori="+"; }elsif($flag==16){$ori="-"; }
           foreach $name(keys%{$gff{$chr}}){
              next if $ori{$chr}{$name} ne $ori;
              if($gff{$chr}{$name}=~/^(\d+)\s+(\d+)/ && $mid>=$1 && $mid<=$2){ #print "$chr\t$name\t$1\t$2\n";
                   $Count{$name}++;
              }
           }
       }elsif(/^(.+?)\s+(16)\s+$chr\s+(\d+)\s+.+?\s+(.+?)\s+.+?\s+.+?\s+.+?\s+(.+?)\s/){
           $id=$1; $flag=$2; $start=$3; $match=$4; $sequence=$5;
           $temp=$match; $match_len=0;
           while($temp=~/^(\d+)([MDSIHN])(.+)/){  $temp=$3; $MM=$2; $l=$1;  if($MM=~/M|D|N/){$match_len+=$l;}  }
           $mid=$start+int($match_len/2);
           if($flag==0){$ori="+"; }elsif($flag==16){$ori="-"; }
           foreach $name(keys%{$gff{$chr}}){
              next if $ori{$chr}{$name} ne $ori;
              if($gff{$chr}{$name}=~/^(\d+)\s+(\d+)/ && $mid>=$1 && $mid<=$2){ #print "$chr\t$name\t$1\t$2\n";
                   $Count{$name}++; #last;
              }
           }
      }
   }
   foreach(keys%Count){
        print $zmwp_out "$_\t$Count{$_}\t$Gene_len{$_}\n";
   }
   $hash=\%hash1;
   return($hash);
}


