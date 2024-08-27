#!/usr/bin/perl -w
use threads;
use File::Spec;
use threads::shared;
use Fcntl qw(:flock);
open IN1, "<$ARGV[0]"; #input.sam
open(IN2,"<$ARGV[1]"); #chr.txt
@IN2=<IN2>; close IN2;
$dir=$ENV{'PWD'};
#chr_001 226     854    g18410  .   +
foreach(@IN2){
   if(/^(chr_\d+)\s+/){ $chr{$1}=0; }
}
print "MKDIR directory\n\n";
$ARGV[0]=~/^(\w+).sam/; $sample=$1;
$cmd=`rm -r $dir/$sample`; print $cmd;
$cmd=`mkdir -p $dir/$sample`; print $cmd;
$cmd=`rm -r $dir/temp`; print $cmd;
$cmd=`mkdir -p $dir/temp`; print $cmd;

foreach my$i(keys%chr){ $Inputi=$i;
   open($Inputi,">","$dir/$sample/$Inputi") or die;
}
print "START Split bed to 181 chr\n";
while(<IN1>){ if(/(chr_\d+)\s+(.+?)\s/){ $chr=$1; $Inputi=$chr;  print $Inputi "$_" if fileno($Inputi); $I++;} }

foreach my$i(keys%chr){ $Inputi=$i;  close  $Inputi; }
print "CLOSE IO\n";
#$C=200737460; $I=174763830;
#---------------------------
print "START.....\n";
@chr=keys%chr;
foreach my $i ( 1 .. 180){
    $p=shift@chr;  $Inputp=$p;
    if( fork() ){  wait if($i+1 >= 20);
    }else{
       $hash2=&MATCH_INPUT($Inputp,$dir,$sample); #%sigmaC=%{$hash2};
       print "END:threads$p\n";
       exit();
    }
    sleep 0.1;
}

until(wait() == -1){
 $cmd=`cat $dir/temp/*txt > $ARGV[0]_polya.txt`;
 print $cmd;
}

#------------------------
sub MATCH_INPUT{
   my($zmwp,$dir,$sample)=@_; print "START: threads:$zmwp\n";
   my$zmwp_rep=$zmwp; $zmwp_in=$zmwp."_in"; $zmwp_out=$zmwp."_out"; $chr=$zmwp;
   open($zmwp_in,"<$dir/$sample/$zmwp") or die;
   open($zmwp_out,">$dir/temp/$zmwp.txt") or die;
   while(<$zmwp_in>){ 
       if(/^(.+?)\s+(0)\s+$chr\s+(\d+)\s+.+?\s+(.+?)\s+.+?\s+.+?\s+.+?\s+(.+?)\s/){ #6167d4cf-9149-4179-8de5-258a6d61c8b5    0       chr_133 432544  60      27M2I4M1I29M5D5M1D2M1D24M1I48M2D53M2I43M1D21M20S        *       0       0       TGCCCTTATGCCCACAATTAAGTATGCATATTCAAAAG
           $id=$1; $flag=$2; $start=$3; $match=$4; $sequence=$5;
           $temp=$match; $match_len=0; undef $poly_seq_filter;
           if($match=~/(\d+)S$/){ $clip=$1; }else{$clip=99;}
           next if $clip<3;
           while($temp=~/^(\d+)([MDSIHN])(.+)/){  $temp=$3; $MM=$2; $l=$1;  if($MM=~/M|D|N/){$match_len+=$l;}  }
           $polya_seq=substr($sequence,-$clip);
           $polya_len= $polya_seq =~tr/A/A/;
           next if $polya_len/length($polya_seq)<0.5;
           $end=$start+$match_len;
           #$polya_seq =~/A{3,}|()(.+)/;
           if(length($polya_seq)>15){$max_non_a = int(length($polya_seq) * 0.4);}else{ $max_non_a =3; }
           @matches = $polya_seq =~ /(A(?:[^A]{0,$max_non_a}A)*)/g;
           $longest_match = "";
           foreach $a(@matches){
                 if (length($a) > length($longest_match)) { $longest_match = $a;}     
           }
           $longest=length($longest_match);
           next if $longest<5;
           $hash1{"$chr\t$start\t$end\t$id\t$longest\t+\t$polya_seq\tpolya:$longest_match\n"}=1;
           print $zmwp_out "$chr\t$start\t$end\t$id\t$longest\t+\t$polya_seq\tpolya:$longest_match\n";
       }elsif(/^(.+?)\s+(16)\s+$chr\s+(\d+)\s+.+?\s+(.+?)\s+.+?\s+.+?\s+.+?\s+(.+?)\s/){
           $id=$1; $flag=$2; $start=$3; $match=$4; $sequence=$5;
           $temp=$match; $match_len=0;
           if($match=~/^(\d+)S/){ $clip=$1; }else{$clip=99;}
           next if $clip<3;
           while($temp=~/^(\d+)([MDNSIH])(.+)/){  $temp=$3; $MM=$2; $l=$1;  if($MM=~/M|D|N/){$match_len+=$l;}  }
           $end=$start+$match_len;
           $polya_seq=substr($sequence,1,$clip);
           $polya_len= $polya_seq =~tr/T/T/;
           $polya_reverse=reverse$polya_seq; $polya_reverse=~tr/atcguATCGU/tagcaTAGCA/;
           next if $polya_len/length($polya_seq)<0.5;
           #print "$id\t$polya_seq\t$polya_reverse\n";
           if(length($polya_seq)>15){$max_non_a = int(length($polya_seq) * 0.4);}else{ $max_non_a =3; }
           @matches = $polya_reverse =~ /(A(?:[^A]{0,$max_non_a}A)*)/g;
           $longest_match = "";
           foreach $a(@matches){
                 if (length($a) > length($longest_match)) { $longest_match = $a;}
           }
           $longest=length($longest_match);
           next if $longest<5;
           $hash1{"$chr\t$start\t$end\t$id\t$longest\t-\t$polya_reverse\tpolya:$longest_match\n"}=1;
           print $zmwp_out "$chr\t$start\t$end\t$id\t$longest\t-\t$polya_reverse\tpolya:$longest_match\n";
      }
   }
   $hash=\%hash1;
   return($hash);
}


