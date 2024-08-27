# perl 01Merge_orf1_to_UTR-gff.pl stringtie.gtf total_transcripts_exon.fasta_orf-1_ATG2TGA
use strict;
my(@GFF,@list);
open (IN1, "<", "$ARGV[0]") or die "Can't open IN1 : $!";
open (IN2, "<", "$ARGV[1]") or die "Can't open IN2 : $!";
open (IN3, "<", "$ARGV[2]") or die "Can't open IN3 : $!";
open (OUT1, ">", "$ARGV[0]_all-UTR.gff3");
open (OUT2, ">", "$ARGV[0]_Not_in_orffinder_mRNA.txt");
{
    @GFF=<IN1>;  close IN1;
    @list=<IN2>; close IN2;
    @Ori=<IN3>; close IN3;
}
foreach(@Ori){ if(){}  }
#STRG.13475.1    ORF6    488     288     -199
my(%list2,%gene,$ori,$n,$mRNA,%mRNA,$ORF, %ORF_list, @ORF_list, %list);
foreach(@list){  if(/(STRG.+?)\s+(ORF\d+)\s+(\d+)\s+(\d+)\s+(.+?)\n/){ $ORF_list{"$1\t$2"}=$1; if($4>$3){$ori="+";  $list{$1}="$3\t$4\t$ori\t$2";  $list2{"$1\t$2"}="$3\t$4\t$ori\t$2";}#else{$ori="-"; $list{$1}="$4\t$3\t$ori\t$2";  $list2{"$1\t$2"}="$4\t$3\t$ori\t$2";} } 
              }
$n=keys%list2;  print "$n\n";

my($transcript_id, $transcrpt_ori, $orf_ori);
foreach(@GFF){	if(/transcript.+?transcript_id\s+"(.+?)";.+?FPKM "(.+?)";/){  $transcript_id=$1;  }elsif(/exon.+?transcript_id\s+"($transcript_id)";/){	$mRNA{$transcript_id}.=$_; }  }
$n=keys%mRNA;  print "$n\n";

my($chr,$gene_start,$gene_end,$mRNA_id,$gene_id,@exon,@exon_coord,$mark,@CDS_coord,$L1,$L2,$i,@n,@m,$start,$end,$length);
foreach(@GFF){
	if(/(chr_\d+)\s.+?transcript\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?gene_id\s+"(.+?)";\s+transcript_id\s+"(.+?)";.+?FPKM "(.+?)";/){
         $chr=$1; $gene_start=$2; $gene_end=$3; $mRNA_id=$6; $gene_id=$5; $transcrpt_ori=$4;
         if(exists $list{$mRNA_id}){  # print "\n$mRNA_id\n";
           @exon=split/\n/,$mRNA{$mRNA_id};
           undef @exon_coord; undef@CDS_coord; undef$length; 
           foreach(@exon){ if(/(chr_\d+)\s.+?exon\s+(\d+)\s+(\d+)\s.+?transcript_id\s+"($mRNA_id)";/){
           		         push(@exon_coord,($2,$3)); $length=$length+$3-$2+1;
           		         push(@CDS_coord,($2,$3));}
            } my@temp_CDS=@CDS_coord;
            #if($list2{$mRNA_id}=~/(\d+)\t(\d+)\t(\+|-)\t(.+)/){ $L1=$1; $L2=$2; $orf_ori=$3; $ORF=$4; }
            undef@ORF_list; foreach(keys%ORF_list){ if($ORF_list{$_} eq $mRNA_id){push@ORF_list,$list2{$_}; } } $n=@ORF_list; print"$n\t";
            foreach my$x(@ORF_list){
            if($x=~/(\d+)\t(\d+)\t(\+|-)\t(.+)/){ $L1=$1; $L2=$2; $orf_ori=$3; $ORF=$4; } 
            if($transcrpt_ori eq "+" && $orf_ori eq "+"){ $ori="+"; }elsif($transcrpt_ori eq "+" && $orf_ori eq "-"){ $ori="-"; 
            }elsif($transcrpt_ori eq "-" && $orf_ori eq "+" ){ $ori="-"; }elsif($transcrpt_ori eq "-" && $orf_ori eq "-"){ $ori="+"; my$L3=$2;; $L2=$length-$L1-1;  $L1=$length-$L3-1;}
            if($transcrpt_ori eq "."){$ori=$orf_ori; } 
            $i=0; undef@n; undef @m; $mark=0; @CDS_coord=@temp_CDS;
            if($ori eq "-"){ my$a=$L1; my$b=$L2; $L1=$length-$b-1; $L2=$length-$a-1; } 
            while($i<@exon_coord){ #print "i\t$i\n";
            	$start=$exon_coord[$i]+$L1; #print "start\t$start\n";
            	if($exon_coord[$i+1]>$start ){ 
            		if($L1>0){ push @m,$exon_coord[$i],$start-1;
            		           shift @CDS_coord; unshift @CDS_coord,$start;                         
                                   $L1=0; #print "m:\t@m\n CDS_corrd:\t@CDS_coord\n";
                    }
            	}else{ push @m,$exon_coord[$i],$exon_coord[$i+1];
            		   shift @CDS_coord;
            		   shift @CDS_coord;
            		   $L1=$L1-($exon_coord[$i+1]-$exon_coord[$i])-1;# print "m:\t@m\n CDS_corrd:\t@CDS_coord\n";
            	}
            	$end=$exon_coord[$i]+$L2; #print "end\t$end\n";
            	if($exon_coord[$i+1]>$end ){ 
            		if($L2>0){ push @n,$end+1,$exon_coord[$i+1];
            		           pop @CDS_coord; push @CDS_coord,$end;                         
                                   $L2=0; #print "n:\t@n\n CDS_corrd:\t@CDS_coord\n";
                    }elsif($L2<=0){ push @n,$exon_coord[$i],$exon_coord[$i+1];
                    	            pop @CDS_coord;
                    	            pop @CDS_coord; #print "n:\t@n\n CDS_corrd:\t@CDS_coord\n";
                    }
            	}elsif($exon_coord[$i+1]==$end){ $L2=0; $mark=1;
                }else{   my$len=$exon_coord[$i+1]-$exon_coord[$i]+1; 
            		$L2=$L2-($exon_coord[$i+1]-$exon_coord[$i])-1; #print "len\t$len\tL2\t$L2\n";
            	}
            	$i+=2;
            } 
            if($mark==0){pop @CDS_coord; push @CDS_coord,$n[0]-1;}
            if($ori eq "+"){
                print OUT1 "$chr\tAUGUSTUS\tgene\t$gene_start\t$gene_end\t.\t$ori\t.\tID=${mRNA_id}_$ORF;Name=${mRNA_id}_$ORF;Note=\n";
                print OUT1 "$chr\tAUGUSTUS\tmRNA\t$gene_start\t$gene_end\t.\t$ori\t.\tID=${mRNA_id}_$ORF;Parent=${mRNA_id}_$ORF\n";
                $i=0; while ( $i < @m ) { print OUT1 "$chr\tAUGUSTUS\tfive_prime_UTR\t$m[$i]\t$m[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.utr5;Parent=${mRNA_id}_$ORF\n"; $i+=2; }
                $i=0; while ($i<@CDS_coord) {   print OUT1 "$chr\tAUGUSTUS\tCDS\t$CDS_coord[$i]\t$CDS_coord[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.cds;Parent=${mRNA_id}_$ORF\n"; $i+=2; } 
                $i=0; while ( $i < @n ) { print OUT1 "$chr\tAUGUSTUS\tthree_prime_UTR\t$n[$i]\t$n[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.utr3;Parent=${mRNA_id}_$ORF\n"; $i+=2; }
            }elsif($ori eq "-"){
            	print OUT1 "$chr\tAUGUSTUS\tgene\t$gene_start\t$gene_end\t.\t$ori\t.\tID=${mRNA_id}_$ORF;Name=${mRNA_id}_$ORF;Note=\n";
                print OUT1 "$chr\tAUGUSTUS\tmRNA\t$gene_start\t$gene_end\t.\t$ori\t.\tID=${mRNA_id}_$ORF;mRNA_id=${mRNA_id}_$ORF\n";
            	$i=0; while ( $i < @m ) { print OUT1 "$chr\tAUGUSTUS\tthree_prime_UTR\t$m[$i]\t$m[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.utr3;Parent=${mRNA_id}_$ORF\n"; $i+=2; }
                $i=0; while ($i<@CDS_coord) {  print OUT1 "$chr\tAUGUSTUS\tCDS\t$CDS_coord[$i]\t$CDS_coord[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.cds;Parent=${mRNA_id}_$ORF\n"; $i+=2; } 
                $i=0; while ( $i < @n ) { print OUT1 "$chr\tAUGUSTUS\tfive_prime_UTR\t$n[$i]\t$n[$i+1]\t.\t$ori\t.\tID=${mRNA_id}_$ORF.t1.utr5;Parent=${mRNA_id}_$ORF\n"; $i+=2; }
            }
        }#else{print OUT2 "$mRNA_id\n";}
        }
    }
}
