#perl 04Merge_orf1_to_UTR-gff.pl Split_CDS-mRNA-AS-NAT/mRNA_need_coding.gtf Predict_mRNA2CDS_GFF/mRNA_need_coding.fasta_orf1_ATG2TGA_longest
my(@GFF,@list);
open (IN1, "<", "$ARGV[0]") or die "Can't open IN1 : $!";
open (IN2, "<", "$ARGV[1]") or die "Can't open IN2 : $!";
open (OUT1, ">", "$ARGV[1]_UTR.gff3");
open (OUT2, ">", "3Not_in_orffinder_mRNA.txt");
{
    @GFF=<IN1>;  close IN1;
    @list=<IN2>; close IN2;
}

#>WT_C5_1_rmdup_STRG.6611.1      ORF1    335     640     306
my(%list2,%gene,$ori,$n,$mRNA,%mRNA);
foreach(@list){
        if(/^>(.+?)\s+(ORF\d+)\s+(\d+)\s+(\d+)\s+(\d+)\n/){ if($4>$3){$ori="+"; $list2{$1}="$3\t$4\t$ori";} }
}
$n=keys%list2;
print "$n\n";

my($transcript_id,);
foreach(@GFF){
	if(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s.+?ID=(.+?);ManInfo/){  $transcript_id=$4; 
	}elsif(/exon.+?Parent=(.+);/){	$mRNA{$1}.=$_;	}
}
$n=keys%mRNA;
print "$n\n";

my($chr,$gene_start,$gene_end,$mRNA_id,$gene_id,@exon,@exon_coord,$mark,@CDS_coord,$L1,$L2,$i,@n,@m,$start,$end,$length);
foreach(@GFF){
	if(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?);/){
         $chr=$1; $gene_start=$2; $gene_end=$3; $mRNA_id=$5; $ori=$4;
         if(exists $list2{$mRNA_id}){ #print "\n$mRNA_id\n";
           @exon=split/\n/,$mRNA{$mRNA_id};
           undef @exon_coord; undef@CDS_coord; undef$length; 
           foreach(@exon){ if(/(chr_\d+)\s.+?exon\s+(\d+)\s+(\d+)\s.+?Parent=($mRNA_id);/){
           		         push(@exon_coord,($2,$3)); $length=$length+$3-$2+1;
           		         push(@CDS_coord,($2,$3));
                          }else{print "$_\n";}
            }
            @temp=sort{$a<=>$b}@exon_coord; @exon_coord=@temp;
            @temp=sort{$a<=>$b}@CDS_coord; @CDS_coord=@temp;
            #print "gene length:$length\tExon coord:@exon_coord\n";
            if($list2{$mRNA_id}=~/(\d+)\t(\d+)\t(\+|-)/){$L1=$1; $L2=$2; }
            $i=0; undef@n; undef @m; $mark=0;
            if($ori eq "-"){ my$a=$L1; my$b=$L2; $L1=$length-$b-1; $L2=$length-$a-1;} 
            #print "L1:$L1\tL2:$L2\n";
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
            		   $L1=$L1-($exon_coord[$i+1]-$exon_coord[$i])-1; #print "m:\t@m\n CDS_corrd:\t@CDS_coord\n";
            	}
            	$end=$exon_coord[$i]+$L2; #print "end\t$end\n";
            	if($exon_coord[$i+1]>$end ){ 
            		if($L2>0){ push @n,$end+1,$exon_coord[$i+1];
            		           pop @CDS_coord; push @CDS_coord,$end;                         
                                   $L2=0; #print "n:\t@n\n CDS_corrd:\t@CDS_coord\n";
                    }elsif($L2<=0){ push @n,$exon_coord[$i],$exon_coord[$i+1];
                    	            pop @CDS_coord;
                    	            pop @CDS_coord; print "$mRNA_id\n" if $L2==0;
                    }#elsif($L2==0){
                                     
                    #}}
            	}elsif($exon_coord[$i+1]==$end){ $L2=0; $mark=1; 
                }else{   #print "L2\t$L2\n"; 
                         my$len=$exon_coord[$i+1]-$exon_coord[$i]+1; 
            		$L2=$L2-($exon_coord[$i+1]-$exon_coord[$i])-1; #print "len\t$len\tL2\t$L2\n";
                        #if($L2==0){$L2=1; print "len\t$len\tL2\t$L2\n";}
            	}
            	$i+=2;
            } 
            if($mark==0){pop @CDS_coord; push @CDS_coord,$n[0]-1;}
            if($mRNA_id=~/^(.+)\.t/){ $gene_id=$1;}else{ $gene_id=$mRNA_id; $mRNA_id.=".t1"; }
            if($ori eq "+"){
                print OUT1 "$chr\tManual\tgene\t$gene_start\t$gene_end\t.\t$ori\t.\tID=$gene_id;Name=$mRNA_id;Note=\n";
                print OUT1 "$chr\tManual\tmRNA\t$gene_start\t$gene_end\t.\t$ori\t.\tID=$mRNA_id;Parent=$gene_id\n";
                $i=0; while ( $i < @m ) { print OUT1 "$chr\tManual\tfive_prime_UTR\t$m[$i]\t$m[$i+1]\t.\t$ori\t.\tID=${mRNA_id}.utr;Parent=$mRNA_id\n"; $i+=2; }
                $i=0; while ($i<@CDS_coord) {   print OUT1 "$chr\tManual\tCDS\t$CDS_coord[$i]\t$CDS_coord[$i+1]\t.\t$ori\t.\tID=${mRNA_id}.cds;Parent=$mRNA_id\n"; $i+=2; } 
                $i=0; while ( $i < @n ) { print OUT1 "$chr\tManual\tthree_prime_UTR\t$n[$i]\t$n[$i+1]\t.\t$ori\t.\tID=${mRNA_id}.utr;Parent=$mRNA_id\n"; $i+=2; }
            }elsif($ori eq "-"){
            	print OUT1 "$chr\tManual\tgene\t$gene_start\t$gene_end\t.\t$ori\t.\tID=$gene_id;Name=$mRNA_id;Note=\n";
                print OUT1 "$chr\tManual\tmRNA\t$gene_start\t$gene_end\t.\t$ori\t.\tID=$mRNA_id;Parent=$gene_id\n";
            	$i=0; while ( $i < @m ) { print OUT1 "$chr\tManual\tthree_prime_UTR\t$m[$i]\t$m[$i+1]\t.\t$ori\t.\tID=${mRNA_id}.utr;Parent=$mRNA_id\n"; $i+=2; }
                $i=0; while ($i<@CDS_coord) {  print OUT1 "$chr\tManual\tCDS\t$CDS_coord[$i]\t$CDS_coord[$i+1]\t.\t$ori\t.\tID=${mRNA_id}.cds;Parent=$mRNA_id\n"; $i+=2; } 
                $i=0; while ( $i < @n ) { print OUT1 "$chr\tManual\tfive_prime_UTR\t$n[$i]\t$n[$i+1]\t.\t$ori\t.\tID=${mRNA_id}.utr;Parent=$mRNA_id\n"; $i+=2; }
            }
        }else{print OUT2 "$mRNA_id\n";}
    }
}
