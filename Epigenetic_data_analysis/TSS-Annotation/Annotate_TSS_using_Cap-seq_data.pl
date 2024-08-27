$dir=$ENV{'PWD'};
open(GENE,"<","$dir/TGD2024_renamed.gff3");
open(IN2,"<","$dir/start_TSS.bed");
open(IN3,"<","$dir/TGD2024_renamed_addAS-NAT-NC.gff3");
#chr_001 Manual  gene    256     1182    .       +       .       ID=g18410;Name=TTHERM_00161861;UserNamed=Unnamed;TGDAnno=Unknow;info=hypothetical protein;
#chr_001 Manual  mRNA    256     1182    .       +       .       ID=g18410.t1;Parent=g18410;
#chr_001 Manual  five_prime_UTR  256     333     .       +       .       ID=g18410.t1.utr;Parent=g18410.t1;
#chr_001 Manual  CDS     334     342     .       +       .       ID=g18410.t1.cds;Parent=g18410.t1;
#chr_001 Manual  CDS     394     453     .       +       .       ID=g18410.t1.cds;Parent=g18410.t1;
#chr_001 Manual  CDS     505     855     .       +       .       ID=g18410.t1.cds;Parent=g18410.t1;
#chr_001 Manual  three_prime_UTR 856     1182    .       +       .       ID=g18410.t1.utr;Parent=g18410.t1;
foreach$line(<IN3>){
   if($line=~/^(chr_\d+)\s+Manual\s+gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?);Name=(.+?);/){
     if(grep /^$1$/,@chr){  }else{ push@chr,$1;  } $Name_gene{$Tname}++;
     $chr=$1; $start=$2; $end=$3; $ori=$4; $gname=$5; $Tname=$6; $Tname2chr{$Tname}=$chr;
     push@{$sorted_name{$chr}},$Tname; $Tname2gname{$Tname}=$gname; $ori{$Tname}=$ori;
     $mid_gene{$Tname}=int(($start+$end)/2); $start_gene{$Tname}=$start; $end_gene{$Tname}=$end;
     push@{$gff{$chr}{$Tname}{"gene"}},$line;
   }elsif($line=~/^(chr_\d+)\s+Manual\s+mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?.t1);Par/){
     push@{$gff{$chr}{$Tname}{"mRNA"}},$line;
   }elsif($line=~/^(chr_\d+)\s+Manual\s+(CDS)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?.);Par.+?t1/){
     $region=$2; push@{$gff{$chr}{$Tname}{"CDS"}},$line; push@{$cds_pos{$chr}{$Tname}{"CDS"}},$3,$4; push@{$exon_pos{$chr}{$Tname}{"exon"}},$3,$4;
   }elsif($line=~/^(chr_\d+)\s+Manual\s+(three.+?)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?.);Par.+?t1/){
     $region=$2; push@{$gff{$chr}{$Tname}{"three"}},$line; push@{$exon_pos{$chr}{$Tname}{"three"}},$3,$4; push@{$exon_pos{$chr}{$Tname}{"exon"}},$3,$4;
   }elsif($line=~/^(chr_\d+)\s+Manual\s+(five.+?)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?.);Par.+?t1/){
     $region=$2; push@{$gff{$chr}{$Tname}{"five"}},$line; push@{$exon_pos{$chr}{$Tname}{"five"}},$3,$4; push@{$exon_pos{$chr}{$Tname}{"exon"}},$3,$4;
   }elsif($line=~/^(chr_\d+)\s+Manual\s+(.+?)\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+ID=(.+?.);Par.+?(t\d+);/ && $7 ne "t1"){
     push@{$gff{$chr}{$Tname}{"ASNC"}},$line;
   }
}
$n=@chr;
@line=@{$gff{"chr_001"}{"TTHERM_00161579"}{"five"}};
#print "@line";
#print "\nchr number:$chr\n";
#chr_001 247     527
#chr_001 2017    2324
foreach(<IN2>){
   if(/^(chr_\d+)\s+(\d+)\s+(\d+)\s/){
       $chr=$1; $start=$2; $end=$3; $mid=($start+$end)/2;
       $TSS{$1}{start}=$2; $TSS{$1}{end}=$3;
       undef$mark;
       foreach $Tname(keys%{$gff{$chr}}){
           @temp=@{$gff{$chr}{$Tname}{"gene"}};  
           @CDS=@{$cds_pos{$chr}{$Tname}{"CDS"}}; 
           $mark=0;
           if($temp[0]=~/^(chr_\d+)\s+Manual\s+gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.\s+(ID=(.+?);.+)\n/ && exists$Name_gene{$Tname}){
               $chr=$1; $gene_start=$2; $gene_end=$3; $ori=$4; $Info=$5; $id=$6;
               if($ori eq "+"  && $end>$gene_start && $end<$gene_end && $start<$CDS[0]){ $TSS=$start;  
                   $line="$chr\tManual\tgene\t$TSS\t$gene_end\t.\t$ori\t.\t$Info\n"; undef@{$gff{$chr}{$Tname}{"gene"}};  push@{$gff{$chr}{$Tname}{"gene"}},$line;
                   $line="$chr\tManual\tmRNA\t$TSS\t$gene_end\t.\t$ori\t.\tID=$id.t1;Parent=$id;\n"; undef@{$gff{$chr}{$Tname}{"mRNA"}};  push@{$gff{$chr}{$Tname}{"mRNA"}},$line;
                   @exon=@{$exon_pos{$chr}{$Tname}{"exon"}};
                   @UTR3=@{$exon_pos{$chr}{$Tname}{"three"}}; undef@{$gff{$chr}{$Tname}{"three"}};
                   @UTR5=@{$exon_pos{$chr}{$Tname}{"five"}};  undef@{$gff{$chr}{$Tname}{"five"}};
                   if(@UTR5>0){ $UTR5[0]=$TSS; }else{ $CDS[0]=$TSS; } 
                   undef@{$gff{$chr}{$Tname}{"CDS"}};
                   $i=0; while($i<@UTR5){ $line="$chr\tManual\tfive_prime_UTR\t$UTR5[$i]\t$UTR5[$i+1]\t.\t$ori\t.\tID=$id.t1.utr;Parent=$id.t1;\n"; push@{$gff{$chr}{$Tname}{"five"}},$line;   $i+=2; }
                   $i=0; while($i<@CDS){  $line="$chr\tManual\tCDS\t$CDS[$i]\t$CDS[$i+1]\t.\t$ori\t.\tID=$id.t1.cds;Parent=$id.t1;\n";   push@{$gff{$chr}{$Tname}{"CDS"}},$line;    $i+=2; }
                   $i=0; while($i<@UTR3){ $line="$chr\tManual\tthree_prime_UTR\t$UTR3[$i]\t$UTR3[$i+1]\t.\t$ori\t.\tID=$id.t1.utr;Parent=$id.t1;\n"; push@{$gff{$chr}{$Tname}{"three"}},$line;  $i+=2; }
                   delete$Name_gene{$Tname};
               }elsif($ori eq "-" && $start>$gene_start && $start<$gene_end && $end>$CDS[-1]){ $TSS=$end; 
                   $line="$chr\tManual\tgene\t$gene_start\t$TSS\t.\t$ori\t.\t$Info\n"; undef@{$gff{$chr}{$Tname}{"gene"}};  push@{$gff{$chr}{$Tname}{"gene"}},$line;
                   $line="$chr\tManual\tmRNA\t$gene_start\t$TSS\t.\t$ori\t.\tID=$id.t1;Parent=$id;\n"; undef@{$gff{$chr}{$Tname}{"mRNA"}};  push@{$gff{$chr}{$Tname}{"mRNA"}},$line;
                   @exon=@{$exon_pos{$chr}{$Tname}{"exon"}}; 
                   @UTR3=@{$exon_pos{$chr}{$Tname}{"three"}}; undef@{$gff{$chr}{$Tname}{"three"}};
                   @UTR5=@{$exon_pos{$chr}{$Tname}{"five"}};  undef@{$gff{$chr}{$Tname}{"five"}};
                   if(@UTR5>0){ $UTR5[-1]=$TSS; }else{ $CDS[-1]=$TSS; }
                   undef@{$gff{$chr}{$Tname}{"CDS"}};
                   $i=0; while($i<@UTR3){ $line="$chr\tManual\tthree_prime_UTR\t$UTR3[$i]\t$UTR3[$i+1]\t.\t$ori\t.\tID=$id.t1.utr;Parent=$id.t1;\n"; push@{$gff{$chr}{$Tname}{"three"}},$line;  $i+=2;}
                   $i=0; while($i<@CDS){  $line="$chr\tManual\tCDS\t$CDS[$i]\t$CDS[$i+1]\t.\t$ori\t.\tID=$id.t1.cds;Parent=$id.t1;\n";   push@{$gff{$chr}{$Tname}{"CDS"}},$line;   $i+=2;}
                   $i=0; while($i<@UTR5){ $line="$chr\tManual\tfive_prime_UTR\t$UTR5[$i]\t$UTR5[$i+1]\t.\t$ori\t.\tID=$id.t1.utr;Parent=$id.t1;\n"; push@{$gff{$chr}{$Tname}{"five"}},$line; $i+=2;}
                   delete$Name_gene{$Tname};
               }
           }
       }
   }
}
@line=@{$gff{"chr_001"}{"TTHERM_00161579"}{"five"}};
#print "@line";
open(OUT1,">","$dir/TGD2024_renamed_TSSanno_addAS-NAT-NC.gff3");
foreach $chr(@chr){
    foreach $Tname(@{$sorted_name{$chr}}){
       if(!exists$Name_gene{$Tname}){
           #print "$Tname-----\n";
           @gene=@{$gff{$chr}{$Tname}{"gene"}}; print OUT1 "$gene[0]";
           @mRNA=@{$gff{$chr}{$Tname}{"mRNA"}}; print OUT1 "$mRNA[0]";
           if($ori{$Tname} eq "+"){
             @UTR5=@{$gff{$chr}{$Tname}{"five"}}; foreach $line(@UTR5){ print OUT1 "$line"; }   
             @CDS=@{$gff{$chr}{$Tname}{"CDS"}}; foreach $line(@CDS){ print OUT1 "$line"; }   
             @UTR3=@{$gff{$chr}{$Tname}{"three"}}; foreach $line(@UTR3){ print OUT1 "$line"; }      
           }elsif($ori{$Tname} eq "-"){
             @UTR3=@{$gff{$chr}{$Tname}{"three"}}; foreach $line(@UTR3){ print OUT1 "$line"; }
             @CDS=@{$gff{$chr}{$Tname}{"CDS"}}; foreach $line(@CDS){ print OUT1 "$line"; }
             @UTR5=@{$gff{$chr}{$Tname}{"five"}}; foreach $line(@UTR5){ print OUT1 "$line"; }
           }
           #print "\n";
       }else{
           #print "$Tname-----\n";
           @gene=@{$gff{$chr}{$Tname}{"gene"}}; print OUT1 "$gene[0]";
           @mRNA=@{$gff{$chr}{$Tname}{"mRNA"}}; print OUT1 "$mRNA[0]";
          # @exon=@{$gff{$chr}{$Tname}{"exon"}}; foreach $line(@exon){ print "$line"; }
           if($ori{$Tname} eq "+"){
             @UTR5=@{$gff{$chr}{$Tname}{"five"}}; foreach $line(@UTR5){ print OUT1 "$line"; }   
             @CDS=@{$gff{$chr}{$Tname}{"CDS"}}; foreach $line(@CDS){ print OUT1 "$line"; }   
             @UTR3=@{$gff{$chr}{$Tname}{"three"}}; foreach $line(@UTR3){ print OUT1 "$line"; }      
           }elsif($ori{$Tname} eq "-"){
             @UTR3=@{$gff{$chr}{$Tname}{"three"}}; foreach $line(@UTR3){ print OUT1 "$line"; }
             @CDS=@{$gff{$chr}{$Tname}{"CDS"}}; foreach $line(@CDS){ print OUT1 "$line"; }
             @UTR5=@{$gff{$chr}{$Tname}{"five"}}; foreach $line(@UTR5){ print OUT1 "$line"; }
           }
           #print "\n";
       }
       @ASNC=@{$gff{$chr}{$Tname}{"ASNC"}}; foreach $line(@ASNC){ print OUT1 "$line"; }
    }
}

foreach(<IN3>){
   if(//){

   }
}

