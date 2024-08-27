$dir=$ENV{'PWD'};
open(IN1,"<","$dir/Rename_gene/TGD2024_need_rename_gene.gff3");
open(IN2,"<","$dir/Rename_gene/2-upd-Genome-GFF3-latest-2.gff3");
open(IN3,"<","$dir/Rename_gene/TGD2024blastTGD2021_otf7");
open(IN4,"<","$dir/Rename_gene/TGD2024_duplicated_nedd_check.txt_checked_cluster");
open(IN5,"<","$dir/LowQuality/TGD_database_wiki_named_gene.txt");
open(IN6,"<","$dir/Interpro_anno/TGD2024_need_rename_gene_pep.fasta_annotation.txt");
open(IN7,"<","$dir/LowQuality/Low_Quality_gene_GENE.fasta_orf1_ATG2TGA_longest_UTR.gff3");
open(IN8,"<","$dir/LowQuality/Low_Quality_gene_LQ.gff3");
foreach(<IN5>){ @temp=split(/\t/); $Tname=$temp[0]; $Username{$Tname}=$temp[1]; if($temp[2]=~/info:;/){}elsif($temp[2]=~/info:(.+);/){$Anno{$Tname}=$1; }  }
foreach(<IN6>){ chomp; @temp=split(/\t/); $name=$temp[0]; $database=$temp[1]; $annotation=$temp[2]; push@{$InterPro{$name}{$database}},$annotation;  }
foreach $name(keys%InterPro){ undef $line; 
         foreach $database(keys%{$InterPro{$name}}){ @temp=@{$InterPro{$name}{$database}}; @uniq=grep{++$ha{$_}<2}@temp;  $line.="DB-$database:".join(";",@temp).";";  } 
         $InterPro_Anno{$name}=$line; 
}
foreach(<IN4>){
      @temp=split(/\t/);
      foreach$line(@temp){ 
          if($line=~/^(.+?):gene/){ $name=$1; $Rmdup_2024gene{$1}++;  
          }elsif($line=~/^(.+?):rename/){ $New_Tname{$name}=$1; $Named{$name}=$1; $Used_Tname{$1}++; }elsif($line=~/^(.+?):(.+)/){ $NeedDeleteGene{$1}=$2;   }else{ print "$line\n";  }  
       }
}
$n=keys%Used_Tname;
print "After Input TGD2024_duplicated_nedd_check.txt_checked_cluster.\nThe Used name number:$n\n";
undef @temp; undef$name;
#IN1
#chr_065 Manual  gene    1738186 1741894 .       -       .       ID=chr_065-GSAman22373;Name=chr_065-GSAman22373.t1;Note=
#chr_065 Manual  mRNA    1738186 1741894 .       -       .       ID=chr_065-GSAman22373.t1;mRNA_id=chr_065-GSAman22373
#chr_065 Manual  three_prime_UTR 1738186 1738295 .       -       .       ID=chr_065-GSAman22373.t1.utr;Parent=chr_065-GSAman22373.t1
#chr_065 Manual  CDS     1738296 1741793 .       -       .       ID=chr_065-GSAman22373.t1.cds;Parent=chr_065-GSAman22373.t1
#chr_065 Manual  five_prime_UTR  1741794 1741894 .       -       .       ID=chr_065-GSAman22373.t1.utr;Parent=chr_065-GSAman22373.t1
open(RAWGFF,">","$dir/Rename_gene/TGD2024_mRNA_rmdup.gtf");
#chr_025 GSAman  mRNA    47137   49020   .       -       .       ID=GSAman01311291gene;ManInfo=LQ;
#chr_025 GSAman  exon    47137   47257   .       -       .       Parent=GSAman01311291gene;
foreach(<IN8>){
    if(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?);.+?=(.+?)\s/){
        $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5;    $ori2024{$name}=$ori;
        $line_gene="$chr\tGSAman\tgene\t$start\t$end\t.\t$ori\t.\tID=$name;Name=$name.t1;Note=\n";
        $line="$chr\tGSAman\tmRNA\t$start\t$end\t.\t$ori\t.\tID=$name.t1;Parent=$name\n";
        #print "$line";
        $mid_point_2024{$chr}{$name}=int(($start+$end)/2); push@{$gff2024{$chr}{$name}},$line_gene,$line;
        if($name=~/GSAman(\d+)/ && length($1)==8){ $New_Tname{$name}="TTHERM_$1"; $Tname="TTHERM_$1"; $Used_Tname{$Tname}++; $Named{$name}=$Tname; $LQ{$name}=$Tname; }else{print "$name\n";}
        print RAWGFF "$line";
    }elsif(/(chr_\d+)\s.+?exon\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?Parent=(.+?);\s/){
        $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5; $line="$chr\tGSAman\tCDS\t$start\t$end\t.\t$ori\t.\tID=$name.t1.cds;Parent=$name.t1\n";
         push@{$gff2024{$chr}{$name}},$line; push@{$CDS2024{$name}},$start,$end;
    }elsif(/(chr_\d+)\s.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?Parent=(.+?);\s/){
         $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5; $line="$chr\tGSAman\tCDS\t$start\t$end\t.\t$ori\t.\tID=$name.t1.cds;Parent=$name.t1\n";
         push@{$gff2024{$chr}{$name}},$line; push@{$CDS2024{$name}},$start,$end;
    }else{print "$_";}
}
$n=keys%Used_Tname;
print "After Low_Quality_gene_LQ.gff3. The Used name:$n\n";
foreach(<IN7>){
    if(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1;.+?=(.+?)\s/){
        $chr=$1; $start=$2; $end=$3; $ori=$4; $name="$chr-$5";    $ori2024{$name}=$ori;
        $mid_point_2024{$chr}{$name}=int(($start+$end)/2);
        if($name=~/GSAman(\d+)/ && length($1)==8){ $New_Tname{$name}="TTHERM_$1"; $Tname="TTHERM_$1"; $Used_Tname{$Tname}++; $Named{$name}=$Tname; }
        $line="$chr\tGSAman\tmRNA\t$start\t$end\t.\t$ori\t.\tID=$name.t1;Parent=$name\n";
        print RAWGFF "$line"; push@{$gff2024{$chr}{$name}},$line;
    }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?);Name=(.+?.t1);Note=/){ $name="$1-$5";
        $line="$1\tGSAman\tgene\t$2\t$3\t.\t$4\t.\tID=$name;Name=$name.t1;Note=\n";
         push@{$gff2024{$1}{$name}},$line;
    }elsif(/(chr_\d+)\s.+?five_prime_UTR\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1.utr;Parent=(.+?).t1\s/){$name="$1-$5";
         $line="$1\tGSAman\tfive_prime_UTR\t$2\t$3\t.\t$4\t.\tID=$name.t1.utr;Parent=$name.t1\n"; push@{$gff2024{$1}{$name}},$line;
    }elsif(/(chr_\d+)\s.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1.cds;Parent=(.+?).t1\s/){ $name="$1-$5";
         $line="$1\tGSAman\tCDS\t$2\t$3\t.\t$4\t.\tID=$name.t1.cds;Parent=$name.t1\n";
          push@{$gff2024{$1}{$name}},$line; push@{$CDS2024{$name}},$2,$3;
    }elsif(/(chr_\d+)\s.+?three_prime_UTR\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1.utr;Parent=(.+?).t1\s/){ $name="$1-$5";
         $line="$1\tGSAman\tthree_prime_UTR\t$2\t$3\t.\t$4\t.\tID=$name.t1.utr;Parent=$name.t1\n";  push@{$gff2024{$1}{$name}},$line;
    }else{print "$_";}
}
$n=keys%Used_Tname;
print "After Low_Quality_gene_GENE.fasta_orf1_ATG2TGA_longest_UTR.gff3.\nUsed name:$n\n";
@IN1=<IN1>;close IN1;
foreach(@IN1){
    if(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1;.+?=(.+?)\s/){
        next if exists$NeedDeleteGene{$5};
        print RAWGFF "$_";
        $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5;    $ori2024{$name}=$ori;
        $mid_point_2024{$chr}{$name}=int(($start+$end)/2); push@{$gff2024{$chr}{$name}},$_;
        if($name=~/(chr_\d+)-GSAman(\d+)/ && length($2)==8){ $New_Tname{$name}="TTHERM_$2"; $Tname="TTHERM_$2"; $Used_Tname{$Tname}++; $Named{$name}=$Tname; }
    }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s.+?ID=(.+?);Name=(.+?.t1);Note=/){
        next if exists$NeedDeleteGene{$4};
        $name=$4; push@{$gff2024{$1}{$name}},$_;
    }elsif(/(chr_\d+)\s.+?five_prime_UTR\s+(\d+)\s+(\d+)\s.+?ID=(.+?).t1.utr;Parent=(.+?).t1\s/){
        next if exists$NeedDeleteGene{$4}; $name=$4; push@{$gff2024{$1}{$name}},$_;
    }elsif(/(chr_\d+)\s.+?CDS\s+(\d+)\s+(\d+)\s.+?ID=(.+?).t1.cds;Parent=(.+?).t1\s/){
        next if exists$NeedDeleteGene{$5}; $name=$4; push@{$gff2024{$1}{$name}},$_; push@{$CDS2024{$name}},$2,$3;
    }elsif(/(chr_\d+)\s.+?three_prime_UTR\s+(\d+)\s+(\d+)\s.+?ID=(.+?).t1.utr;Parent=(.+?).t1\s/){
        next if exists$NeedDeleteGene{$5}; $name=$4; push@{$gff2024{$1}{$name}},$_;
    }else{print "$_";}
}
close RAWGFF;
$cmd=`sort -k1.5,1.7n -k4.1,4.10n Rename_gene/TGD2024_mRNA_rmdup.gtf | awk '{match(\$9,/ID=(.+?).t1/,a); print a[1];}' > Rename_gene/TGD2024_need_rename_gene_sorted.txt`;
print $cmd;
#IN2
#chr_001 AUGUSTUS        gene    226     854     .       +       .       ID=g18410;Name=TTHERM_00161861;Note="hypothetical protein"
#chr_001 AUGUSTUS        mRNA    226     854     .       +       .       ID=g18410.t1;Parent=g18410
#chr_001 AUGUSTUS        CDS     226     342     .       +       0       ID=g18410.t1.cds;Parent=g18410.t1
foreach(<IN2>){
    if(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?);Name=(T.+?);Note="(.+)"\s/){
         $chr=$1; $start=$2; $end=$3; $ori=$4; $gname=$5; $Tname=$6;   $ori2021{$Tname}=$ori; $Tname2g{$Tname}=$gname;
         $mid_point_2021{$chr}{$Tname}=int(($start+$end)/2);      $info{$Tname}=$7;    $Tname2chr{$Tname}=$chr; $Used_gname{$gname}++;
         push@{$gff2021{$chr}{$Tname}},$_; $len_Tname{$Tname}=$end-$start+1; 
         $Tname=~/TTHERM_(\d+)/; if($1==9){$Unregular_Tname{$Tname}++; }
         $mark=1;
         foreach $name(keys%{$mid_point_2024{$chr}}){$mid=$mid_point_2024{$chr}{$name}; if($mid>=$start && $mid<=$end){  push@{$mid_in_2024{$Tname}},$name; } }
    }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?);Name=(g.+?);Note="(.+)"\s/){
          push@{$gff2021{$chr}{$Tname}},$_; $mark=0;
    }elsif(/(chr_\d+)\s.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1.cds;Parent=g.+?\s/ && $mark==1){
          push@{$gff2021{$chr}{$Tname}},$_;
    }elsif(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?)t1;/ && $mark==1){
          push@{$gff2021{$chr}{$Tname}},$_;
    }
}

foreach(@IN1){
   if(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1;.+?=(.+?)\s/){
      next if exists$NeedDeleteGene{$5};
      $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5;
      foreach $Tname(keys%{$mid_point_2021{$chr}}){$mid=$mid_point_2021{$chr}{$Tname}; if($mid>=$start && $mid<=$end){  push@{$mid_in_2021{$name}},$Tname; } }     
   }
}

##IN3
#chr_065-GSAman01355278gene   TTHERM_00151130 35.714  70      40      2       97      165     622     687     2.9     28.9
foreach(<IN3>){
   next if $_=~/^#/;
   @temp=split(/\s/,$_);
   next if $temp[1]=~/g\d+/;
   next if exists$NeedDeleteGene{$temp[0]};
   $temp[0]=~/((chr_\d+)-.+)/; $name=$1; $name_chr=$2;
   $Tname_chr=$Tname2chr{$temp[1]};
   push@{$blast_2021T2024{$name}},$temp[1] if ($name_chr eq $Tname_chr) or !exists$ori2021{$temp[1]};  
}
foreach $name(sort keys%blast_2021T2024){
     @temp=@{$blast_2021T2024{$name}};
     @uniq=grep{++$ha{$_}<2}@temp; #rmdup
     @{$blast_2021T2024_uniq{$name}}=@uniq;   
}

if(exists$Used_Tname{"TTHERM_00293350"}){ print "AAA!"; }
close IN1; close IN2; close IN3; close IN4;
open(OUT1,">","$dir/Rename_gene/TGD2024-2021_uniq_name_map.txt");
open(OUT2,">","$dir/Rename_gene/TGD2024-2021_ori-changed_map.txt");
open(OUT3,">","$dir/Rename_gene/TGD2024_new-gene_need_rename.txt");

foreach $name(keys%New_Tname){
    $Tname=$New_Tname{$name};
    print OUT1 "$name\t$Tname\t$Tname2g{$Tname}\tManual_Named\n";
    $Named{$name}=$New_Tname{$name};   $Used_Tname{$Tname}++;
}
$n=keys%Named;
print "Named gene:$n\n";

foreach $name(sort keys%mid_in_2021){
   @temp_2021=@{$mid_in_2021{$name}};
   if(exists$Named{$name}){ print OUT1 "$name\t$Named{$name}\t$Tname2g{$Named{$name}}\tManual_Named\n" ; $Tname=$Named{$name}; $Used_Tname{$Tname}++; next; }
   if(exists$New_Tname{$name}){ print OUT1 "$name\t$New_Tname{$name}\t$Tname2g{$New_Tname{$name}}\tManual_Named\n"; $Tname=$New_Tname{$name}; $Used_Tname{$Tname}++; next; }
   foreach $Tname(@temp_2021){ @temp_2024=@{$mid_in_2024{$Tname}}; 
      if(@temp_2021==1 && @temp_2024==1){
        $pos_uniq_match{$name}=$Tname; 
        $gname=$Tname2g{$Tname}; 
        next if exists$Used_Tname{$Tname};
        if($Tname=~/TTHERM_(\d+)/ && length($1)==8){print OUT1 "$name\t$Tname\t$gname\tone2one\n"; $Named{$name}=$Tname; $Used_Tname{$Tname}++;
        }elsif($Tname=~/TTHERM_(\d+)/ && length($1)==9){ print OUT1 "$name\t$Tname\t$gname\told-unregular-name\n"; $Named{$name}=$Tname; $Used_Tname{$Tname}++;  }
        if($ori2021{$Tname} ne $ori2024{$name}){ print OUT2 "$name:$ori2024{$name}\t$Tname:$ori2021{$Tname}\n";  }
     }elsif(@temp_2021==1 && @temp_2024>1){
        #next if $temp_2021[0]=~/TTHERM_(\d+)/ && length($1)!=8;
        foreach $line(@temp_2024){  
            next if exists$Used_Tname{$line};
            print OUT1"$name\t$temp_2021[0]\t$Tname2g{$temp_2021[0]}\t2021-one:2024-more\n";
            $Named{$name}=$temp_2021[0]; $Used_Tname{$Tname}++;
        }
     }elsif(@temp_2021>1 && @temp_2024==1){
        $line=join("\t",@temp_2021); #print "$temp_2024[0]\t$line\n";
        undef$max;
        foreach $line(@temp_2021){
            next if exists$Used_Tname{$line};
            #if($line=~/TTHERM_\d+0$/){ shift@max; push@max,$line; last; }
            if($len_Tname{$line}>$max){$max=$len_Tname{$line}; shift@max; push@max,$line; }
        }
        if(@max==1){ print OUT1 "$name\t$max[0]\t$Tname2g{$max[0]}\t2021-more:2024-one\n"; 
             $Named{$name}=$max[0]; $Used_Tname{$max[0]}++; }
     }elsif(@temp_2021==1 or @temp_2024==1){
           $line1=join("\t",@temp_2021); $line2=join("\t",@temp_2024);
     }else{
          $line1=join("\t",@temp_2021); $line2=join("\t",@temp_2024);
     }  
   }
}

$n=keys%ori2024; print "gene number:$n\n";
$i=0;
foreach $name(sort keys%ori2024){
     #$i++;
     next if exists$Named{$name};
     @temp=@{$blast_2021T2024_uniq{$name}};
     if(@temp==0){ print OUT3 "$name\n"; $unnamed_gene{$name}++;  next;  }
     undef$new_Tname; undef$gname;
     foreach $Tname(@temp){
         if(!exists$Used_Tname{$Tname}){  $new_Tname=$Tname; 
         do{srand($i); $num=int(rand(30000)); $gname="g$num"; $i++;}while(exists$Used_gname{$gname}); $Used_gname{$gname}++;  last;   }
     }
      #print OUT1 "$name\t$new_Tname\t$gname\tBlast\n"; $Named{$name}=$new_Tname;
     if($new_Tname=~/TTHERM/){print OUT1 "$name\t$new_Tname\t$gname\tBlast\n"; $Named{$name}=$new_Tname; $Used_Tname{$new_Tname}++; }
}
 close OUT2; close OUT3;
open(IN,"<","$dir/Rename_gene/TGD2024_need_rename_gene_sorted.txt") or die;
print "-----\n";
foreach(<IN>){  if(/^((chr_\d+)-.+?)\s/){ push@{$name_sorted{$2}},$1; }   }
$j=0;
foreach$chr(sort{$a cmp $b}keys%name_sorted){
   $i=0; @temp=@{$name_sorted{$chr}}; 
   while($i<=$#temp){  
      $name=$temp[$i];
       if(exists$Named{$name} && $Named{$name}=~/(\d+)/ && length($1)==8){ $i++; $j++; $new_Tname=$Named{$name}; next; 
       }elsif(exists$Named{$name} && $Named{$name}=~/(0(\d+))/ && length($2)==8 && !exists$Used_Tname{"TTHERM_$2"}){   $new_Tname="TTHERM_$2"; $Named{$name}=$new_Tname; $Used_Tname{$new_Tname}++; }
      $mark_before=0; $mark_after=0;
      $before=$i-1; $n=0; if($before<0){ $mark_before=-1; }else{ until($Named{$temp[$before]}=~/TTHERM_(\d+)/ && length($1)==8){ $before--; $n++;  if($before<0){$mark_before=-1; last;} } }
      $after=$i+1;  $m=0; if($after>$#temp){ $mark_after=-1; }else{ until($Named{$temp[$after]}=~/TTHERM_(\d+)/ && length($1)==8){ $after++; $m++; if($after>$#temp){$mark_after=-1; last;} } }
      $num=$n+$m+1; $Named{$temp[$before]}=~/TTHERM_(\d+)/; $id_before=$1; $Named{$temp[$after]}=~/TTHERM_(\d+)/; $id_after=$1;
      if($mark_before==-1){ until(!exists$Used_Tname{$new_Tname}){ $id_after=$id_after-2*$num;  $new_Tname="TTHERM_".sprintf("%08s",$id_after);$mark_before=0; } 
      }elsif($mark_after==-1){ until(!exists$Used_Tname{$new_Tname}){ $id_before=$id_before+2*$num;    $new_Tname="TTHERM_".sprintf("%08s",$id_before); $mark_after=0; } 
      }else{ $bin=int(($id_after-$id_before)/($num+1)); if($bin==0){$bin=2;} $b=$id_before; until(!exists$Used_Tname{$new_Tname}){ $b=$b+$bin;  $new_Tname="TTHERM_".sprintf("%08s",$b);   }    }
      do{srand($j); $num=int(rand(30000)); $gname="g$num"; $j++;}while(exists$Used_gname{$gname}) ; $Used_gname{$gname}++;   
      print  OUT1 "$name\t$new_Tname\t$gname\tCalculated\n"; $Used_Tname{$new_Tname}++; $Named{$name}=$new_Tname;
      $i++; $j++;
   } 
}
close OUT1;

$n=keys%Named;
print "Named gene:$n\nj num:$j\n";
open(LQ,"<","$dir/Rename_gene/LowQuality.txt");
foreach(<LQ>){ chomp; $LowQuality{$_}++;  } close LQ;
open(MAP,"<","$dir/Rename_gene/TGD2024-2021_uniq_name_map.txt");
open(OUT1,">","$dir/Rename_gene/TGD2024_renamed.gff3");
open(LowGene,">","$dir/Rename_gene/LowQuality_Tname2gname.txt");
open(T2G,">","$dir/Rename_gene/TGD2024_Tname2gname.txt");
undef%Used_gname;
foreach(<MAP>){
   #print "$_";
   @temp=split(/\t/);
   $name=$temp[0]; $Tname=$temp[1]; $gname=$temp[2];
    #print "error:$name\tUsed_Tname:$Used_Tname{$name}\tTname:$Tname\n" if $Used_Tname{$name} ne $Tname;
   if($gname=~/g\d+/){ $Used_gname{$gname}++; $MAP_Tname2gname{$Tname}=$gname;  } 
}
$n=keys%MAP_Tname2gname;
print "gene has gname:$n\n";
open(OneLineAnno,">","$dir/TGD2024_Tname2Annotation.xlxs");
#if(exists $gff2024{"chr_001"}{"chr_001-GSAman00161850gene"}){print "existst!\n";}else{print "error\n";}
#if(exists $Named{"chr_001-GSAman00161850gene"}){print "existst!\n";}else{print "error\n";}
foreach $chr(sort{$a cmp $b} keys%name_sorted){
   $i=0; @gene=@{$name_sorted{$chr}};
   while($i<=$#gene){   
      $name=$gene[$i];
      next if !exists$Named{$name};
      @temp=@{$gff2024{$chr}{$name}};
      foreach $line(@temp){
         if($line=~/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1;.+?=(.+?)\s/){
             $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5; $Tname=$Named{$name}; $gname=$MAP_Tname2gname{$Tname};
             print OUT1 "$chr\tManual\tmRNA\t$start\t$end\t.\t$ori\t.\tID=$gname.t1;Parent=$gname;\n";
             #print OUT1 "$chr\tManual\tmRNA\t$start\t$end\t.\t$ori\t.\tID=$Tname.t1;Parent=$Tname;\n";
         }elsif($line=~/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?);Name=(.+?.t1);Note=/){
             $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5; $Tname=$Named{$name}; 
             if(exists$Username{$Tname}){$UserName=$Username{$Tname};}else{$UserName="Unnamed";} if(exists$Anno{$Tname}){$Anno=$Anno{$Tname};}else{$Anno="Unknow";}
             if(!exists$MAP_Tname2gname{$Tname}){
                do{srand($j); $num=int(rand(30000)); $gname="g$num"; $j++;}while(exists$Used_gname{$gname}) ; $Used_gname{$gname}++;
                 $MAP_Tname2gname{$Tname}=$gname;
             }else{ $gname=$MAP_Tname2gname{$Tname}; }
             print OUT1 "$chr\tManual\tgene\t$start\t$end\t.\t$ori\t.\tID=$gname;Name=$Tname;UserNamed=$UserName;TGDAnno=$Anno;info=$info{$Tname};\n" if !exists$InterPro_Anno{$name};
             print OUT1 "$chr\tManual\tgene\t$start\t$end\t.\t$ori\t.\tID=$gname;Name=$Tname;UserNamed=$UserName;TGDAnno=$Anno;info=$info{$Tname};Annotation=$InterPro_Anno{$name};\n" if exists$InterPro_Anno{$name};
             print OneLineAnno "$Tname\t$UserName\t$gname\tTGDAnno=$Anno\tInterProAnno=$InterPro_Anno{$name}\n";
             #print OUT1 "$chr\tManual\tgene\t$start\t$end\t.\t$ori\t.\tID=${Tname};Name=$Tname;UserNamed=$UserName;TGDAnno=$Anno;info=$info{$Tname};Annotation=$InterPro_Anno{$name};\n";
             print T2G "$chr\t$Tname\t$gname\t$name\n";
         }elsif($line=~/(chr_\d+)\s.+?five_prime_UTR\s+(\d+)\s+(\d+)\s+.\s+(.)\s+.+?ID=(.+?).t1.utr;/){
             $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5; $Tname=$Named{$name}; $gname=$MAP_Tname2gname{$Tname};
             print OUT1 "$chr\tManual\tfive_prime_UTR\t$start\t$end\t.\t$ori\t.\tID=$gname.t1.utr;Parent=$gname.t1;\n";
             #print OUT1 "$chr\tManual\tfive_prime_UTR\t$start\t$end\t.\t$ori\t.\tID=$Tname.t1.utr;Parent=$Tname;\n";
         }elsif($line=~/(chr_\d+)\s.+?CDS\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1.cds;/){
             $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5; $Tname=$Named{$name}; $gname=$MAP_Tname2gname{$Tname};
             print OUT1 "$chr\tManual\tCDS\t$start\t$end\t.\t$ori\t.\tID=$gname.t1.cds;Parent=$gname.t1;\n";
             #print OUT1 "$chr\tManual\tCDS\t$start\t$end\t.\t$ori\t.\tID=$Tname.t1.cds;Parent=$Tname;\n";
         }elsif($line=~/(chr_\d+)\s.+?three_prime_UTR\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1.utr;/){
             $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5; $Tname=$Named{$name}; $gname=$MAP_Tname2gname{$Tname};
             print OUT1 "$chr\tManual\tthree_prime_UTR\t$start\t$end\t.\t$ori\t.\tID=$gname.t1.utr;Parent=$gname.t1;\n";
             #print OUT1 "$chr\tManual\tthree_prime_UTR\t$start\t$end\t.\t$ori\t.\tID=$Tname.t1.utr;Parent=$Tname;\n";
         }else{print "$_";}
      }
      $i++;
   }
   foreach $Tname(keys%{$gff2021{$chr}}){
     if(exists$LowQuality{$Tname} && !exists$MAP_Tname2gname{$Tname}){ 
        @temp=@{$gff2021{$chr}{$Tname}};
        #foreach $line(@temp){ print OUT1"$line"; }
        print LowGene "$Tname\t$Tname2g{$Tname}\tLowQuality\n";
     }
   }
}



















