$dir=$ENV{'PWD'};
open(IN1,"<","$dir/Rename_gene/TGD2024_duplicated_nedd_check.txt_checked") or die;
open(OUT1,">","$dir/Rename_gene/TGD2024_duplicated_nedd_check.txt_checked_cluster") or die;
foreach(<IN1>){
    if(/info:gene\s+query:(.+?)\s/){  $gene{$1}++; print OUT1"$1:gene\n";
    }elsif(/info:gene-(\d+)\s+query:(.+?)\s/){ $gene{$2}++; $gene_name{$2}="TTHERM_$1"; print OUT1 "$2:gene\tTTHERM_$1:rename\n";
    }elsif(/info:(AS.+?)\s+query:(.+?)\s/){ $AS{$2}=$1; print OUT1 "$2:$1\n";
    }elsif(/info:(NAT.+?)\s+query:(.+?)\s/){ $NAT{$2}=$1; print OUT1 "$2:$1\n";
    }elsif(/info:(.+?delete)\s+query:(.+?)\s/){ $delete{$2}++; print OUT1 "$2:$1\n";
    }elsif(/query:.+?\s+(.+?)\|(.+?)\s/){ $one2one{$1}=$2;
    }else{print "$_";}
}
close IN1;

open(IN2,"<","$dir/Rename_gene/TGD2024_need_rename_gene.gff3") or die;
@IN2=<IN2>;close IN2;
foreach(@IN2){
    if(/(chr_\d+)\s.+?mRNA\s+(\d+)\s+(\d+)\s+.\s+(.)\s.+?ID=(.+?).t1;.+?=(.+?)\s/){
        $chr=$1; $start=$2; $end=$3; $ori=$4; $name=$5;    $ori2024{$name}=$ori;
        $mid_point_2024{$chr}{$name}=int(($start+$end)/2); push@{$gff2024{$chr}{$name}},$_;
        $gene_len{$name}=$end-$start+1;
    }elsif(/(chr_\d+)\s.+?gene\s+(\d+)\s+(\d+)\s.+?ID=(.+?);Name=(.+?.t1);Note=/){
        $name=$4; push@{$gff2024{$1}{$name}},$_; 
    }elsif(/(chr_\d+)\s.+?five_prime_UTR\s+(\d+)\s+(\d+)\s.+?ID=(.+?).t1.utr;Parent=(.+?).t1\s/){
        $name=$4; push@{$gff2024{$1}{$name}},$_; $exon_num{$name}++;
    }elsif(/(chr_\d+)\s.+?CDS\s+(\d+)\s+(\d+)\s.+?ID=(.+?).t1.cds;Parent=(.+?).t1\s/){
        $name=$4; push@{$gff2024{$1}{$name}},$_; push@{$CDS2024{$name}},$2,$3; $CDS_len{$name}+=$3-$2+1; $CDS_num{$name}++; $exon_num{$name}++;
    }elsif(/(chr_\d+)\s.+?three_prime_UTR\s+(\d+)\s+(\d+)\s.+?ID=(.+?).t1.utr;Parent=(.+?).t1\s/){
        $name=$4; push@{$gff2024{$1}{$name}},$_; $exon_num{$name}++;
    }else{print "$_";}
}
foreach $name1(keys%one2one){
    $name2=$one2one{$name1};
    if($ori2024{$name1} ne $ori2024{$name2}){
        if($CDS_len{$name1}>$CDS_len{$name2}){$info{$name1}="gene";   $info{$name2}="NATnT3"; print OUT1 "$name1:gene\t$name2:NATnT3\n";
        }else{ $info{$name1}="NATnT3";   $info{$name2}="gene";  print OUT1 "$name2:gene\t$name1:NATnT3\n";  } 
    }else{
        if($CDS_len{$name1}>$CDS_len{$name2}){$info{$name1}="gene";   
           if($CDS_num{$name1}>$CDS_num{$name2}){$info{$name2}="ASnIR";}elsif($CDS_num{$name1}<$CDS_num{$name2}){ $info{$name2}="ASnES"; }else{$info{$name2}="delete";}
           print OUT1 "$name1:$info{$name1}\t$name2:$info{$name2}\n";
        }elsif($CDS_len{$name1}<$CDS_len{$name2}){ $info{$name2}="gene"; 
           if($CDS_num{$name2}>$CDS_num{$name1}){$info{$name1}="ASnIR";}elsif($CDS_num{$name2}<$CDS_num{$name1}){ $info{$name1}="ASnES"; }else{$info{$name1}="delete";}
           print OUT1 "$name1:$info{$name1}\t$name2:$info{$name2}\n";
        }else{ if($exon_num{$name2}>$exon_num{$name1}){ $info{$name1}="ASnES";   $info{$name2}="gene";  
               }elsif($exon_num{$name2}<$exon_num{$name1}){ $info{$name1}="gene";   $info{$name2}="ASnES";   
               }else{ if($gene_len{$name1}>$gene_len{$name2}){$info{$name1}="gene";   $info{$name2}="delete";
                      }else{ $info{$name2}="gene";   $info{$name1}="delete";   } 
               } 
               print OUT1 "$name1:$info{$name1}\t$name2:$info{$name2}\n"; 
             }
    }
}

