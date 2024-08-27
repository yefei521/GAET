use strict;
my$dir=$ENV{'PWD'};
my($gene_name, %SB210, %con, %veg, %star, $m);
open(IN1,"<$dir/05EpiOri_ne-gene_Checked.txt") or die; my@EpiOri_ne_gene=<IN1>; close IN1;
open(IN3,"<$dir/SB210_DRS_DNA_Q7_maped_sorted.gtf_UTR.gff3_Filtered_sorted") or die; my@SB210_filter_geng=<IN3>; close IN3;
open(IN4,"<$dir/conjugation_assembly_accepted_hits_unique_sorted_rmdup.gtf_UTR.gff3_Filtered_sorted") or die; my@con_filter_geng=<IN4>; close IN4;
open(IN5,"<$dir/veg_assembly_accepted_hits_unique_sorted_rmdup.gtf_UTR.gff3_Filtered_sorted") or die; my@veg_filter_geng=<IN5>; close IN5;
open(IN6,"<$dir/starvation_assembly_accepted_hits_unique_sorted_rmdup.gtf_UTR.gff3_Filtered_sorted") or die; my@star_filter_geng=<IN6>; close IN6;

foreach my$m(@EpiOri_ne_gene){ if($m=~/^(SB210_DRS_DNA_Q7_maped_sorted.gtf_UTR.gff3_Filtered_sorted)\t(.+?)\n/){  $SB210{$2}=0;  }  }
foreach my$m(@EpiOri_ne_gene){ if($m=~/^(conjugation_assembly_accepted_hits_unique_sorted_rmdup.gtf_UTR.gff3_Filtered_sorted)\t(.+?)\n/){  $con{$2}=0;  }   }
foreach my$m(@EpiOri_ne_gene){ if($m=~/^(veg_assembly_accepted_hits_unique_sorted_rmdup.gtf_UTR.gff3_Filtered_sorted)\t(.+?)\n/){  $veg{$2}=0;  }   }
foreach my$m(@EpiOri_ne_gene){ if($m=~/^(starvation_assembly_accepted_hits_unique_sorted_rmdup.gtf_UTR.gff3_Filtered_sorted)\t(.+?)\n/){ $star{$2}=0;  }   }

foreach my$line2(@SB210_filter_geng){
       if($line2=~/gene\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+);Name/){ if(exists$SB210{$3}){ $gene_name=$3; $SB210{$3}=$line2; $m=1; }else{ $m=0; }
       }elsif($line2=~/five_prime_UTR\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){  $SB210{$3}.=$line2;
       }elsif($line2=~/CDS\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){  $SB210{$3}.=$line2;
       }elsif( $line2=~/three_prime_UTR\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){ $SB210{$3}.=$line2;
    }
}
foreach(keys%SB210){print "$_\n$SB210{$_}"; }
foreach my$line2(@con_filter_geng){
       if($line2=~/gene\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+);Name/){ if(exists$con{$3}){ $gene_name=$3; $con{$3}.=$line2; $m=1; }else{ $m=0; }
       }elsif($line2=~/five_prime_UTR\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){  $con{$3}.=$line2;
       }elsif($line2=~/CDS\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){  $con{$3}.=$line2;
       }elsif( $line2=~/three_prime_UTR\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){ $con{$3}.=$line2;
    }
}
foreach(keys%con){print "$con{$_}"; }
foreach my$line2(@veg_filter_geng){
       if($line2=~/gene\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+);Name/){ if(exists$veg{$3}){ $gene_name=$3; $veg{$3}.=$line2; $m=1; }else{$m=0;}
       }elsif($line2=~/five_prime_UTR\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){ $veg{$3}.=$line2;
       }elsif($line2=~/CDS\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){  $veg{$3}.=$line2;
       }elsif( $line2=~/three_prime_UTR\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){ $veg{$3}.=$line2;
    }
}
foreach(keys%veg){print "$veg{$_}"; }
foreach my$line2(@star_filter_geng){
       if($line2=~/gene\t(\d+)\t(\d+)\t.+?ID=((.+?)_ORF\d+);Name/){ if(exists$star{$3}){ $gene_name=$3; $star{$3}.=$line2; $m=1; }else{$m=0; }
       }elsif($line2=~/five_prime_UTR\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){  $star{$3}.=$line2;
       }elsif($line2=~/CDS\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){  $star{$3}.=$line2;
       }elsif( $line2=~/three_prime_UTR\t(\d+)\t(\d+)\t.+?ID=($gene_name)/ && $m==1){  $star{$3}.=$line2;
    }
}
foreach(keys%star){print "$star{$_}"; }

