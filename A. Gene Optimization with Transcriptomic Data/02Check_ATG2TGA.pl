use strict;
my$dir=$ENV{'PWD'};
my(%Error, $name, %fasta, %mark, );
open(IN1,"<$dir/Thermophila_tetrahymena_2021_UTR_manualcheck_220224_cds.fasta"); #input CDS fasta
open(IN2,"<$dir/01Error_gene.txt"); #error gene name
open(OUT1,">$dir/02ATG2TGA_gene_result.txt");

while(<IN2>){ if(/.\t(.+?)\s+/){ $Error{$1}="Error gene";  } }
my$n=keys%Error; print "Error gene number\t$n\n";
while(<IN1>){
   if(/>(.+?)\n/){ $name=$1;
   }else{ chomp; $fasta{$name}.=$_;}
}
my$n=keys%fasta; print "Gene number\t$n\n";
$n=0;my$m=0;
foreach(keys%fasta){
    if($fasta{$_}=~/^ATG.+?TGA$/){ $mark{$_}="$_\tATG-TGA\tyes\t$Error{$_}\n"; $m++;
    }else{ $fasta{$_}=~/^(\w{3}).+?(\w{3})$/; $mark{$_}="$_\t$1-$2\tno\t$Error{$_}\n"; $n++; }
}

print "gene ATG2TGA number\t$m\n";
print "gene error number\t$n\n";
#output
foreach(keys%mark){
    print OUT1"$mark{$_}";
}
