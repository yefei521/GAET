#conda activate DBSCAN
##sh multi-P_DBscan.sh > log 2>&1 &
#
$dir=$ENV{'PWD'};
@files=glob"$dir/gene_TES_cluster/*";
#open(IN1,"<","$dir/");

foreach $file(@files){
   open(IN1,"<","$file");
   foreach $line(<IN1>){
      if($line=~/.+?,"(\d+)",.+?\/(.+?t1).+?,"(\d+)"/ && $1>0){ #"1","1","gene_TES/YF00000040.t1","127390"
        push@{$TES{$2}{$1}},$3;
      }
   }
}

open(OUT1,">","$dir/Manual_check-total-gene.gff3_Right_UTR.gff3_exon-add-Low-Q-1271_TES.txt");
foreach $name(keys%TES){
    foreach $num(keys%{$TES{$name}}){
       @temp=@{$TES{$name}{$num}};
       print OUT1 "$name\t$num\t@temp\n";
    }
}
