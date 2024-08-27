open(IN1,"<","$ARGV[0]") or die;
my @string = <IN1>;
#my $key = "AATAAA";
#my$key="ATTAAA";
my$key="GTGT";
my(%fasta);
foreach(@string){
   chomp; 
   if(/>(.+)/){ $name=$1; 
   }else{ $fasta{$name}=$_; $num++;
      $pos = index( $fasta{$name}, $key, 1 );
      next if $pos<0;
      $pos+=2;
      $PAS{$pos}++;
      #print "$name\t$pos\n";
   }
}
#my(%PAS);
open(OUT1,">PAS-${key}_position.txt");


#foreach my$name(keys%fasta){ 
#  my ( $pos, $now ) = ( 0, -1 );
#  until ( $pos == -1 ) {
#     $pos = index( $fasta{$name}, $key, $now + 1 );
#     if($pos==$now){ $PAS{$name}=$fasta{$name}; }
#     $now = $pos ;
#  } 
#}

foreach(keys%PAS){
   $ratio=$PAS{$_}/$num;
   print OUT1"$_\t$PAS{$_}\t$ratio\n";
}
