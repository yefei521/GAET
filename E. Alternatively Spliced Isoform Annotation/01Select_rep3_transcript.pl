use strict;
my(%gff, $name);
open(IN1,"<$ARGV[0]") or die;
open(OUT1,">$ARGV[0]_common") or die;

foreach(<IN1>){
   if(/transcript_id "(.+?)";.+? class_code "j";.+? num_samples "2";/){  $name=$1; $gff{$name}=$_;
   }elsif(/transcript_id "($name)";/){$gff{$name}.=$_;  }
}

foreach(keys%gff){
  print OUT1"$gff{$_}";
}
