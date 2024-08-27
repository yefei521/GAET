#! /usr/bin/perl -w
##gffread Final.evm.update.gff3 -g ../run/ATCC.contig.fa -x ATCC.cds.fa
use strict;
die "#usage;perl $0 <input.cds.fa><out.pep.fa>\n" unless @ARGV==2;
my $incds=shift;
my $outpep=shift;
cds2pep($incds,$outpep);
sub cds2pep{
	my ($infile,$outfile)=@_;
	open IN,'<',$infile||die;
	open OUT,'>',$outfile||die;
	my $p=code();
	$/=">";<IN>;$/="\n";
	while(<IN>){
		chomp;
		my $head=$_;
		$/=">";
		chomp(my $seq=<IN>);
		$/="\n";
		$seq=~s/\n+//g;
		my $out;
		for(my $i=0;$i<length$seq;$i+=3){
			my $codon=uc(substr($seq,$i,3));
			last if (length$codon <3);
			$out.= exists $p->{"standard"}{$codon} ? $p->{"standard"}{$codon} : "X";
		}
		$out =~ s/U$//;
		my $len=length$out;
		$out =~ s/([A-Z]{50})/$1\n/g;
		chop $out unless $len % 50;
		print OUT ">$head [translate_table: 6]\n$out\n"
	}
	close OUT;
}

sub code{
	my $p={
        	"standard" =>
                 	{  
               'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',
               'TGC' => 'C', 'TGT' => 'C', 
               'GAC' => 'D', 'GAT' => 'D', 
               'GAA' => 'E', 'GAG' => 'E',
               'TTC' => 'F', 'TTT' => 'F',
               'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',
               'CAC' => 'H', 'CAT' => 'H',
               'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I', 
               'AAA' => 'K', 'AAG' => 'K', 
               'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',
               'ATG' => 'M', 
               'AAC' => 'N', 'AAT' => 'N',
               'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',
               'CAA' => 'Q', 'CAG' => 'Q',
               'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',
               'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',
               'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',
               'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',
               'TGG' => 'W',
               'TAC' => 'Y', 'TAT' => 'Y',
               'TAA' => 'Q', 'TAG' => 'Q', 'TGA' => '*'
          }
	};
	return $p;
}
__END__
