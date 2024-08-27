#!/usr/bin/perl

use strict;
use warnings;
##gffread Final.evm.update.gff3 -g ../run/ATCC.contig.fa -x ATCC.cds.fa
#use strict;
die "#usage;perl $0 <input.exon.fa>\n" unless @ARGV==1;

# 定义新的密码子表
my %codon_table = (
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
);

# 读取 DNA 序列文件
#my $filename = shift or die "Usage: $0 <filename> \n";
open(my $fh, '<', $ARGV[0]) or die "Can't open $ARGV[0]: $!";

my $cds_output_file = "$ARGV[0]_cds.fasta";
my $protein_output_file = "$ARGV[0]_pep.fasta";

open(my $cds_fh, '>', $cds_output_file) or die "Can't open $cds_output_file for writing: $!";
open(my $protein_fh, '>', $protein_output_file) or die "Can't open $protein_output_file for writing: $!";

my %longest_cds;
my %longest_protein;

# 读取每个序列并预测 ORF
my$seq_id; my%input_seq;
while(<$fh>){
  if(/^>(.+)\s+/){
     $seq_id=$1;
   }else{
    chomp; $input_seq{$seq_id}.=$_;
   }
}
foreach $seq_id(keys%input_seq) {
        my $sequence = $input_seq{$seq_id};
        chomp $sequence;
        # 预测 ORF
        my ($cds, $protein, $start_pos, $end_pos) = predict_orf($sequence);

        # 更新最长的 CDS 和蛋白序列
        if (!exists $longest_cds{$seq_id} || length($cds) > length($longest_cds{$seq_id}{cds})) {
            $longest_cds{$seq_id}{cds} = $cds;
            $longest_cds{$seq_id}{length} = length($cds);
            $longest_cds{$seq_id}{protein} = $protein;
            $longest_cds{$seq_id}{start_pos} = $start_pos;
            $longest_cds{$seq_id}{end_pos} = $end_pos;
        }
}

close($fh);

# 输出最长的 CDS 和蛋白序列到文件
foreach my $seq_id (keys %longest_cds) {
    my $cds_length = length($longest_cds{$seq_id}{cds});
    my $start_pos = $longest_cds{$seq_id}{start_pos};
    my $end_pos = $longest_cds{$seq_id}{end_pos};
    print $cds_fh ">$seq_id:$start_pos-$end_pos:$cds_length\n";
    print $cds_fh "$longest_cds{$seq_id}{cds}\n";

    print $protein_fh ">$seq_id\n";
    print $protein_fh "$longest_cds{$seq_id}{protein}\n";
}

close($cds_fh);
close($protein_fh);

# 预测 ORF 的子程序
sub predict_orf {
    my ($sequence) = @_;

    my $longest_cds = '';
    my $longest_protein = '';
    my $start_pos = 0;
    my $end_pos = 0;

    # 寻找起始密码子 'ATG'
    while ($sequence =~ /ATG/g) {
        my $temp_start_pos = pos($sequence) - 2;
        my $temp_end_pos = 0;
        my $cds = '';
        my $protein = '';

        # 记录起始位置
        if ($start_pos == 0) {
            $start_pos = $temp_start_pos;
        }

        # 从起始密码子开始翻译
        for (my $i = $temp_start_pos - 1; $i < length($sequence) - 2; $i += 3) {
            my $codon = substr($sequence, $i, 3);
            if (exists $codon_table{$codon}) {
                my $amino_acid = $codon_table{$codon};
                $protein .= $amino_acid;
                $cds .= $codon;
                if ($amino_acid eq '*') {  # 如果遇到终止密码�止翻译
                    $temp_end_pos = $i + 3;
                    last;
                }
            } else {
                last;  # 如果遇到非法密码�止翻译
            }
        }

        # 记录结束位置
        if ($temp_end_pos > $end_pos) {
            $end_pos = $temp_end_pos;
        }

        # 更新最长的 CDS 和蛋白序列
        if (length($cds) > length($longest_cds)) {
            $longest_cds = $cds;
            $longest_protein = $protein;
        }
    }

    return ($longest_cds, $longest_protein, $start_pos, $end_pos);
}

