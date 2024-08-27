#Introduction
LGC characterizes and identifies lncRNAs based on the relationship between ORF (open reading frame) Length and GC content. LGC is able to accurately distinguish lncRNAs from protein-coding RNAs in a cross-species manner without species-specific adjustments, and is robustly effective in discriminating lncRNAs from protein-coding RNAs across species that range from plants to mammals.

#Installation
Install Python 2.7 and Biopython:
Get Python 2.7 at https://www.python.org/download/ or install with your operating system’s package manager.
Get Biopython at http://biopython.org/wiki/Biopython .
You can follow the guide ( http://biopython.org/DIST/docs/install/Installation.html#sec4) to install Python and Biopython step by step.

$ tar zxf LGC-1.0.tar.gz # Depress LGC-1.0.tar.gz
$ python2.7 LGC-1.0.py input.fasta output.txt # Run LGC

Successful run of LGC will print as following:
$ Input: input.fasta # Input file
$ Output: output.txt # Output file
$ Scan ORF ... # Scan ORF and calculate coding potential score
$ Done # LGC runs to completion
$ Computation time XXX # Computation time of LGC

#Input
Fasta format.

#Output
There are nine columns in the output file.

Sequence name: name of transcript sequence
ORF Length: length of the longest ORF
GC Content: GC content of the longest ORF
Coding Potential Score: coding potential score for a transcript, which is protein-coding RNA if greater than 0 or ncRNA if smaller than 0. '0' indicates that mRNA probability equals lncRNA probability. Also, if the ORF length is shorter than 100nt, '0' is output.
Coding Label: "Coding" represents mRNA and “Non-coding” represents lncRNA.
pc: Probability of ORF for coding sequence
pnc: Probability of ORF for non-coding sequence
fc: Stop-codon probability for coding sequence
fnc:Stop-codon probability for non-coding sequence