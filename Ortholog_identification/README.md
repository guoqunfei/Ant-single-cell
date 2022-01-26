**The pipeline of ortholog identification**
1)  Obtain the genome files and protein sequence files of the two species aligned.
2)  Obtain an annotation file whose ID is protein ID and only retains CDS information.
3)  Extract the start codon and stop codon sequences of mRNA according to the annotation file, and judge the correctness of the open reading frame; if the information cannot be extracted from some mRNA, it needs to be filtered.
4)  According the correctness of the open reading frame or the length of the protein, the optimal mRNA in a locus is screened.
5)  Generate gene annotation files based on the results of the screening, including GFF， CDS and PEP.
6)  Merge PEP file of all species and use BLAST for library building alignment.
7)  According to the comparison results, obtain the results of the pairwise comparison.
