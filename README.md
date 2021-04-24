# vbind
In what follows, we will describe how to use our software tool to compute matchings between a gene and a nucleotide pool. We require the gene sequence to be specified in a text file without any line-breaks. The pool should also be specified in a text file wherein each new line identifies a nucleotide sequence composed of “A”, “T”, “G” and “C” characters. Note that nucleotide sequences with “N” will be ignored. Once the gene and the pool are specified in their respective text files, eg. gene.txt and pool.txt, we can provide additional settings for the matching problem. These include (i) tolerance: the maximum number of mismatches allowed, (ii) topology of the gene, i.e., an integer that takes the value 0 for linear matching and 1 for circular matching, and (iii) the number of cores to be used by the software. For example,
“gene.txt pool.txt 1 1 1” is a complete input specification, indicating that the problem of computing bindings, in the circular topology while allowing for at most one mismatch, between the sequence in pool.txt with those in the pool.txt. All such instances of the matching problem can be gathered in a text file, placed in vbind/data/input.

To solve the instances of the matching problem in a text file “inputs.txt”, we need to run the following command:
./vbind.sh inputs.txt
and to solve only a particular instance, we can specify its line number: x, by
./vbind.sh inputs.txt x .
