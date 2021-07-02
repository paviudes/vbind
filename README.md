# vbind

vbind is a userfriedly tool for RNA sequencing. In particular, it is used to compute and visualize the bindings between a pool of nuleotides and a gene sequence.

## Usage

In what follows, we will describe how to use our software tool to compute matchings between a gene and a nucleotide pool.

1. Key ingriedients: The __gene__ and the __pool__.
  * __gene__: The gene sequence should be a string of characters from the alphabet {`A`, `T`, `G`, `C`}. We require the gene sequence to be specified in a text file without any line-breaks.
  * __pool__: The pool is a collection of several sRNA nucleotides, each represented as a string of characters from the alphabet {`A`, `T`, `G`, `C`}. Ideally, we require a text file, wherein each new line identifies a nucleotide sequence. However, since a pool is typically derived from a gene bank, we accomodate some other formats for specifying the pool. In particular, we accept a text file containing the nueotides of the pool, interleaved by other information that is irrelevant for sequencing. Hence, the pool should be specified by a text file, wherein every `k`-th line, for some `k > 0`.
2. Other settings:
  * __lines to skip__: number of lines to skip (denoted above by `k`) in the text file describing the pool, before reading a valid nuleoide sequence.
  * __tolerance__: the maximum number of mismatches allowed
  * __topology of the gene__: an integer that takes the value 0 for linear matching and 1 for circular matching
  * __cores__: the number of cores to be used by the software.
 
For example, “gene.txt pool.txt 4 1 1 1” is a complete input specification, indicating that the problem of computing bindings, in the circular topology while allowing for at most one mismatch, between the sequence in gene.txt with those in the pool.txt. Furthermore, every fourth sequence in pool.txt is a valid sRNA nuleotide. All such instances of the matching problem can be gathered in a text file, placed in vbind/data/input.

To solve the instances of the matching problem in a text file “example.txt”, we need to run the following command:
`./vbind.sh example.txt`
and to solve only a particular instance, we can specify its line number: `x`, by
`./vbind.sh example.txt x` .

## Installation

vbind is entriely written in Python. The following packages are required for its seamless functioning.

| Package         | Description                                                                  |
|-----------------|------------------------------------------------------------------------------|
| numpy           | [Linear algebra](https://www.numpy.org)                                      |
| scipy           | [Linear algebra](https://www.scipy.org)                                      |
| tqdm            | [Progress bar](https://tqdm.github.io)                                       |
| multiprocessing | [Parallel execution](https://docs.python.org/3/library/multiprocessing.html) |
| datetime        | [Date and Time](https://docs.python.org/3/library/datetime.html)             |


## Contact

Please contact the following developers for any additional help or requests for customization.
* Pavithran Iyer: pavithran.iyer@uwaterloo.ca .
* Charith Adkar: Charith.Raj.Adkar.Purushothama@USherbrooke.ca .
