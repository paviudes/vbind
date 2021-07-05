## Welcome to GitHub Pages

![Logo](https://github.com/paviudes/vbind/blob/master/logo.png?raw=true)

__vbind__ is a userfriedly tool for RNA sequencing. In particular, it can be used to compute and visualize the bindings between a pool of sRNA nuleotides and a genome sequence.

## Usage

In what follows, we will describe the two parts of our software tool. The first part computes matchings between a gene and a nucleotide pool. Second, presents some tools to visualize the solution to the mapping problem.

### Part one: Computing matchings between a gene and a nuleotide pool
1. Key ingriedients: The __gene__ and the __pool__.
  * __gene__: The gene sequence should be a string of characters from the alphabet {`A`, `T`, `G`, `C`}. We require the gene sequence to be specified in a text file without any line-breaks.
  * __pool__: The pool is a collection of several sRNA nucleotides, each represented as a string of characters from the alphabet {`A`, `T`, `G`, `C`}. Ideally, we require a text file, wherein each new line identifies a nucleotide sequence. However, since a pool is typically derived from a gene bank, we accomodate some other formats for specifying the pool. In particular, we accept a text file containing the nueotides of the pool, interleaved by other information that is irrelevant for sequencing. Hence, the pool should be specified by a text file, wherein every `k`-th line, for some `k > 0`.
2. Other settings:
  * __lines to skip__: number of lines to skip (denoted above by `k`) in the text file describing the pool, before reading a valid nuleoide sequence.
  * __tolerance__: the maximum number of mismatches allowed
  * __topology__ of the gene: an integer that takes the value 0 for linear matching and 1 for circular matching

A problem instance is specified as: `<gene> <pool> <lines to skip> <tolerance> <topology> <cores>`.

   For example, `gene.txt pool.txt 4 1 1 1` is a complete input specification, indicating that the problem of computing bindings, in the circular topology while allowing for at most one mismatch, between the sequence in `gene.txt` with those in the `pool.txt`. Furthermore, every fourth sequence in `pool.txt` is a valid sRNA nuleotide.

All such instances of the matching problem can be gathered in a text file, placed in vbind/data/input.

To solve the instances of the matching problem in a text file “example.txt”, we need to run the following command:
`./vbind.sh example.txt`
and to solve only a particular instance, we can specify its line number: `x`, by
`./vbind.sh example.txt x` .

### Part two: Postprocessing -- Visualize matching results

Matching results computed from part one can be visualized using plots in our software tool. In addition to plotting, we also offer the capability of normalizing to reads per million. The post processing steps can be executed from `./analyze.sh`.

## Installation

vbind is entriely developed and meant to be used with the _Python 3 interpretter_. The following packages are required for its seamless functioning.

| Package         | Description                                                                  |
|-----------------|------------------------------------------------------------------------------|
| numpy           | [Linear algebra](https://www.numpy.org)                                      |
| scipy           | [Linear algebra](https://www.scipy.org)                                      |
| tqdm            | [Progress bar](https://tqdm.github.io)                                       |
| multiprocessing | [Parallel execution](https://docs.python.org/3/library/multiprocessing.html) |
| datetime        | [Date and Time](https://docs.python.org/3/library/datetime.html)             |

## License

This software is licensed under the [BSD Clause 3](https://opensource.org/licenses/BSD-3-Clause) license.

## Support or Contact

Having trouble with vbind? Check out our [documentation](https://github.com/paviudes/vbind) or [contact us](mailto:pavithran.iyer@uwaterloo.ca,charith.adkar@usherbrooke.ca?subject=[vbind%20querry]) and we’ll help you sort it out.
