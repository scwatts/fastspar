# FastSpar
FastSpar is a c++ implementation of the SparCC algorithm published here: Friedman, J. & Alm, E. J. Inferring correlation networks from genomic survey data. PLoS Comput. Biol. 8, e1002687 (2012). FastSpar is up to several thousand times faster than the original Python2 implementation and uses much less memory.

Many of the processes of FastSpar were ported from original SparCC implementation found here: https://bitbucket.org/yonatanf/sparcc

Additionally, SparCC's method of p-value has been replaced with exact p-value calculation. This code was adopted from the permp function from the statmod R package and from an example provided by Scott Ritchie. The advantages of this alternative method are discussed here: Phipson, B. & Smyth, G. K. Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn Stat Appl Genet Mol Biol. 9: Article 39 (2010).

FastSpar is currently in development and may lack certain features expected in complete programs. While FastSpar generally works without issue, extra care must be taken to ensure that correctly formatted input data are provided.


## Installing
FastSpar can be installed using a package containing pre-compiled binaries or installed from source. If running Ubuntu, it is recommended you install using the package.


### GNU/Linux
For most 64-bit linux distributions (e.g. Ubuntu, Debian, RedHat, etc) the easiest way to obtain FastSpar is via statically compiled binaries on the GitHub release page. These binaries can be downloaded and run immediately without any setup as they have no dependencies.


### macOS
To install FastSpar and required dependencies, homebrew can be used. The GitHub repository contains a brew recipe which handles dependency resolution, compilation and installation for FastSpar. To use the brew recipe, the below command can be run:
```bash
brew install https://raw.githubusercontent.com/scwatts/fastspar/master/scripts/fastspar.rb
```

### Compiling from source
Compilation from source requires these libraries and software:
```
c++11 (gcc-4.9.0+, clang-4.9.0+, etc)
OpenMP 4.0+
Gfortran
Armadillo 6.7+
LAPACK
OpenBLAS
GNU Scientific Library 2.1+
GNU getopt
GNU make
```

After meeting the above requirements, install FastSpar from source can be done with these commands:
```bash
git clone https://github.com/scwatts/fastspar.git
cd fastspar
./configure --prefix=/usr/
make
make install
```
Once done, the FastSpar executables can be run from the command line


## Usage
### Correlation inference
To run FastSpar, you must have absolute OTU counts in BIOM tsv format file (with no metadata). The `fake_data.txt` (from the original SparCC implementation) will be used as an example:

```bash
fastspar --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
```

The number of iterations (each iteration re-estimated fractions by drawing from a dirichlet distrubition) can also be changed. The number of exclusion iterations (the number of times highly correlation OTU pairs are discovered and excluded) can also be tweaked. Here's an example:

```bash
fastspar --iterations 500 --exclude_iterations 20 --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
```

Further, the minimum threshold to exclude correlated OTU pairs can be increased:
```bash
fastspar --threshold 0.2 --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
```


### Calculation of exact p-values
To calculate the p-value of the infered correlations, bootstraping can be performed. This process involves infering correlation from random permutations of the original OTU count data. The p-values are then calculated from the bootstrap correlations. In the below example, we calculate p-values from 1000 bootstrap correlations.


First we generate the 1000 boostrap counts:

```bash
mkdir bootstrap_counts
fastspar_bootstrap --otu_table fake_data.txt --number 1000 --prefix bootstrap_counts/fake_data
```

And then infer correlations for each bootstrap count (running in parallel with all processes available):

```bash
mkdir bootstrap_correlation
parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
```

From these correlations, the p-values are then calculated:
```bash
fastspar_exactpvalues --otu_table fake_data.txt --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv
```


### Threading
If FastSpar is compiled with OpenMP, threading can be used by invoking `--threads <thread_number>` at the command line for several tools:
```bash
fastspar --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv --iterations 10000 --threads 32
```


## License
This project is licensed under GNU GPLv3
