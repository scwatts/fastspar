# FastSpar
FastSpar is a C++ implementation of the [SparCC](https://bitbucket.org/yonatanf/sparcc) algorithm [published](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687) by Jonathan Friedman and Eric Alm. FastSpar is up to several thousand times faster than the original Python2 implementation and uses much less memory. FastSpar offers threading support and an *p*-value estimator which accounts for the possibility of repeating data permutations (see [this](https://arxiv.org/pdf/1603.05766.pdf) paper for further details). Code for the *p*-value estimator was adopted from the permp function from the [statmod](https://cran.r-project.org/web/packages/statmod/index.html) R package and from an example provided by [Scott Ritchie](https://github.com/sritchie73).


## Citation
If you use this tool, please cite the FastSpar preprint and original SparCC paper:
* [Watts, S.C., Ritchie, S.C., Inouye, M. & Holt, K.E. (2018). FastSpar: Rapid and scalable correlation estimation for compositional data. bioRxiv doi: 10.1101/222190.](https://www.biorxiv.org/content/early/2018/03/03/272583)
* [Friedman, J. & Alm, E.J. (2017). Inferring correlation networks from genomic survey data. PLoS Comput. Biol. 8, e1002687.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687)


## Installing
FastSpar can be installed via prebuilt distribution packages, from source, or using statically linked binaries.


### GNU/Linux
For most 64-bit linux distributions (e.g. Ubuntu, Debian, RedHat, etc) the easiest way to obtain FastSpar is via statically compiled binaries on the GitHub release page. These binaries can be downloaded and run immediately without any setup as they have no dependencies.

Alternative, there are also Ubuntu and Debian packages available.


### macOS
On macOS, FastSpar can be installed using homebrew. The GitHub repository contains a brew recipe which handles dependency resolution, compilation and installation for FastSpar. To use the brew recipe (on systems with homebrew installed), the below command can be run:
```bash
brew install https://raw.githubusercontent.com/scwatts/fastspar/master/scripts/fastspar.rb
```

### Compiling from source
Compilation from source requires these libraries and software:
```
C++11 (gcc-4.9.0+, clang-4.9.0+, etc)
OpenMP 4.0+
Gfortran
Armadillo 6.7+
LAPACK
OpenBLAS
GNU Scientific Library 2.1+
GNU getopt
GNU make
```

After meeting the above requirements, installing FastSpar from source can be done with these commands:
```bash
git clone https://github.com/scwatts/fastspar.git
cd fastspar
./configure --prefix=/usr/
make
make install
```
Once completed, the FastSpar executables can be run from the command line....


## Usage
### Correlation inference
To run FastSpar, you must have absolute OTU counts in BIOM tsv format file (with no metadata). The `fake_data.txt` (from the original SparCC implementation) will be used as an example:

```bash
fastspar --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
```

The number of iterations (each iteration re-estimated fractions by drawing from a dirichlet distrubition) can also be changed. The number of exclusion iterations (the number of times highly correlation OTU pairs are discovered and excluded) can also be tweaked. Here's an example:

```bash
fastspar --iterations 50 --exclude_iterations 20 --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
```

Further, the minimum threshold to exclude correlated OTU pairs can be increased:
```bash
fastspar --threshold 0.2 --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
```


### Calculation of exact *p*-values
To calculate the *p*-value of the infered correlations, bootstraping can be performed. This process involves infering correlation from random permutations of the original OTU count data. The *p*-values are then calculated from the bootstrap correlations. In the below example, we calculate *p*-values from 1000 bootstrap correlations.


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

From these correlations, the *p*-values are then calculated:
```bash
fastspar_exactpvalues --otu_table fake_data.txt --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_fake_data_ --permutations 1000 --outfile pvalues.tsv
```


### Threading
If FastSpar is compiled with OpenMP, threading can be used by invoking `--threads <thread_number>` at the command line for several tools:
```bash
fastspar --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv --iterations 50 --threads 10
```

## License
GNU GPLv3
