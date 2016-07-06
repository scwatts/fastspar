# SparCpp
A c++ implementation of the SparCC algorithm published here: Friedman, J. & Alm, E. J. Inferring correlation networks from genomic survey data. PLoS Comput. Biol. 8, e1002687 (2012). There is an approximately 50x speed up compared to the original Python2 implementation.

Much of the implmentation was ported from original SparCC implmentation found here: https://bitbucket.org/yonatanf/sparcc

## Installing
SparCpp has been written using c++ with Armadillo, GNU Scientific Library (GSL), and getopt. Compilation will require these libraries.


### Prerequisities
For compilation the following is required:
```
C++11
Armadillo 6.7+
GNU Scientific Library 2.1+
GNU getopt
```

### Installing
Compile SparCpp, using O2 optimisation and linking against Armadillo and GSL:
```bash
g++ -Wall -O2 -std=c++11 -larmadillo -lgsl -lgslcblas -lm -o sparcpp sparcpp.cpp
```

## Usage
To run SparCpp, you must have a BIOM tsv format file (with no metadata). The `fake_data.txt` (from the original SparCC implementation) will be used as an example:

```bash
./sparcpp --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
```

The number of iterations (each iteration re-estimated fractions by drawing from a dirichlet distrubition) can also be changed. The number of exclusion iterations (the number of times highly correlation OTU pairs are discovered and excluded) can also be tweaked. Here's an example:

```bash
./sparcpp --iterations 500 --exclude_iterations 20 --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
```

Further, the minimum threshold to exclude correlated OTU pairs can be increased:
```bash
./sparcpp --threshold 0.2 --otu_table fake_data.txt --correlation median_correlation.tsv --covariance median_covariance.tsv
```

## License
This project is licensed under the GNU GPLv3 Licence
