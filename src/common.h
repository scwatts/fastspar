#ifndef __COMMON_H__
#define __COMMON_H__


#include <fstream>
#include <iomanip>
#include <string>
#include <thread>
#include <vector>


#include <gsl/gsl_randist.h>
#include "omp.h"


#include "config.h"
#include "local_config.h"


#if defined(FASTSPAR_RPACKAGE)
    #include "RcppArmadillo.h"
#else
    #include "armadillo"
#endif


// Object to hold information about an OTU table and method to load from file
struct OtuTable {
    std::vector<std::string> sample_names;
    std::vector<std::string> otu_ids;
    arma::Mat<float> counts;
    int otu_number = 0;
    int sample_number = 0;

    // Load OTU table
    void load_otu_file(std::string filename);
};


// Load a correlation (or covariance table) from file
arma::Mat<float> load_correlation_file(std::string &filename);

// Save an square OTU matrix (e.g. correlation matrix) to file
void write_out_square_otu_matrix(arma::Mat<float> &matrix, OtuTable &otu_table, std::string filename);

// Convert character to integer (for commandline argument parsing)
int int_from_optarg(const char *optarg);

// Convert character to float (for commandline argument parsing)
float float_from_optarg(const char *optarg);

void directory_exists(const std::string &optarg);

#endif
