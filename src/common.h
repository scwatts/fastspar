#ifndef __COMMON_H__
#define __COMMON_H__


#include <string>
#include <vector>


#include <armadillo>


// Object to hold information about an OTU table and method to load from file
struct OtuTable {
    std::vector<std::string> sample_names;
    std::vector<std::string> otu_ids;
    arma::Mat<double> counts;
    int otu_number = 0;
    int sample_number = 0;

    // Load OTU table
    void load_otu_file(std::string filename);
};


// Load a correlation (or covariance table) from file
arma::Mat<double> load_correlation_file(std::string& filename);

// Save an square OTU matrix (e.g. correlation matrix) to file
void write_out_square_otu_matrix(arma::Mat<double> & matrix, OtuTable & otu_table, std::string filename);


// Convert character to integer (for commandline argument parsing)
int get_int_from_char(const char* optarg);


// Convert character to float (for commandline argument parsing)
float get_float_from_char(const char* optarg);


#endif
