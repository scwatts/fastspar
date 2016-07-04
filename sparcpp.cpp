#include <fstream>
#include <iostream>
#include <sstream>
#include <time.h>
#include <vector>

#include <armadillo>
#include <gsl/gsl_randist.h>


struct OtuTable {
    std::vector<std::string> sample_names;
    std::vector<std::string> otu_ids;
    std::vector<int> otu_observations;
    int otu_number;
    int sample_number;
};


struct OtuTable LoadOtuFile(std::string filename) {
    // Used to store strings from file prior to assignment
    std::string line;
    std::string ele;
    std::stringstream line_stream;
    // Other variables
    struct OtuTable otu_table;
    int otu_number = 0;
    int sample_number = 0;
    bool id;
    // Open file stream
    std::ifstream otu_file;
    otu_file.open(filename);
    // Process header
    std::getline(otu_file, line);
    line_stream.str(line);
    // Iterate header columns
    while(std::getline(line_stream, ele, '\t')) {
        //TODO: Add assertion here
        // Skip the OTU_id column (first column)
        if (ele == "OTU_id") {
            continue;
        }
        // Store samples
        otu_table.sample_names.push_back(ele);
        ++sample_number;
    }
    // Process sample counts, need to get OTU IDS first
    while(std::getline(otu_file, line)) {
        // TODO: Is there an alternate design pattern to loop variables as below
        // Loop variables
        // (Re)sets variables for loop
        id = true;
        line_stream.clear();
        // Add current line to line stream and then split by tabs
        line_stream.str(line);
        while (std::getline(line_stream, ele, '\t')) {
            // Grab the OTU_id
            if (id) {
                otu_table.otu_ids.push_back(ele);
                id = false;
                continue;
            }
            // Add current element to OTU count after converting to int
            otu_table.otu_observations.push_back(std::stoi(ele));
        }
        ++otu_number;
    // TODO: Check if growing std::vector is sustainable for large tables
    }
    // Add counts to otu_table struct
    otu_table.otu_number = otu_number;
    otu_table.sample_number = sample_number;
    return otu_table;
}


arma::Mat<double> EstimateComponentFractions(const struct OtuTable& otu_table, arma::Mat<int>& counts) {
    // Set up rng environment and seed
    const gsl_rng_type * rng_type;
    gsl_rng_env_setup();
    // gsl_rng_default is a global
    rng_type = gsl_rng_default;
    gsl_rng * p_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(p_rng, time(NULL));
    // TODO: check if it's more efficient to gather fractions and then init arma::Mat (instead of init elements)
    // Estimate fractions by drawing from dirichlet distribution
    arma::Mat<double> fractions(otu_table.sample_number, otu_table.otu_number);
    for(int i = 0; i < otu_table.sample_number; ++i) {
        // Get arma row and add pseudo count (then convert to double vector for rng function)
        arma::Row<int> row_pseudocount = counts.row(i) + 1;
        std::vector<double> row_pseudocount_vector = arma::conv_to<std::vector<double>>::from(row_pseudocount);
        std::vector<double> * pr = &row_pseudocount_vector;
        // Draw from dirichlet dist, storing results in theta double array
        size_t row_size = row_pseudocount_vector.size();
        double theta[row_size];
        // The function takes double arrays and it seems that you must pass the address of the first element to function
        gsl_ran_dirichlet(p_rng, row_size, &row_pseudocount_vector[0], theta);
        // Create arma::Row from double[] and update fractions row
        arma::Mat<double> estimated_fractions_row(theta, 1, 50);
        fractions.row(i) = estimated_fractions_row;
    }
    // Free rng memory
    gsl_rng_free(p_rng);
    return fractions;
}


arma::Mat<double> FromFileGetComponentFractions(const struct OtuTable& otu_table, std::string filename) {
    // NOTE: This is a temp function to provide a consistent fraction component estimatation. In final build,
    // the function FileGetComponentFractions will be used.
    // Defining some variables
    std::string line;
    std::string ele;
    std::stringstream line_stream;
    std::vector<double> fraction_elements;
    // Open file stream
    std::ifstream fractions_file;
    fractions_file.open(filename);
    // Load all fraction elements into a vector
    while (std::getline(fractions_file, line)) {
        // Add current line to line stream and then split by tabs
        line_stream.clear();
        line_stream.str(line);
        while (std::getline(line_stream, ele, '\t')) {
            fraction_elements.push_back(std::stod(ele));
        }
    }
    arma::Mat<double> fractions(fraction_elements);
    fractions.reshape(otu_table.otu_number, otu_table.sample_number);
    arma::inplace_trans(fractions);
    return fractions;
}


arma::Mat<double> CalculateLogRatioVariance(const arma::Mat<double>& fractions) {
    // TODO: Given this is a square mat, check that we're required to iterate over all as in SparCC (thinking only half)
    arma::Mat<double> variance(fractions.n_cols, fractions.n_cols);
    for (int i = 0; i < fractions.n_cols - 1; ++i) {
        for (int j = i + 1; j < fractions.n_cols; ++j) {
            arma::Col<double> log_ratio_col = arma::log(fractions.col(i) / fractions.col(j));
            double col_variance = arma::var(log_ratio_col);
            variance(i, j) = variance(j, i) = col_variance;
        }
    }
    return variance;
}


int main() {
    // Load the OTU table from file
    std::string otu_filename;
    otu_filename = "fake_data.txt";
    struct OtuTable otu_table = LoadOtuFile(otu_filename);
    // Construct count matrix
    arma::Mat<int> counts(otu_table.otu_observations);
    counts.reshape(otu_table.sample_number, otu_table.otu_number);

    // STEP 1: Estimate component fractions and get log ratio variance
    // TODO: *** TEMP *** Will load in a pre-calculated component fraction matrix so that results/steps can be checked
    //arma::Mat<double> fractions = EstimateComponentFractions(&otu_table, &counts);
    arma::Mat<double> fractions = FromFileGetComponentFractions(otu_table, "fractions.tsv");
    arma::Mat<double> variance = CalculateLogRatioVariance(fractions);
}