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


struct VarianceResults {
    arma::Col<double> basis_variance;
    arma::Mat<double> mod;
};


struct BasisResults {
    arma::Mat<double> basis_correlation;
    arma::Mat<double> basis_covariance;
};


struct OtuTable loadOtuFile(std::string filename) {
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


arma::Mat<double> estimateComponentFractions(const struct OtuTable& otu_table, arma::Mat<int>& counts) {
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


arma::Mat<double> fromFileGetComponentFractions(const struct OtuTable& otu_table, std::string filename) {
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


arma::Mat<double> calculateLogRatioVariance(const arma::Mat<double>& fractions) {
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


struct VarianceResults calculateComponentVariance(arma::Mat<double> variance,
                                                  arma::Mat<double>& mod,
                                                  std::vector<arma::uword>& excluded,
                                                  const struct OtuTable& otu_table) {
    // NOTE: The variance matrix is passed as value so that we can modify it only in this scope
    struct VarianceResults variance_results;
    // If any pairs have been excluded, set variance to zero
    if (!excluded.empty()) {
        variance((arma::Col<arma::uword>)excluded).fill(0.0);
    }
    arma::Col<double> variance_vector = arma::sum(variance, 1);
    // Using double type as we'll need to get the inverse which fails when using an int mat
    variance_results.mod = mod;
    variance_results.basis_variance = variance_results.mod.i() * variance_vector;
    // Set variances less than 0 to minimum variance
    variance_results.basis_variance(arma::find(variance_results.basis_variance < 0)).fill(1e-4);
    return variance_results;
}


struct BasisResults calculateCorAndCov(const arma::Mat<double>& variance, const arma::Col<double>& basis_variance,
                                       const struct OtuTable& otu_table) {
    // TODO: Determine if basis_cor_el and basis_cov_el can all be calculated and then have arma::Mat initialised
    // Initialise matrices and set diagonals
    struct BasisResults basis_results;
    // TODO: We only need to set diagonal on final table (we otherwise set diag to zero then check for high correlates)
    std::vector<double> basis_cor_diag(otu_table.otu_number, 1);
    basis_results.basis_correlation = arma::diagmat((arma::Col<double>) basis_cor_diag);
    basis_results.basis_covariance = arma::diagmat((arma::Col<double>) basis_variance);
    // Calculate correlation and covariance for each element add set in basis matrices
    for (int i = 0; i < otu_table.otu_number - 1; ++i) {
        for (int j = i + 1; j < otu_table.otu_number; ++j) {
            // Calculate cor and cov
            double basis_cov_el = 0.5 * (basis_variance(i) + basis_variance(j) - variance(i, j));
            double basis_cor_el = basis_cov_el / sqrt(basis_variance(i)) / sqrt(basis_variance(j));
            // Check if we got a valid correlation
            if (abs(basis_cor_el) > 1) {
                if (std::signbit(basis_cor_el) == 0) {
                    basis_cor_el = -1;
                }
                else if (std::signbit(basis_cor_el) == 1) {
                    basis_cor_el = 1;
                }
                // Recalculate correlation after setting cor
                double basis_cov_el = basis_cor_el * sqrt(basis_variance(i)) * sqrt(basis_variance(j));
            }
            // TODO: Check if we can avoid repetition here as well
            // Set basis_correlation and basis_covariance matrices
            basis_results.basis_correlation(i, j) = basis_results.basis_correlation(j, i) = basis_cor_el;
            basis_results.basis_covariance(i, j) = basis_results.basis_covariance(i, j) = basis_cov_el;
        }
    }
    return basis_results;
}


void findAndAddExcludedPairs(struct BasisResults basis_results, struct VarianceResults& variance_results,
                             const struct OtuTable& otu_table, std::vector<arma::uword>& excluded) {
    // NOTE: BasisResults are passed by value as we need to makes changes only in this function
    basis_results.basis_correlation.diag().zeros();
    basis_results.basis_correlation = arma::abs(basis_results.basis_correlation);
    // Set previously excluded correlations to zero
    basis_results.basis_correlation((arma::Col<arma::uword>)excluded).fill(0.0);
    // Get all elements with the max value
    double max_correlate = basis_results.basis_correlation.max();
    arma::Col<arma::uword> max_correlate_idx = arma::find(basis_results.basis_correlation == max_correlate);
    // If max correlation is above a threshold, subtract one from the appropriate mod matrix positions
    if (max_correlate > 0.1) {
        // For each max correlate pair
        for (arma::Col<arma::uword>::iterator it = max_correlate_idx.begin(); it != max_correlate_idx.end(); ++it) {
            // Substract from mod matrix
            unsigned int diagonal_idx = *it % otu_table.otu_number;
            variance_results.mod.diag()[diagonal_idx] -= 1;
            // TODO: Check if it's quicker to select all elements at once and then subtract one
            variance_results.mod(*it) -= 1;
            // Also add excluded indices to running list
            excluded.push_back(*it);
        }
    }
}


int main() {
    // Set some parameters
    const int iterations = 20;
    const int exclude_iterations = 10;

    // Define some out-of-loop variables
    std::vector<arma::Mat<double>> correlations;
    std::vector<arma::Mat<double>> covariances;

    // Load the OTU table from file
    std::string otu_filename;
    otu_filename = "fake_data.txt";
    struct OtuTable otu_table = loadOtuFile(otu_filename);
    // Construct count matrix
    arma::Mat<int> counts(otu_table.otu_observations);
    counts.reshape(otu_table.sample_number, otu_table.otu_number);

    for (int i = 0; i < iterations; ++i) {
        std::cout << "Running iteration: " << i << std::endl;
        // STEP 1: Estimate component fractions and get log ratio variance
        // TEMP: Will load in a pre-calculated component fraction matrix so that results/steps can be checked
        //arma::Mat<double> fractions = EstimateComponentFractions(&otu_table, &counts);
        arma::Mat<double> fractions = fromFileGetComponentFractions(otu_table, "fractions.tsv");
        arma::Mat<double> variance = calculateLogRatioVariance(fractions);

        // STEP 2: Calculate component variation
        std::vector<arma::uword> excluded;
        // Initialise mod matrix
        std::vector<double> mod_diag(otu_table.otu_number, otu_table.otu_number - 2);
        arma::Mat<double> mod = arma::diagmat((arma::Col<double>) mod_diag);
        mod = mod + 1;
        // And now calculate
        struct VarianceResults variance_results = calculateComponentVariance(variance, mod, excluded, otu_table);

        // STEP 3: Calculate correlation and covariance from basis variation and estimated fractions
        struct BasisResults basis_results = calculateCorAndCov(variance, variance_results.basis_variance, otu_table);

        // STEP 4: Exclude highly correlated pairs, repeating correlation/ covariance calculation each iteration
        std::cout << "Running exclusion: " << i << std::endl;
        for (int j = 0; j < exclude_iterations; ++j) {
            findAndAddExcludedPairs(basis_results, variance_results, otu_table, excluded);
            // NOTE: Make sure the previous iteration correlation basis is used
            // Recalculate component variation after pair exclusion
            variance_results = calculateComponentVariance(variance, variance_results.mod, excluded, otu_table);
            // Recalculate correlation and covariance after pair exclusion
            basis_results = calculateCorAndCov(variance, variance_results.basis_variance, otu_table);
        }

        // SparCC performs Pearson product-moment correlation coefficients (numpy.corrcoef) after apply a center
        // log ratio if the maximum correlation magnitude is greater than 1
        // TODO: Implement this (:
        if (arma::abs(basis_results.basis_correlation).max() > 1.0){
            std::string error_str = "Input triggered condition to perform clr correlation, this is not yet implemented";
            std::cerr << error_str << std::endl;
        }

        // Finally add the calculated correlation and covariances to a running vector
        //correlations.push_back(basis_results.basis_correlation);
        //covariances.push_back(basis_results.basis_covariance);
    }
}