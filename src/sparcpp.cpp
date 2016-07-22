#include <string>
#include <time.h>
#include <vector>


#include <armadillo>
#include <gsl/gsl_randist.h>
#include <getopt.h>


#include "sparcpp.h"
#include "common.h"


// Initialise a SparCpp object (must be parsed a pointer to an OTU table struct and other paramters)
SparCpp::SparCpp(const OtuTable * otu_table, int iterations, int exclusion_iterations, int exclusion_threshold, gsl_rng * p_rng) {
    SparCpp::otu_table = otu_table;
    SparCpp::iterations = iterations;
    SparCpp::exclusion_iterations = exclusion_iterations;
    SparCpp::exclusion_threshold = exclusion_threshold;
    SparCpp::p_rng = p_rng;
}


// Run the SparCpp algorithm
void SparCpp::infer_correlation_and_covariance() {
    for (int i = 0; i < iterations; ++i) {
        // Create a SparCppIteration object
        SparCppIteration sparcpp_iteration(otu_table, exclusion_iterations, exclusion_threshold);
        std::cout << "Running iteration: " << i + 1 << std::endl;

        // TODO: Refactor as a method for SparCppIterations

        // STEP 1: Estimate component fractions and get log ratio variance
        sparcpp_iteration.estimate_component_fractions(p_rng);
        sparcpp_iteration.calculate_log_ratio_variance();

        // STEP 2: Calculate component variation
        sparcpp_iteration.calculate_component_variance();

        // STEP 3: Calculate correlation and covariance from basis variation and estimated fractions
        sparcpp_iteration.calculate_correlation_and_corvariance();

        // STEP 4: Exclude highly correlated pairs, repeating correlation/ covariance calculation each iteration
        for (int j = 0; j < exclusion_iterations; ++j) {
            // Find and exclude pairs
            sparcpp_iteration.find_and_exclude_pairs(exclusion_threshold);
            // Recalculate component variation after pair exclusion
            sparcpp_iteration.calculate_component_variance();
            // Recalculate correlation and covariance after pair exclusion
            sparcpp_iteration.calculate_correlation_and_corvariance();
        }

        // SparCC performs Pearson product-moment correlation coefficients (numpy.corrcoef) after apply a center
        // log ratio if the maximum correlation magnitude is greater than 1
        // TODO: Implement this (:
        if (arma::abs(sparcpp_iteration.basis_correlation).max() > 1.0){
            std::string error_str = "Input triggered condition to perform clr correlation, this is not yet implemented";
            std::cerr << error_str << std::endl;
	    exit(0);
        }

        // Finally add the calculated correlation and covariances to a running vector
        correlation_vector.push_back(sparcpp_iteration.basis_correlation);
        covariance_vector.push_back(sparcpp_iteration.basis_covariance.diag());
    }
}


// Initialise a SparCpp object (must be parsed a pointer to an OTU table struct and other paramters)
SparCppIteration::SparCppIteration(const OtuTable * otu_table, int exclusion_iterations, int exclusion_threshold) {
    SparCppIteration::otu_table = otu_table;
    SparCppIteration::exclusion_iterations = exclusion_iterations;
    SparCppIteration::exclusion_threshold = exclusion_threshold;

    // We also have to setup the mod matrix
    std::vector<double> mod_diag(otu_table->otu_number, otu_table->otu_number - 2);
    arma::Mat<double> mod = arma::diagmat((arma::Col<double>) mod_diag);
    SparCppIteration::mod = mod + 1;
}


void SparCppIteration::estimate_component_fractions(gsl_rng * p_rng) {
    // TODO: check if it's more efficient to gather fractions and then init arma::Mat (instead of init elements)
    // Estimate fractions by drawing from dirichlet distribution
    arma::Mat<double> temp_fractions(otu_table->sample_number, otu_table->otu_number);
    fractions.fill(1 / otu_table->otu_number);
    for(int i = 0; i < otu_table->sample_number; ++i) {
        // Get arma row and add pseudo count (then convert to double vector for rng function)
        arma::Row<double> row_pseudocount = otu_table->counts.row(i) + 1;
        std::vector<double> row_pseudocount_vector = arma::conv_to<std::vector<double>>::from(row_pseudocount);
        // Draw from dirichlet dist, storing results in theta double array
        size_t row_size = row_pseudocount_vector.size();
        double theta[row_size];
        // The function takes double arrays and it seems that you must pass the address of the first element to function
        gsl_ran_dirichlet(p_rng, row_size, &row_pseudocount_vector[0], theta);
        // Create arma::Row from double[] and update fractions row
        arma::Mat<double> estimated_fractions_row(theta, 1, otu_table->otu_number);
        temp_fractions.row(i) = estimated_fractions_row;
    }
    // Move temp_fractions to SparCppIteration.factions
    fractions = std::move(temp_fractions);
}


void SparCppIteration::calculate_log_ratio_variance() {
    // TODO: Given this is a square mat, check that we're required to iterate over all as in SparCC (thinking only half)
    // TODO: !IMPORTANT we're not assigning all elements in the variance matrix, is this correct?
    arma::Mat<double> temp_variance(fractions.n_cols, fractions.n_cols, arma::fill::zeros);
    for (unsigned int i = 0; i < fractions.n_cols - 1; ++i) {
        for (unsigned int j = i + 1; j < fractions.n_cols; ++j) {
            arma::Col<double> log_ratio_col = arma::log(fractions.col(i) / fractions.col(j));
            double col_variance = arma::var(log_ratio_col);
            temp_variance(i, j) = temp_variance(j, i) = col_variance;
        }
    }
    // Move temp_variance to SparCppIteration.variance
    variance = std::move(temp_variance);
}


void SparCppIteration::calculate_component_variance() {
    // NOTE: We make a copy of the variance matrix here to restrict modifications outside this function
    arma::Mat<double> variance_munge = variance;
    // If any pairs have been excluded, set variance to zero
    if (!excluded.empty()) {
        variance_munge((arma::Col<arma::uword>)excluded).fill(0.0);
    }
    arma::Col<double> variance_vector = arma::sum(variance_munge, 1);
    // Using double type as we'll need to get the inverse which fails when using an int mat
    basis_variance = mod.i() * variance_vector;
    // Set variances less than 0 to minimum variance
    basis_variance(arma::find(basis_variance < 0)).fill(1e-4);
}


void SparCppIteration::calculate_correlation_and_corvariance(){
    // TODO: Determine if basis_cor_el and basis_cov_el can all be calculated and then have arma::Mat initialised
    // Initialise matrices and set diagonals
    // TODO: We only need to set diagonal on final table (we otherwise set diag to zero then check for high correlates)
    std::vector<double> basis_cor_diag(otu_table->otu_number, 1);
    arma::Mat<double> temp_basis_correlation = arma::diagmat((arma::Col<double>) basis_cor_diag);
    arma::Mat<double> temp_basis_covariance = arma::diagmat((arma::Col<double>) basis_variance);
    // Calculate correlation and covariance for each element add set in basis matrices
    for (int i = 0; i < otu_table->otu_number - 1; ++i) {
        for (int j = i + 1; j < otu_table->otu_number; ++j) {
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
            temp_basis_correlation(i, j) = temp_basis_correlation(j, i) = basis_cor_el;
            temp_basis_covariance(i, j) = temp_basis_covariance(i, j) = basis_cov_el;
        }
    }
    // Move temp_basis_correlation and temp_basis_covariance to SparCppIteration
    basis_correlation = std::move(temp_basis_correlation);
    basis_covariance = std::move(temp_basis_covariance);
}


void SparCppIteration::find_and_exclude_pairs(float threshold) {
    // NOTE: BasisResults are passed by value as we need to makes changes only in this function
    basis_correlation.diag().zeros();
    basis_correlation = arma::abs(basis_correlation);
    // Set previously excluded correlations to zero
    basis_correlation((arma::Col<arma::uword>)excluded).fill(0.0);
    // Get all elements with the max value
    double max_correlate = basis_correlation.max();
    arma::Col<arma::uword> max_correlate_idx = arma::find(basis_correlation == max_correlate);
    // If max correlation is above a threshold, subtract one from the appropriate mod matrix positions
    if (max_correlate > threshold) {
        // For each max correlate pair
        for (arma::Col<arma::uword>::iterator it = max_correlate_idx.begin(); it != max_correlate_idx.end(); ++it) {
            // Substract from mod matrix
            unsigned int diagonal_idx = *it % otu_table->otu_number;
            mod.diag()[diagonal_idx] -= 1;
            // TODO: Check if it's quicker to select all elements at once and then subtract one
            mod(*it) -= 1;
            // Also add excluded indices to running list
            excluded.push_back(*it);
        }
    }
}


void SparCpp::calulcate_median_correlation_and_covariance() {
    // Get median of all i,j elements across the iterations for correlation
    // Add correlation matrices to arma Cube so that we can get views of all i, j of each matrix
    arma::Cube<double> correlation_cube(otu_table->otu_number, otu_table->otu_number, correlation_vector.size());
    correlation_cube.fill(0.0);

    // Fill cube with correlation matrix slices
    int cube_slice = 0;
    for (std::vector<arma::Mat<double>>::iterator it = correlation_vector.begin();
         it != correlation_vector.end(); ++it) {
        correlation_cube.slice(cube_slice) = *it;
        ++cube_slice;
    }

    // Get median value for each i, j element across all n iterations
    arma::Mat<double> temp_median_correlation(otu_table->otu_number, otu_table->otu_number);
    for (int i = 0; i < otu_table->otu_number; ++i) {
        for (int j = 0; j < otu_table->otu_number; ++j) {
            arma::Row<double> r = correlation_cube.subcube(arma::span(i), arma::span(j), arma::span());
            temp_median_correlation(i, j) = arma::median(r);
        }
    }

    // Get median for diagonal elements across iterations for covariance
    // Add covariance diagonals to arma Mat so that we can get row views for all i, j elements
    arma::Mat<double> covariance_diagonals(otu_table->otu_number, iterations);

    // Fill matrix will covariance diagonals
    int matrix_column = 0;
    for (std::vector<arma::Col<double>>::iterator it = covariance_vector.begin(); it != covariance_vector.end(); ++it) {
        covariance_diagonals.col(matrix_column) = *it;
        ++matrix_column;
    }

    // Get the median of each i, j element
    arma::Col<double> median_covariance_diag = arma::median(covariance_diagonals, 1);

    // Split into coordinate meshed grid
    arma::Mat<double> median_covariance_y(otu_table->otu_number, otu_table->otu_number);
    arma::Mat<double> median_covariance_x(otu_table->otu_number, otu_table->otu_number);
    median_covariance_y.each_col() = median_covariance_diag;
    median_covariance_x.each_row() = arma::conv_to<arma::Row<double>>::from(median_covariance_diag);

    // Calculate final median covariance
    arma::Mat<double> temp_median_covariance(otu_table->otu_number, otu_table->otu_number);
    temp_median_covariance = temp_median_correlation % arma::pow(median_covariance_x, 0.5) % arma::pow(median_covariance_y, 0.5);

    // Move the temp matrices to SparCpp
    median_correlation = std::move(temp_median_correlation);
    median_covariance = std::move(temp_median_covariance);
}


void printHelp() {
    std::cerr << "Program: SparCpp (c++ implementation of SparCC)" << std::endl;
    std::cerr << "Version: 0.1a" << std::endl;
    std::cerr << "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage:" << std::endl;
    std::cerr << "  sparcpp [options] --otu_table <of> --correlation <rf> --covariance <vf>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -t <of>, --otu_table <of>" << std::endl;
    std::cerr << "                OTU input table" << std::endl;
    std::cerr << "  -r <rf>, -correlation <rf>" << std::endl;
    std::cerr << "                Correlation output table" << std::endl;
    std::cerr << "  -v <vf>, --covariance <vf>" << std::endl;
    std::cerr << "                Covariance output table" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -h, --help    show this help message and exit" << std::endl;
    std::cerr << "  -i <int>, --iterations <int>" << std::endl;
    std::cerr << "                Number of interations to perform (50 default)" << std::endl;
    std::cerr << "  -x <int>, --exclusion_iterations <int>" << std::endl;
    std::cerr << "                Number of exclusion interations to perform (10 default)" << std::endl;
    std::cerr << "  -t <float>, --exclusion_iterations <float>" << std::endl;
    std::cerr << "                Correlation strength exclusion threshold (0.1 default)" << std::endl;

}


int main(int argc, char **argv) {
    // Set some default parameters
    int iterations = 20;
    int exclude_iterations = 10;
    float threshold = 0.1;

    // Declare some important variables
    std::string otu_filename;
    std::string correlation_filename;
    std::string covariance_filename;

    // Commandline arguments (for getlongtops)
    struct option long_options[] =
        {
            {"otu_table", required_argument, NULL, 't'},
            {"correlation", required_argument, NULL, 'r'},
            {"covariance", required_argument, NULL, 'v'},
            {"iterations", required_argument, NULL, 'i'},
            {"exclude_iterations", required_argument, NULL, 'x'},
            {"threshold", required_argument, NULL, 'e'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, 0, 0}
        };

    // Check if have an attemp at arguments, else print help
    if (argc < 2) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: options -t/--otu_table, -r/--correlation, and -v/--covariance are required" << std::endl;
        exit(0);
    }

    // Parse commandline arguments
    while (1) {
        int option_index = 0;
        int c;
        c = getopt_long (argc, argv, "ht:r:v:i:x:e:", long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch(c) {
            // TODO: do we need case(0)?
            case 't':
                otu_filename = optarg;
                break;
            case 'r':
                correlation_filename = optarg;
                break;
            case 'v':
                covariance_filename = optarg;
                break;
            case 'i':
                iterations = get_int_from_char(optarg);
                break;
            case 'x':
                exclude_iterations = get_int_from_char(optarg);
                break;
            case 'e':
                threshold = get_float_from_char(optarg);
                break;
            case 'h':
                printHelp();
                exit(0);
            default:
                exit(1);
        }
    }
    // Abort execution if given unknown arguments
    if (optind < argc){
        std::cerr << argv[0] << " invalid argument: " << argv[optind++] << std::endl;
    }
    // Make sure we have filenames
    if (otu_filename.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -t/--otu_table is required" << std::endl;
        exit(1);
    }
    if (correlation_filename.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -r/--correlation is required" << std::endl;
        exit(1);
    }
    if (covariance_filename.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -v/--covariance is required" << std::endl;
        exit(1);
    }
    // Ensure threshold is less than 100
    if (threshold > 1) {
        std::cerr << "Threshold cannot be greather than 1.0\n";
        exit(1);
    }
    // Check that the OTU file exists
    std::ifstream otu_file;
    otu_file.open(otu_filename);
    if (!otu_file.good()) {
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: OTU table file " << otu_filename << " does not exist" << std::endl;
        exit(1);
    }

    // Set up rng environment and seed
    const gsl_rng_type * rng_type;
    gsl_rng_env_setup();
    // gsl_rng_default is a global
    rng_type = gsl_rng_default;
    gsl_rng * p_rng = gsl_rng_alloc(rng_type);
    gsl_rng_set(p_rng, time(NULL));

    // Load the OTU table from file and construct count matrix
    struct OtuTable otu_table;
    otu_table.load_otu_file(otu_filename);

    // Initialise a SparCpp object
    SparCpp sparcpp(&otu_table, iterations, exclude_iterations, threshold, p_rng);

    // Run SparCpp iterations
    sparcpp.infer_correlation_and_covariance();

    // Calculate the final SparCpp correlation and covariances
    sparcpp.calulcate_median_correlation_and_covariance();

    // Write median correlation and covariance matrices
    write_out_square_otu_matrix(sparcpp.median_correlation, otu_table, correlation_filename);
    write_out_square_otu_matrix(sparcpp.median_covariance, otu_table, covariance_filename);

    // Free rng memory
    gsl_rng_free(p_rng);
}
