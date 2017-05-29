#include "fastspar.h"


extern "C" void openblas_set_num_threads(int num_threads);


///////////////////////////////
//      FastSpar entry       //
///////////////////////////////

#if defined(FASTSPAR_CPACKAGE)
int main(int argc, char **argv) {
    // Starting message
    fprintf(stdout, "Starting FastSpar\n");

    // Set the number of threads openblas uses to 1
    openblas_set_num_threads(1);

    // Get commandline options
    FastsparOptions fastspar_options = get_commandline_arguments(argc, argv);

    // Load the OTU table from file and construct count matrix
    OtuTable otu_table;
    otu_table.load_otu_file(fastspar_options.otu_filename);


    // Check that OTUs have variance in their counts
    std::vector<int> invariant_otus;
    for (int i = 0; i < otu_table.otu_number; ++i) {
        arma::Col<double> unique_counts = arma::unique(otu_table.counts.col(i));
        if (unique_counts.n_elem == 1) {
            invariant_otus.push_back(i);
        }
    }

    // Print warnings
    if (! invariant_otus.empty()) {
        fprintf(stderr, "\nWarning: the following OTUs have only one unique permutation and it is recommended to remove them from this analysis:\n");
        for (auto &i : invariant_otus) {
            fprintf(stderr, "\t%s (row %d)\n", otu_table.otu_ids[i].c_str(), i);
        }
        fprintf(stderr, "\n");

        // Prompt user for continue
        if (! fastspar_options.assume_yes) {
            if (continue_exit_prompt() == false) {
                return 0;
            }
        }
    }


    // Initialise a FastSpar object
    FastSpar fastspar(&otu_table, fastspar_options.iterations, fastspar_options.exclude_iterations,
                      fastspar_options.threshold, fastspar_options.threads);


    // Run FastSpar iterations
    fprintf(stdout, "Running SparCC iterations\n");
    fastspar.infer_correlation_and_covariance();

    // Calculate the final FastSpar correlation and covariances
    fprintf(stdout, "Calculating final SparCC correlations and covariances\n");
    fastspar.calculate_median_correlation_and_covariance();

    // Write median correlation and covariance matrices
    fprintf(stdout, "Writing out median correlation and covariance matrices\n");
    write_out_square_otu_matrix(fastspar.median_correlation, otu_table, fastspar_options.correlation_filename);
    write_out_square_otu_matrix(fastspar.median_covariance, otu_table, fastspar_options.covariance_filename);

    return 0;
}
#endif




///////////////////////////////
//      Misc Function(s)     //
///////////////////////////////
// Continue or exit prompt
bool continue_exit_prompt() {
    // Define output and input strings
    std::string response;
    std::string prompt_message("Do you want to continue [Y/n]?");

    // Print message and get response
    fprintf(stderr, "%s", prompt_message.c_str());
    std::getline(std::cin, response);

    // Process response
    if (response.front() == 'Y' || response.front() == 'y' || response.empty()) {
        return true;
    } else if (response.front() == 'N' || response.front() == 'n') {
        return false;
    } else {
        fprintf(stderr, "Did not recognise %s (accepting 'y' or 'n')\n", response.c_str());
        return false;
    }
}




///////////////////////////////
//   Object Initialisation   //
///////////////////////////////

// Initialise a FastSpar object (must be parsed a pointer to an OTU table struct and other paramters)
FastSpar::FastSpar(const OtuTable *_otu_table, unsigned int _iterations, unsigned int _exclusion_iterations, unsigned int _exclusion_threshold, unsigned int _threads) {
    otu_table = _otu_table;
    iterations = _iterations;
    exclusion_iterations = _exclusion_iterations;
    exclusion_threshold = _exclusion_threshold;
    threads = _threads;

    // OpenMP function from omp.h. This sets the number of threads in a more reliable way but also ignores OMP_NUM_THREADS
    omp_set_num_threads(threads);
}


// Initialise a FastSpar object (must be parsed a pointer to an OTU table struct and other paramters)
FastSparIteration::FastSparIteration(const OtuTable *_otu_table, unsigned int _exclusion_iterations, unsigned int _exclusion_threshold) {
    otu_table = _otu_table;
    exclusion_iterations = _exclusion_iterations;
    exclusion_threshold = _exclusion_threshold;
    components_remaining = otu_table->otu_number;

    // We also have to setup the mod matrix
    std::vector<double> mod_diag(otu_table->otu_number, otu_table->otu_number - 2);
    arma::Mat<double> _mod = arma::diagmat((arma::Col<double>) mod_diag);
    mod = _mod + 1;
}




//////////////////////////////
// Core algorithm
//////////////////////////////

// Run the correlation algorithm
void FastSpar::infer_correlation_and_covariance() {
    // Make sure we have at least 4 components to run SparCC correlation calculation
    if (otu_table->otu_number < 4) {
        fprintf(stderr, "The SparCC algorithm requires at least 4 components to run\n");
        exit(0);
    }

    // Set up rng environment and seed
    FastSpar::p_rng = get_default_rng_handle();

#pragma omp parallel for schedule(static, 1)
    for (unsigned int i = 0; i < iterations; ++i) {
        // Create a FastSparIteration object
        FastSparIteration fastspar_iteration(otu_table, exclusion_iterations, exclusion_threshold);
        fprintf(stdout, "\tRunning iteration: %i\n", i + 1);

        // TODO: Refactor as a method for FastSparIterations

        // STEP 1: Estimate component fractions and get log ratio variance
        fastspar_iteration.estimate_component_fractions(p_rng);
        fastspar_iteration.calculate_fraction_log_ratio_variance();

        // STEP 2: Calculate basis variance
        fastspar_iteration.calculate_basis_variance();

        // STEP 3: Calculate correlation and covariance
        fastspar_iteration.calculate_correlation_and_covariance();

        // STEP 4: Exclude highly correlated pairs, repeating correlation/ covariance calculation for each iteration
        for (unsigned int j = 0; j < exclusion_iterations; ++j) {
            // The SparCC algorithm is only valid for 4 or more components, after exclusion if we
            // have left with fewer than 3 we exit early
            if (fastspar_iteration.components_remaining < 4) {
                fprintf(stderr, "Exclusion iterations ending early; we need at least 4 non-excluded pairs for SparCC\n");
                break;
            }

            // Find and exclude pairs
            fastspar_iteration.find_and_exclude_pairs(exclusion_threshold);

            // Recalculate component variation after pair exclusion
            fastspar_iteration.calculate_basis_variance();

            // Recalculate correlation and covariance after pair exclusion
            fastspar_iteration.calculate_correlation_and_covariance();
        }

        // SparCC performs Pearson product-moment correlation coefficients (numpy.corrcoef) after apply a center
        // log ratio if the maximum correlation magnitude is greater than 1
        // TODO: Implement this (:
        if (arma::abs(fastspar_iteration.basis_correlation).max() > 1.0){
            fprintf(stderr, "Input triggered condition to perform clr correlation, this is not yet implemented\n");
            exit(0);
        }

        // Finally add the calculated correlation and covariances to a running vector
#pragma omp critical
	{
        correlation_vector.push_back(fastspar_iteration.basis_correlation);
        covariance_vector.push_back(fastspar_iteration.basis_covariance.diag());
	}
    }

    // Free rng memory
    gsl_rng_free(p_rng);
}




////////////////////////////////
// Methods for core algorithm //
////////////////////////////////

// For each sample, estimate fractions by drawing from dirichlet distribution parameterised by sample pseudocount
void FastSparIteration::estimate_component_fractions(gsl_rng *p_rng) {
    // Set size of fraction matrix
    fraction_estimates.set_size(otu_table->sample_number, otu_table->otu_number);

    // Set row size
    size_t row_size = static_cast<size_t>(otu_table->otu_number);

    // Estimate fractions
    for(int i = 0; i < otu_table->sample_number; ++i) {
        // Get arma row and add pseudo count (then convert to double vector for rng function)
        arma::Row<double> row_pseudocount = otu_table->counts.row(i) + 1;

        // Draw from dirichlet dist, storing results in theta double array
        double *theta = new double[row_size];

        // The function takes double arrays and it seems that you must pass the address of the first element to function
        gsl_ran_dirichlet(p_rng, row_size, row_pseudocount.memptr(), theta);

        // Create arma::Row from double[] and update fractions row
        arma::Mat<double> estimated_fractions_row(theta, 1, otu_table->otu_number);
        fraction_estimates.row(i) = estimated_fractions_row;

        // Free dynamic memory
        delete[] theta;
    }
}


// Calculate the log-ration variance for each combination of sample fraction estimate
void FastSparIteration::calculate_fraction_log_ratio_variance() {
    // TODO: Test the amount of memory pre-computing the log matrix is consuming
    // Log fraction matrix and initialise a zero-filled matrix (diagonals must be initialised)
    arma::Mat<double> log_fraction_estimates = arma::log(fraction_estimates);
    arma::Mat<double> temp_fraction_variance(fraction_estimates.n_cols, fraction_estimates.n_cols, arma::fill::zeros);

    // Calculate log-ratio variance for fraction estimates
    for (unsigned int i = 0; i < fraction_estimates.n_cols - 1; ++i) {
         for (unsigned int j = i + 1; j < fraction_estimates.n_cols; ++j) {
             // Calculate variance of log fractions
             double col_variance = arma::var(log_fraction_estimates.col(i) - log_fraction_estimates.col(j));

             // Add to matrix
             temp_fraction_variance(i, j) = col_variance;
         }
    }

    // Reflect lower triangle to upper and move temp_fraction_variance to FastSparIteration.variance
    fraction_variance = arma::symmatu(temp_fraction_variance);
}


// Calculate the basis variance
void FastSparIteration::calculate_basis_variance() {
    // NOTE: We make a copy of the variance matrix here to restrict modifications outside this function
    arma::Mat<double> fraction_variance_munge = fraction_variance;

    // If any pairs have been excluded, set variance to zero
    if (!excluded_pairs.empty()) {
        fraction_variance_munge((arma::Col<arma::uword>)excluded_pairs).fill(0.0);
    }

    // Calculate the component variance
    arma::Col<double> component_variance= arma::sum(fraction_variance_munge, 1);

    // Solve Ax = b where A is the mod matrix and b is the component variance
    basis_variance = arma::solve(mod, component_variance, arma::solve_opts::fast);

    // Set variances less than 0 to minimum variance
    basis_variance(arma::find(basis_variance < 0)).fill(1e-4);
}


// Calculate the correlation and covariance
void FastSparIteration::calculate_correlation_and_covariance(){
    // Initialise matrices and vectors
    std::vector<double> basis_cor_diag(otu_table->otu_number, 1);
    arma::Mat<double> temp_basis_correlation = arma::diagmat((arma::Col<double>) basis_cor_diag);
    arma::Mat<double> temp_basis_covariance = arma::diagmat((arma::Col<double>) basis_variance);

    // Calculate correlation and covariance for each element add set in basis matrices
    for (int i = 0; i < otu_table->otu_number - 1; ++i) {
        for (int j = i + 1; j < otu_table->otu_number; ++j) {
            // Calculate cor and cov
            double basis_cov_el = 0.5 * (basis_variance(i) + basis_variance(j) - fraction_variance(i, j));
            double basis_cor_el = basis_cov_el / std::sqrt(basis_variance(i)) / std::sqrt(basis_variance(j));

            // Check if we got a valid correlation
            if (abs(basis_cor_el) > 1) {
                if (std::signbit(basis_cor_el) == 0) {
                    basis_cor_el = -1;
                }
                else if (std::signbit(basis_cor_el) == 1) {
                    basis_cor_el = 1;
                }
                // Recalculate correlation after setting cor
                basis_cov_el = basis_cor_el * std::sqrt(basis_variance(i)) * std::sqrt(basis_variance(j));
            }

            // TODO: We can simply keep the upper triangle and process later
            // Set basis_correlation and basis_covariance matrices
            temp_basis_correlation(i, j) = temp_basis_correlation(j, i) = basis_cor_el;
            temp_basis_covariance(i, j) = temp_basis_covariance(i, j) = basis_cov_el;
        }
    }

    // Move temp_basis_correlation and temp_basis_covariance to FastSparIteration
    basis_correlation = std::move(temp_basis_correlation);
    basis_covariance = std::move(temp_basis_covariance);
}


// Find the highest correlation and exclude this pair if correlation is above threshold
void FastSparIteration::find_and_exclude_pairs(float threshold) {
    // Set diagonal to zero as we're not interesting in excluding self-pairs and get absolute correlation value
    arma::Mat<double> basis_correlation_abs = arma::abs(basis_correlation);
    basis_correlation_abs.diag().zeros();

    // Set previously excluded correlations to zero
    basis_correlation_abs((arma::Col<arma::uword>)excluded_pairs).fill(0.0);

    // Get all elements with the max value
    double max_correlate = basis_correlation_abs.max();
    arma::Col<arma::uword> max_correlate_idx = arma::find(basis_correlation_abs == max_correlate);

    // If max correlation is above a threshold, subtract one from the appropriate mod matrix positions
    if (max_correlate > threshold) {
        // For each max correlate pair
        for (arma::Col<arma::uword>::iterator it = max_correlate_idx.begin(); it != max_correlate_idx.end(); ++it) {
            // Substract from mod matrix
            unsigned int i_idx = *it % otu_table->otu_number;
            unsigned int j_idx = std::floor(*it / otu_table->otu_number);

            // Diagonal elements
            --mod(i_idx, i_idx);
            --mod(j_idx, j_idx);

            // Non-diagonal elements
            --mod(i_idx, j_idx);
            --mod(j_idx, i_idx);

            // Also add excluded indices to running list
            excluded_pairs.push_back(i_idx * otu_table->otu_number + j_idx);
            excluded_pairs.push_back(j_idx * otu_table->otu_number + i_idx);

            // Finally decrease remaining component count
            --components_remaining;

            // The original implementation only excludes a single pair even if multiple pairs of same correlation are
            // identified (as determined by numpy.argmax; defined as the first occurence of the maximum value)
            break;
        }
    }
}


// Calculate the median correlation and covariance for each element across all iterations
void FastSpar::calculate_median_correlation_and_covariance() {
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

    // Move the temp matrices to FastSpar
    median_correlation = std::move(temp_median_correlation);
    median_covariance = std::move(temp_median_covariance);
}
