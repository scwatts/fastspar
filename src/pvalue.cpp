#include "pvalue.h"


extern "C" {
    // TODO: See if we can forward declare in header also
    void c_gausq2(int* n, double* d, double* e, double* z, int* ierr);
}


// Collect bootstrap correlation paths by globbing
std::vector<std::string> get_bootstrap_correlation_paths(std::string glob_path) {
    // Perform glob
    glob_t glob_results;
    glob(glob_path.c_str(), GLOB_TILDE, NULL, &glob_results);

    // Put globbed paths onto string vector
    std::vector<std::string> bootstrap_correlations;
    for(unsigned int i = 0; i < glob_results.gl_pathc; ++i) {
        bootstrap_correlations.push_back(std::string(glob_results.gl_pathv[i]));
    }

    // Free glob memory and return string vector
    globfree(&glob_results);
    return bootstrap_correlations;
}


void count_values_more_extreme(arma::Mat<double> &abs_observed_correlation,
                               arma::Mat<double> &abs_bootstrap_correlation,
                               arma::Mat<int> &extreme_value_counts) {
    // TODO: Check if this is faster than actually just processing the elements
    // Set diagonal to zero to avoid processing self pairs of OTUs
    abs_bootstrap_correlation.diag().zeros();
    // Find values more extreme than observed in bootstrap for each i,j element
    arma::Col<arma::uword> extreme_value_index = arma::find(abs_bootstrap_correlation >= abs_observed_correlation);
    // For each more extreme value, increment count in extreme_value_counts
    for (arma::Col<arma::uword>::iterator it = extreme_value_index.begin(); it != extreme_value_index.end(); ++it) {
#pragma omp atomic
        ++extreme_value_counts(*it);
    }

}


double factorial(double number) {
    double fact = 1;
    for (int i = 1; i <= number; ++i) {
        fact *= i;
    }
    return fact;
}


double calculate_possbile_otu_permutations(std::unordered_map<double, int> &count_frequency, int sample_number) {
    // The total permutations for a single OTU can be calculated by factorial division. We try to
    // simplify it here (ported from R code authored by Scott Ritchie)
    int max = 0;
    double numerator = 1;
    double denominator = 1;
    for (std::unordered_map<double, int>::iterator it = count_frequency.begin(); it != count_frequency.end(); ++it) {
        // Factorial of 1 is 1
        if (it->second == 1){
            continue;
        }
        // This part simplifies factorial division - we cancel largest/remove denominator
        if (it->second > max) {
            // If previous max was greater than 0, multiply denominator by factorial(max)
            if (max > 0){
                denominator *= factorial(max);
            }
            // Assigned new max
            max = it->second;
        } else {
                denominator *= factorial(it->second);
        }
    }
    // Calculate the simpilfied numerator
    for (int i = sample_number; i > max; --i) {
        numerator *= i;
    }
    // Finally calculate possible permutations for the OTU
    // TODO: This method can produce Inf, is this okay?
    if (std::isinf(numerator) || std::isinf(denominator)) {
        return std::numeric_limits<double>::infinity();
    } else {
        return numerator / denominator;
    }
}


double calculate_exact_pvalue(double otu_pair_possible_permutations, int extreme_value_count, unsigned int permutations) {
    // Function adapted and ported from statmod::permp
    // This cast is messy (also can't pass otu_pair as double reference for some reason)
    double prob[(int)otu_pair_possible_permutations];
    double prob_binom_sum = 0;
    for (int i = 0; i < otu_pair_possible_permutations; ++i) {
        prob[i] = (double)(i + 1) / otu_pair_possible_permutations;
    }
    for (int i = 0; i < otu_pair_possible_permutations; ++i) {
        prob_binom_sum += gsl_cdf_binomial_P(extreme_value_count, prob[i], permutations);
    }
    // Return p-value
    return prob_binom_sum / otu_pair_possible_permutations;
}


double calculate_pvalue_with_integral_estimate(double otu_pair_possible_permutations, int extreme_value_count, unsigned int permutations) {
    // Function adapted and ported from statmod::permp and statmod::gaussquad
    // TODO: See if there is a better way to init array elements w/o hard coding
    // Start statmod::gaussquad port
    // Declare variables for c_gausq2 call
    int n = 128;
    double a[n];
    double b[n];
    double z[n];
    int ierr = 0;

    // Initialise array elements for c_gausq2 call
    for (int i = 0; i < n; ++i) {
        a[i] = 0;
        z[i] = 0;
    }
    for (int i = 0; i < (n - 1); ++i){
        int i1 = i + 1;
        b[i] = i1 / sqrt(4 * pow(i1, 2) - 1);
    }
    b[n] = 0;
    z[0] = 1;

    // Make call to Fortran subroutine. Variables a, b, z and ierr are modified
    c_gausq2(&n, a, b, z, &ierr);

    // Further cals
    double u = 0.5 / otu_pair_possible_permutations;
    double weights[n];
    double nodes[n];
    for (int i = 0; i < n; ++i) {
        weights[i] = pow(z[i], 2);
    }
    for (int i = 0; i < n; ++i) {
        nodes[i] = u * (a[i] + 1) / 2;
    }
    // End statmod::gaussquad port

    // Start statmod::permp port
    double weight_prob_product_sum = 0;
    for (int i = 0; i < n; ++i) {
        weight_prob_product_sum += gsl_cdf_binomial_P(extreme_value_count, nodes[i], permutations) * weights[i];
    }
    double integral = 0.5 / (otu_pair_possible_permutations * weight_prob_product_sum);
    // TODO: Check if the double cast correctly done (it is required but maybe adding 1.0 instead of 1 is sufficient)
    // End statmod::permp port with p-value return
    return ((double)extreme_value_count + 1) / ((double)permutations + 1) - integral;
}


arma::Mat<double> calculate_pvalues(OtuTable &otu_table, arma::Mat<double> &observed_correlation, std::vector<std::string> &bootstrap_correlation_fps, unsigned int permutations, unsigned int threads) {
    // Get absolute correlations
    arma::Mat<double> abs_observed_correlation = arma::abs(observed_correlation);

    // OpenMP function from omp.h. This sets the number of threads in a more reliable way but also ignores OMP_NUM_THREADS
    omp_set_num_threads(threads);

    // Loop through vector of bootstrap correlations and count values for i, j elements that are more extreme than obs
    printf("Reading in %zu bootstrap correlations\n", bootstrap_correlation_fps.size());
    arma::Mat<int> extreme_value_counts(otu_table.otu_number, otu_table.otu_number, arma::fill::zeros);

#pragma omp parallel for schedule(static, 1)
    for (int unsigned i = 0; i < bootstrap_correlation_fps.size(); ++i) {
        printf("\tBootstrap correlation %i: %s\n", i, bootstrap_correlation_fps[i].c_str());
        // Load the bootstrap correlation and get absolute values
        arma::Mat<double> bootstrap_correlation = load_correlation_file(bootstrap_correlation_fps[i]);
        arma::Mat<double> abs_bootstrap_correlation = arma::abs(bootstrap_correlation);

        // Count if bootstrap correlation is greater than observed correlation
        count_values_more_extreme(abs_observed_correlation, abs_bootstrap_correlation, extreme_value_counts);
    }

    // Calculate total possible permutations for each OTU
    printf("Calculating %i total permutations\n", otu_table.otu_number);
    arma::Col<double> possible_permutations(otu_table.otu_number, arma::fill::zeros);
    for (int i = 0; i < otu_table.otu_number; ++i) {
        printf("\tTotal permutation %i of %d\n", i, otu_table.otu_number);
        // First we need to get the frequency of each count for an OTU across all samples
        // TODO: Check that after changing from int counts to double (for corrected OTU counts) that this isn't borked
        // Main concern that equality will not be true in some cases where they would be otherwise due to float error
        std::unordered_map<double, int> count_frequency;
        for (arma::Mat<double>::col_iterator it = otu_table.counts.begin_col(i); it != otu_table.counts.end_col(i); ++it) {
            ++count_frequency[*it];
        }

        // Call function to calculate possible permutations and store
        possible_permutations(i) = calculate_possbile_otu_permutations(count_frequency, otu_table.sample_number);
    }

    // Calculate p-values; loop through each i, j element in the extreme counts matrix
    printf("Calculating the %i p-values\n", otu_table.otu_number*otu_table.otu_number);
    arma::Mat<double> pvalues(otu_table.otu_number, otu_table.otu_number, arma::fill::zeros);

#pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < otu_table.otu_number; ++i) {
        printf("\tCalculating p-values for row %i of %d\n", i, otu_table.otu_number);
        for (int j = 0; j < otu_table.otu_number; ++j) {
            // Get the total possible permutations between the OTU pair
            // TODO: Check if this is producing desired results
            double otu_pair_possible_permutations = possible_permutations(i) * possible_permutations(j);
            // The 'if' code in block below was ported from statmod::permp and statmod::gaussquad
            if (otu_pair_possible_permutations <= 10000 ) {
                // Exact p-value calculation
                // If fewer than 10000 possible permutations, we can safely cast double to int
#pragma omp atomic write
                pvalues(i, j) = calculate_exact_pvalue((int)otu_pair_possible_permutations, extreme_value_counts(i, j), permutations);
            } else {
                // Integral approximation for p-value calculation
#pragma omp atomic write
                pvalues(i, j) = calculate_pvalue_with_integral_estimate(otu_pair_possible_permutations, extreme_value_counts(i, j), permutations);
            }
        }
    }

    return pvalues;
}


#if defined(FASTSPAR_CPACKAGE)
int main(int argc, char **argv) {
    // Get commandline arguments
    PvalOptions pval_options = get_commandline_arguments(argc, argv);

    // Collect bootstrap correlation file paths and then make sure we have found the correct number
    std::vector<std::string> bs_cor_paths = get_bootstrap_correlation_paths(pval_options.bootstrap_prefix);
    if (bs_cor_paths.size() != pval_options.permutations) {
	    printf("ERROR: number of permutations, %i, isn't equal to the number of bootstrap correlations found, %zu\n",
			    pval_options.permutations, bs_cor_paths.size());
	    exit(0);
    }

    // Read in otu tables (used to calculate total possible permutations)
    printf("Reading in OTU count table\n");
    OtuTable otu_table;
    otu_table.load_otu_file(pval_options.otu_filename);

    // Read in observed correlation
    printf("Reading in observed correlations\n");
    arma::Mat<double> observed_correlation = load_correlation_file(pval_options.correlation_filename);

    // Calulate p-values
    arma::Mat<double> pvalues = calculate_pvalues(otu_table, observed_correlation, bs_cor_paths,
                                                    pval_options.permutations, pval_options.threads);
    // Write out p-values
    printf("Writing out p-values\n");
    write_out_square_otu_matrix(pvalues, otu_table, pval_options.out_filename);

    return 0;
}
#endif
