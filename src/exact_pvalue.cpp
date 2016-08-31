#include <cmath>
#include <unordered_map>
#include <string>
#include <vector>


#include <armadillo>
#include <gsl/gsl_cdf.h>
#include <getopt.h>
#include <glob.h>


#include "exact_pvalue.h"
#include "common.h"


extern "C" {
    // TODO: See if we can forward declare in header also
    void c_gausq2(int* n, double* d, double* e, double* z, int* ierr);
}


// Collect bootstrap correlation paths by globbing
std::vector<std::string> get_bootstrap_correlation_paths(std::string& glob_path) {
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


void count_values_more_extreme(arma::Mat<double>& abs_observed_correlation,
                               arma::Mat<double>& abs_bootstrap_correlation,
                               arma::Mat<int>& extreme_value_counts) {
    // TODO: Check if this is faster than actually just processing the elements
    // Set diagonal to zero to avoid processing self pairs of OTUs
    abs_bootstrap_correlation.diag().zeros();
    // Find values more extreme than observed in bootstrap for each i,j element
    arma::Col<arma::uword> extreme_value_index = arma::find(abs_bootstrap_correlation >= abs_observed_correlation);
    // For each more extreme value, increment count in extreme_value_counts
    for (arma::Col<arma::uword>::iterator it = extreme_value_index.begin(); it != extreme_value_index.end(); ++it) {
        extreme_value_counts(*it) += 1;
    }

}


double factorial(double number) {
    double fact = 1;
    for (int i = 1; i <= number; ++i) {
        fact *= i;
    }
    return fact;
}


double calculate_possbile_otu_permutations(std::unordered_map<double, int>& count_frequency, int sample_number) {
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


double calculate_exact_pvalue(double otu_pair_possible_permutations, int& extreme_value_count, unsigned int& permutations) {
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


double calculate_pvalue_with_integral_estimate(double& otu_pair_possible_permutations, int& extreme_value_count, unsigned int& permutations) {
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


void printHelp() {
    std::cerr << "Program: SparCpp exact p-values (c++ implementation of SparCC)" << std::endl;
    std::cerr << "Version: 0.1" << std::endl;
    std::cerr << "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage:" << std::endl;
    std::cerr << "  exact_pavlues --otu_table <of> --correlation <rf> --prefix <p> --permutations <n> --output <o>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -t/--otu_table <of>   OTU input table used to generated correlations" << std::endl;
    std::cerr << "  -r/--correlation <rf> Correlation table generated by SparCpp" << std::endl;
    std::cerr << "  -p/--prefix <p>       Prefix output of bootstrap output files" << std::endl;
    std::cerr << "  -n/--permutations <n> Number of permutations/ bootstraps" << std::endl;
    std::cerr << "  -o/--outfile <o>      Output p-value matrix filename" << std::endl;

}


int main(int argc, char **argv) {
    // Define some variables
    std::string otu_filename;
    std::string correlation_filename;
    std::string bootstrap_prefix;
    std::string out_filename;
    unsigned int permutations = 0;

    // Commandline arguments (for getlongtops)
    struct option long_options[] =
        {
            {"otu_table", required_argument, NULL, 't'},
            {"prefix", required_argument, NULL, 'p'},
            {"correlation", required_argument, NULL, 'r'},
            {"permutations", required_argument, NULL, 'n'},
            {"outfile", required_argument, NULL, 'o'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, 0, 0}
        };

    // Check if have an attemp at arguments, else print help
    if (argc < 2) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: options -t/--otu_table, -r/--correlation, -p/--prefix, -n/--permutations and -o/--outfile are required" << std::endl;
        exit(0);
    }

    // Parse commandline arguments
    while (1) {
        int option_index = 0;
        int c;
        c = getopt_long (argc, argv, "ht:p:r:n:o:", long_options, &option_index);
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
            case 'p':
                bootstrap_prefix = optarg;
                break;
            case 'n':
                permutations = get_int_from_char(optarg);
                break;
            case 'o':
                out_filename = optarg;
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

    // Make sure we have filenames and parameters
    if (otu_filename.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -t/--otu_table is required" << std::endl;
        exit(1);
    }
    if (bootstrap_prefix.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -p/--prefix is required" << std::endl;
        exit(1);
    }
    if (correlation_filename.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -r/--correlation is required" << std::endl;
        exit(1);
    }
    if (permutations == 0) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -n/--permutations is required" << std::endl;
        exit(1);
    }
    if (out_filename.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -o/--output is required" << std::endl;
        exit(1);
    }

    // Check that the OTU file exists
    std::ifstream checkfile;
    checkfile.open(otu_filename);
    if (!checkfile.good()) {
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: OTU table file " << otu_filename << " does not exist" << std::endl;
        exit(1);
    }
    checkfile.close();
    // Check that the OTU file exists
    checkfile.open(correlation_filename);
    if (!checkfile.good()) {
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: correlation table file " << correlation_filename << " does not exist" << std::endl;
        exit(1);
    }
    checkfile.close();

    // Make sure our prefix has a trailing '*' for globbing
    if (bootstrap_prefix.back() != '*') {
        bootstrap_prefix += "*";
    }

    // Collect bootstrap correlation file paths and then make sure we have found the correct number
    std::vector<std::string> bs_cor_paths = get_bootstrap_correlation_paths(bootstrap_prefix);
    if (bs_cor_paths.size() != permutations) {
	    printf("ERROR: number of permutations, %i, isn't equal to the number of bootstrap correlations found, %zu\n",
			    permutations, bs_cor_paths.size());
	    exit(0);
    }

    // Read in otu tables (used to calculate total possible permutations)
    printf("Reading in OTU count table\n");
    struct OtuTable otu_table;
    otu_table.load_otu_file(otu_filename);

    // Read in observed correlation
    printf("Reading in observed correlations\n");
    arma::Mat<double> observed_correlation = load_correlation_file(correlation_filename);
    arma::Mat<double> abs_observed_correlation = arma::abs(observed_correlation);

    // Loop through vector of bootstrap correlations and count values for i, j elements that are more extreme than obs
    printf("Reading in %zu bootstrap correlations\n", bs_cor_paths.size());
    arma::Mat<int> extreme_value_counts(otu_table.otu_number, otu_table.otu_number, arma::fill::zeros);

#pragma omp parallel for
    for (int unsigned i = 0; i < bs_cor_paths.size(); ++i) {
        printf("\tBootstrap correlation %i: %s\n", i, bs_cor_paths[i].c_str());
        // Load the bootstrap correlation and get absolute values
        arma::Mat<double> bootstrap_correlation = load_correlation_file(bs_cor_paths[i]);
        arma::Mat<double> abs_bootstrap_correlation = arma::abs(bootstrap_correlation);

        // Count if bootstrap correlation is greater than observed correlation
        count_values_more_extreme(abs_observed_correlation, abs_bootstrap_correlation, extreme_value_counts);
    }

    // Calculate total possible permutations for each OTU
    printf("Calculating %i total permutations\n", otu_table.otu_number);
    arma::Col<double> possible_permutations(otu_table.otu_number, arma::fill::zeros);
    for (int i = 0; i < otu_table.otu_number; ++i) {
	printf("\tTotal permutation %i for %s\n", i, otu_table.otu_ids[i].c_str());
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

#pragma omp parallel for
    for (int i = 0; i < otu_table.otu_number; ++i) {
        printf("\tCalculating p-values for row %i with %s\n", i, otu_table.otu_ids[i].c_str());
        for (int j = 0; j < otu_table.otu_number; ++j) {
            // Get the total possible permutations between the OTU pair
            // TODO: Check if this is producing desired results
            double otu_pair_possible_permutations = possible_permutations(i) * possible_permutations(j);
            // The 'if' code in block below was ported from statmod::permp and statmod::gaussquad
            if (otu_pair_possible_permutations <= 10000 ) {
                // Exact p-value calculation
                // If fewer than 10000 possible permutations, we can safely cast double to int
                pvalues(i, j) = calculate_exact_pvalue((int)otu_pair_possible_permutations, extreme_value_counts(i, j), permutations);
            } else {
                // Integral approximation for p-value calculation
                pvalues(i, j) = calculate_pvalue_with_integral_estimate(otu_pair_possible_permutations, extreme_value_counts(i, j), permutations);
            }
        }
    }

    // Write out p-values
    printf("Writing out p-values\n");
    write_out_square_otu_matrix(pvalues, otu_table, out_filename);
}
