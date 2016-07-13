#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <vector>

#include <armadillo>
#include <gsl/gsl_cdf.h>
#include <getopt.h>
#include <glob.h>


extern "C" {
    void c_gausq2(int* n, double* d, double* e, double* z, int* ierr);
}


struct OtuTable {
    std::vector<std::string> sample_names;
    std::vector<std::string> otu_ids;
    std::vector<int> otu_observations;
    int otu_number;
    int sample_number;
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
        if (ele == "OTU_id" || ele == "#OTU ID") {
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


std::vector<std::string> getBootstrapCorrelationPaths(std::string& glob_path) {
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


arma::Mat<double> loadCorrelation(std::string& filename, struct OtuTable& otu_table) {
    // Used to store strings from file prior to matrix construction
    std::string line;
    std::string ele;
    std::stringstream line_stream;
    std::vector<double> correlations_vector;
    correlations_vector.reserve(otu_table.otu_number * otu_table.otu_number);
    // Other variables
    bool id;
    // Open file stream
    std::ifstream cor_file;
    cor_file.open(filename);
    // Skip header, order SHOULD be the same as input OTU table
    std::getline(cor_file, line);
    line_stream.str(line);
    // Process correlation elements
    while(std::getline(cor_file, line)) {
        // (Re)sets variables for loop
        id = true;
        line_stream.clear();
        // Add current line to line stream and then split by tabs
        line_stream.str(line);
        while (std::getline(line_stream, ele, '\t')) {
            // Skip the OTU ID column
            if (id) {
                id = false;
                continue;
            }
            // Add current element to correlation mat after converting to double
            correlations_vector.push_back(std::stod(ele));
        }
    }
    // Construct matrix and return it
    arma::Mat<double> correlations(correlations_vector);
    correlations.reshape(otu_table.otu_number, otu_table.otu_number);
    return correlations;
}


void countValuesMoreExtreme(arma::Mat<double>& abs_observed_correlation,
                            arma::Mat<double>& abs_bootstrap_correlation,
                            arma::Mat<int>& extreme_value_counts) {
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


double calculatePossiblePermutationsForOTU(std::unordered_map<int, int>& count_frequency, struct OtuTable& otu_table) {
    // The total permutations for a single OTU can be calculated by factorial division. We try to
    // simplify it here (ported from R code authored by Scott Ritchie)
    int max = 0;
    double numerator = 1;
    double denominator = 1;
    for (std::unordered_map<int, int>::iterator it = count_frequency.begin(); it != count_frequency.end(); ++it) {
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
    for (int i = otu_table.sample_number; i > max; --i) {
        numerator *= i;
    }
    // Finally calculate possible permutations for the OTU
    // TODO: This method can produce Inf, is this okay?
    return numerator / denominator;
}


double calculateExactPvalue(double otu_pair_possible_permutations, arma::Mat<int>& extreme_value_counts,
                                int& permutations, int& i, int& j) {
    // Function adapted and ported from statmod::permp
    // This cast is messy (also can't pass otu_pair as double reference for some reason)
    double prob[(int)otu_pair_possible_permutations];
    double prob_binom_sum = 0;
    for (int i = 0; i < otu_pair_possible_permutations; ++i) {
        prob[i] = (double)(i + 1) / otu_pair_possible_permutations;
    }
    for (int i = 0; i < otu_pair_possible_permutations; ++i) {
        prob_binom_sum += gsl_cdf_binomial_P(extreme_value_counts(i, j), prob[i], permutations);
    }
    // Return p-value
    return prob_binom_sum / otu_pair_possible_permutations;
}


double calculatePvalueWithIntegralEstimate(double& otu_pair_possible_permutations, arma::Mat<int>& extreme_value_counts,
                                            int& permutations, int& i, int& j) {
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
    for (int k = 0; k < n; ++k) {
        weight_prob_product_sum += gsl_cdf_binomial_P(extreme_value_counts(i, j), nodes[k], permutations) * weights[i];
    }
    double integral = 0.5 / (otu_pair_possible_permutations * weight_prob_product_sum);
    // TODO: Check if the double cast correctly done (it is required but maybe adding 1.0 instead of 1 is sufficient)
    // End statmod::permp port with p-value return
    return ((double)extreme_value_counts(i, j) + 1) / ((double)permutations + 1) - integral;
}


void writeOutMatrix(arma::Mat<double>& matrix, std::string out_filename, struct OtuTable& otu_table) {
    // Get stream handle
    std::ofstream outfile;
    outfile.open(out_filename);
    // Write out header
    outfile << "#OTU ID";
    for (std::vector<std::string>::iterator it = otu_table.otu_ids.begin(); it != otu_table.otu_ids.end(); ++it) {
        outfile << "\t" << *it;
    }
    outfile << std::endl;
    // Write out values
    for (unsigned int i = 0; i < matrix.n_rows; ++i) {
        for (unsigned int j = 0; j < matrix.n_cols; ++j) {
            // Write the OTU id as first field in row
            if (j == 0) {
                outfile << otu_table.otu_ids[i];
            }
            outfile << "\t" << std::fixed << std::setprecision(5) << matrix(i, j);
        }
        outfile << std::endl;
    }
}


int getIntFromChar(const char* optarg) {
    // Check at most the first 8 characters are numerical
    std::string optstring(optarg);
    std::string string_int = optstring.substr(0, 8);
    for (std::string::iterator it = string_int.begin(); it != string_int.end(); ++it) {
        if (!isdigit(*it)) {
            std::cerr << "This doesn't look like a number: " << optarg << std::endl;;
            exit(1);
        }
    }
    return std::atoi(string_int.c_str());
}


void printHelp() {
    std::cerr << "Program: SparCpp exact p-values (c++ implementation of SparCC)" << std::endl;
    std::cerr << "Version: 0.1a" << std::endl;
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
    int permutations = 0;

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
                permutations = getIntFromChar(optarg);
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


    // Read in otu tables (used to calculate total possible permutations)
    struct OtuTable otu_table = loadOtuFile(otu_filename);
    arma::Mat<int> counts(otu_table.otu_observations);
    counts.reshape(otu_table.sample_number, otu_table.otu_number);

    // Read in observed correlation
    arma::Mat<double> observed_correlation = loadCorrelation(correlation_filename, otu_table);
    arma::Mat<double> abs_observed_correlation = arma::abs(observed_correlation);

    // Loop through vector of bootstrap correlations and count values for i, j elements that are more extreme than obs
    arma::Mat<int> extreme_value_counts(otu_table.otu_number, otu_table.otu_number, arma::fill::zeros);
    std::vector<std::string> bs_cor_paths = getBootstrapCorrelationPaths(bootstrap_prefix);
    for (std::vector<std::string>::iterator it = bs_cor_paths.begin(); it != bs_cor_paths.end(); ++it) {
        arma::Mat<double> bootstrap_correlation = loadCorrelation(*it, otu_table);
        arma::Mat<double> abs_bootstrap_correlation = arma::abs(bootstrap_correlation);
        countValuesMoreExtreme(abs_observed_correlation, abs_bootstrap_correlation, extreme_value_counts);
    }

    // Calculate total possible permutations for each OTU
    arma::Col<double> possible_permutations(otu_table.otu_number);
    for (int i = 0; i < otu_table.otu_number; ++i) {
        // First we need to get the frequency of each count for an OTU across all samples
        std::unordered_map<int, int> count_frequency;
        for (arma::Mat<int>::col_iterator it = counts.begin_col(i); it != counts.end_col(i); ++it) {
            ++count_frequency[*it];
        }
        // Call function to calculate possible permutations and store
        possible_permutations(i) = calculatePossiblePermutationsForOTU(count_frequency, otu_table);
    }


    // Calculate p-values; loop through each i, j element in the extreme counts matrix
    arma::Mat<double> pvalues(otu_table.otu_number, otu_table.otu_number, arma::fill::zeros);
    for (int i = 0; i < otu_table.otu_number; ++i) {
        for (int j = 0; j < otu_table.otu_number; ++j) {
            // Get the total possible permutations between the OTU pair
            // TODO: Check if this is producing desired results
            double otu_pair_possible_permutations = possible_permutations(i) * possible_permutations(j);
            // The 'if' code in block below was ported from statmod::permp and statmod::gaussquad
            if (otu_pair_possible_permutations <= 10000 ) {
                // Exact p-value calculation
                // If fewer than 10000 possible permutations, we can safely cast double to int
                pvalues(i, j) = calculateExactPvalue((int)otu_pair_possible_permutations, extreme_value_counts, permutations, i, j);
            } else {
                // Integral approximation for p-value calculation
                pvalues(i, j) = calculatePvalueWithIntegralEstimate(otu_pair_possible_permutations, extreme_value_counts, permutations, i, j);
            }
        }
    }


    // Write out p-values
    writeOutMatrix(pvalues, out_filename, otu_table);
}
