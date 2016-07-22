#include <time.h>
#include <vector>


#include <armadillo>
#include <getopt.h>
#include <gsl/gsl_randist.h>


#include "common.h"


// Get bootstrap samples
arma::Mat<double> get_bootstrap(OtuTable & otu_table, gsl_rng * p_rng) {
        arma::Mat<double> bootstrap(otu_table.sample_number, otu_table.otu_number);
        for (int j = 0; j < otu_table.otu_number; ++j) {
            // Get a random list of sample indices equal to the number of samples
            std::vector<arma::uword> indices(otu_table.sample_number);
            for (std::vector<arma::uword>::iterator it = indices.begin(); it != indices.end(); ++it){
                // With the random row element, determine position in counts matrix
                int r = gsl_rng_uniform_int(p_rng, otu_table.sample_number);
                *it = arma::sub2ind(arma::size(otu_table.counts), r, j);
            }
            // Need to convert to uvec (Col<uword>) for use with arma subviews
            arma::Col<arma::uword> arma_indices(indices);
            // Gather random elements and set row in bootstrap matrix
            bootstrap.col(j) = otu_table.counts(arma_indices);
        }
        return bootstrap;
}


void printHelp() {
    std::cerr << "Program: SparCpp bootstraps (c++ implementation of SparCC)" << std::endl;
    std::cerr << "Version: 0.1a" << std::endl;
    std::cerr << "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage:" << std::endl;
    std::cerr << "  bootstraps --otu_table <of> --number <n> --prefix <p>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -t/--otu_table <of>   OTU input table" << std::endl;
    std::cerr << "  -n/--number <n>       Number of bootstrap samples to generate" << std::endl;
    std::cerr << "  -p/--prefix <p>       Prefix out bootstrap output files" << std::endl;

}


int main(int argc, char **argv) {
    // Define some variables
    int bootstrap_number = 0;
    std::string otu_filename;
    std::string bootstrap_prefix;

    // Commandline arguments (for getlongtops)
    struct option long_options[] =
        {
            {"otu_table", required_argument, NULL, 't'},
            {"number", required_argument, NULL, 'n'},
            {"prefix", required_argument, NULL, 'p'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, 0, 0}
        };

    // Check if have an attemp at arguments, else print help
    if (argc < 2) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: options -t/--otu_table, -n/--number, and -p/--prefix are required" << std::endl;
        exit(0);
    }

    // Parse commandline arguments
    while (1) {
        int option_index = 0;
        int c;
        c = getopt_long (argc, argv, "ht:p:n:", long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch(c) {
            // TODO: do we need case(0)?
            case 't':
                otu_filename = optarg;
                break;
            case 'p':
                bootstrap_prefix = optarg;
                break;
            case 'n':
                bootstrap_number = get_int_from_char(optarg);
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
    if (bootstrap_number == 0) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -n/--number is required" << std::endl;
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

    // Load the OTU table from file
    struct OtuTable otu_table;
    otu_table.load_otu_file(otu_filename);

    // Get specified number of bootstrap samples
    for (int i = 0; i < bootstrap_number; ++i) {
        // Get the bootstrap
        arma::Mat<double> bootstrap = get_bootstrap(otu_table, p_rng);

        // Transpose matrix in place before writing out
        arma::inplace_trans(bootstrap);
        std::string bootstrap_filename = bootstrap_prefix + "_"  + std::to_string(i) + ".tsv";
        write_out_square_otu_matrix(bootstrap, otu_table, bootstrap_filename);
    }

    // Free rng memory
    gsl_rng_free(p_rng);
}
