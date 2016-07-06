#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <time.h>
#include <vector>

#include <armadillo>
#include <getopt.h>
#include <gsl/gsl_randist.h>


struct OtuTable {
    std::vector<std::string> sample_names;
    std::vector<std::string> otu_ids;
    std::vector<int> otu_observations;
    int otu_number;
    int sample_number;
};


struct OtuTable loadOtuFile(std::string filename) {
    // TODO: Catch now filename or non-existent file
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


void writeOutMatrix(arma::Mat<int>& matrix, std::string out_filename, struct OtuTable& otu_table) {
    // Get stream handle
    std::ofstream outfile;
    outfile.open(out_filename);
    // Write out header
    outfile << "OTU_id";
    for (std::vector<std::string>::iterator it = otu_table.otu_ids.begin(); it != otu_table.otu_ids.end(); ++it) {
        outfile << "\t" << *it;
    }
    outfile << std::endl;
    // Write out values
    for (unsigned int i = 0; i < matrix.n_rows; ++i) {
        for (unsigned int j = 0; j < matrix.n_cols; ++j) {
            // Write the OTU id as first field in row
            if (j == 0) {
                outfile << otu_table.sample_names[i];
            }
            outfile << "\t" << matrix(i, j);
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
        c = getopt_long (argc, argv, "t:p:n:", long_options, &option_index);
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
                bootstrap_number = getIntFromChar(optarg);
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
    struct OtuTable otu_table = loadOtuFile(otu_filename);
    // Construct count matrix
    arma::Mat<int> counts(otu_table.otu_observations);
    counts.reshape(otu_table.sample_number, otu_table.otu_number);

    // Get specified number of bootstrap samples
    for (int i = 0; i < bootstrap_number; ++i) {
        // For each row in the matrix, random select elements (with replacement)
        arma::Mat<int> bootstrap(otu_table.sample_number, otu_table.otu_number);
        std::string bootstrap_filename = bootstrap_prefix + "_"  + std::to_string(i) + ".tsv";
        for (int j = 0; j < otu_table.otu_number; ++j) {
            // Get a random list of sample indices equal to the number of samples
            std::vector<arma::uword> indices(otu_table.sample_number);
            for (std::vector<arma::uword>::iterator it = indices.begin(); it != indices.end(); ++it){
                // With the random row element, determine position in counts matrix
                int r = gsl_rng_uniform_int(p_rng, otu_table.sample_number);
                *it = arma::sub2ind(arma::size(counts), r, j);
            }
            // Need to convert to uvec (Col<uword>) for use with arma subviews
            arma::Col<arma::uword> arma_indices(indices);
            // Gather random elements and set row in bootstrap matrix
            bootstrap.col(j) = counts(arma_indices);
        }
        // Transpose matrix in place before writing out
        arma::inplace_trans(bootstrap);
        writeOutMatrix(bootstrap, bootstrap_filename, otu_table);
    }

    // Free rng memory
    gsl_rng_free(p_rng);
}
