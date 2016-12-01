#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>


#include <armadillo>
#include <getopt.h>


#include "common.h"
#include "reduce.h"


// TODO: Have this as a method? Also I think there is a good amount of code reuse in
// common.cpp load correlation mat
SquareMatrix load_square_matrix(std::string& filename) {
    // Used to store strings from file prior to matrix construction
    int element_count;
    std::string line;
    std::string ele;
    std::stringstream line_stream;

    /* We first want to get the dimensions of the square matrix so we can allocate memory */
    // Open file stream
    std::ifstream square_file;
    square_file.open(filename);
    // Get the header line and count the number of tabs dimensions
    std::getline(square_file, line);
    element_count = std::count(line.begin(), line.end(), '\t');

    /* Now read in elements of sequare matrix */
    std::vector<double> element_vector;
    std::vector<std::string> otu_vector;
    element_vector.reserve(element_count * element_count);
    otu_vector.reserve(element_count);
    // Other variables
    bool id;
    // Process correlation elements
    while(std::getline(square_file, line)) {
        // (Re)sets variables for loop
        id = true;
        line_stream.clear();
        // Add current line to line stream and then split by tabs
        line_stream.str(line);
        while (std::getline(line_stream, ele, '\t')) {
            // Skip the OTU ID column
            if (id) {
                otu_vector.push_back(ele);
                id = false;
                continue;
            }
            // Add current element to correlation mat after converting to double
            element_vector.push_back(std::stod(ele));
        }
    }
    // Initialise SquareMatrix and return
    SquareMatrix square_matrix(element_vector, otu_vector);
    square_matrix.elements.reshape(element_count, element_count);
    return square_matrix;
}


SparseMatrix filter_matrix(SquareMatrix table, arma::Col<arma::uword> filtered_element_indices) {
    // Extract elements which pass filter
    arma::Col<double> filtered_table = table.elements(filtered_element_indices);

    // Record col and row names for OTUs elements in vector
    std::vector<std::vector<std::string>> otus(filtered_element_indices.n_elem, std::vector<std::string>(2));
    for (unsigned int i = 0; i < filtered_element_indices.n_elem; ++i) {
        // Get the how far along the current row this element is (using modulus)
        otus[i][0] = table.otus[filtered_element_indices(i) % table.otus.size()];
        // Get which row this element is in
        otus[i][1] = table.otus[floor(filtered_element_indices(i) / table.otus.size())];
    }

    // Create a SparseMatrix struct and return it
    SparseMatrix sparse_matrix;
    sparse_matrix.elements = filtered_table;
    sparse_matrix.otus = otus;
    return sparse_matrix;
}


void write_sparse_matrix(struct SparseMatrix& matrix, std::string out_filename) {
    // Get stream handle
    std::ofstream outfile;
    outfile.open(out_filename);
    outfile << "otu_1\totu2\tvalue" << std::endl;
    // Write out values
    for (unsigned int i = 0; i < matrix.elements.n_elem; ++i) {
        outfile << matrix.otus[i][0] << '\t' << matrix.otus[i][1] << '\t' << matrix.elements[i] << std::endl;
    }

}


void printHelp() {
    std::cerr << "Program: Filter correlation and p-value table. Output as sparse matrices" << std::endl;
    std::cerr << "Version 0.0.2" << std::endl;
    std::cerr << "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage:" << std::endl;
    std::cerr << "  reduce --correlation_table <cf> --pvalue_table <pf> --correlation <rt> --pvalue <pt> --output_prefix <op>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -c/--correlation_table <cf>       Correlation input table" << std::endl;
    std::cerr << "  -p/--pvalue_table <pf>            Number of bootstrap samples to generate" << std::endl;
    std::cerr << "  -r/--correlation <rt>             Absolute (sign is ignored) correlation threshold (default: 0.1)" << std::endl;
    std::cerr << "  -v/--pvalue <pt>                  P-value threshold (default: 0.05)" << std::endl;
    std::cerr << "  -o/--output_prefix <op>           Output prefix" << std::endl;

}


int main(int argc, char **argv) {
    // Define some variables
    float correlation_threshold = 0.1;
    float pvalue_threshold = 0.05;
    std::string correlation_filename;
    std::string pvalue_filename;
    std::string output_prefix;

    // Commandline arguments (for getlongtops)
    struct option long_options[] =
        {
            {"correlation_table", required_argument, NULL, 'c'},
            {"pvalue_table", required_argument, NULL, 'p'},
            {"correlation", required_argument, NULL, 'r'},
            {"pvalue", required_argument, NULL, 'v'},
            {"output_prefix", required_argument, NULL, 'o'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, 0, 0}
        };

    // Check if have an attemp at arguments, else print help
    // TODO: Fix this. Will fail to catch if given optional args but not all required...
    if (argc < 3) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: options -c/--correlation_table, -p/--pvalue_table and -o/--output_prefix are required" << std::endl;
        exit(0);
    }

    // Parse commandline arguments
    while (1) {
        int option_index = 0;
        int c;
        c = getopt_long (argc, argv, "hc:p:r:v:o:", long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch(c) {
            // TODO: do we need case(0)?
            case 'c':
                correlation_filename = optarg;
                break;
            case 'p':
                pvalue_filename = optarg;
                break;
            case 'r':
                correlation_threshold = get_float_from_char(optarg);
                break;
            case 'v':
                pvalue_threshold = get_float_from_char(optarg);
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'h':
                printHelp();
                exit(0);
            default:
                printHelp();
                exit(1);
        }
    }

    // Abort execution if given unknown arguments
    if (optind < argc){
        std::cerr << argv[0] << " invalid argument: " << argv[optind++] << std::endl;
    }
    // Make sure we have filenames and parameters
    if (correlation_filename.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -c/--correlation_table is required" << std::endl;
        exit(1);
    }
    if (pvalue_filename.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -p/--pvalue_table is required" << std::endl;
        exit(1);
    }
    if (output_prefix.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -o/--output_prefix is required" << std::endl;
        exit(1);
    }
    // Check that the correlation file exists
    std::ifstream checkfile;
    checkfile.open(correlation_filename);
    if (!checkfile.good()) {
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: correlation table file " << correlation_filename << " does not exist" << std::endl;
        exit(1);
    }
    checkfile.close();
    // Check that the pvalue file exists
    checkfile.open(pvalue_filename);
    if (!checkfile.good()) {
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: p-value table file " << pvalue_filename << " does not exist" << std::endl;
        exit(1);
    }

    // Load correlation and p-value tables
    SquareMatrix correlation_table = load_square_matrix(correlation_filename);
    SquareMatrix pvalue_table = load_square_matrix(pvalue_filename);

    // Find index of elements which both pass specified threshold
    arma::Col<arma::uword> all_filtered_element_indices = arma::find(arma::abs(correlation_table.elements) >= correlation_threshold &&
                                                          pvalue_table.elements <= pvalue_threshold);

    // Select indices in the upper triangle of the matrix (to prevent outputing the same
    // correlation or p-value twice)
    std::vector<arma::uword> filtered_element_indices_vector;
    for (unsigned int i = 0; i < all_filtered_element_indices.n_elem; ++i) {
        // Convert linear index (i) to subscript (i,j)
        arma::uword linear_index = all_filtered_element_indices(i);
        arma::Col<arma::uword> subscript_index = arma::ind2sub(arma::size(correlation_table.elements), linear_index);

        // Check it is the in the upper half and then add it to the filtered vector
        if (subscript_index(0) >= subscript_index(1)) {
            filtered_element_indices_vector.push_back(linear_index);
        }
    }

    // Convert the vector of filtered indices so we can use it with Armadillo
    arma::Col<arma::uword> filtered_element_indices(filtered_element_indices_vector);

    // Will have a generalised function for this
    SparseMatrix filtered_correlation_table = filter_matrix(correlation_table, filtered_element_indices);
    SparseMatrix filtered_pvalue_table = filter_matrix(pvalue_table, filtered_element_indices);

    // Write filtered OTU correlates to file
    write_sparse_matrix(filtered_correlation_table, output_prefix + "_filtered_correlation.tsv");
    write_sparse_matrix(filtered_pvalue_table, output_prefix + "_filtered_pvalue.tsv");
}
