#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <time.h>
#include <vector>


#include "armadillo"
#include "getopt.h"


struct SquareMatrix {
    arma::Mat<double> elements;
    std::vector<std::string> otus;

    // Using initaliser lists to enable calling on Mat/vector constructors at the time of struct construction
    SquareMatrix(std::vector<double> element_vector, std::vector<std::string> otu_vector) : elements(element_vector), otus(otu_vector) {}
};


struct SparseMatrix {
    arma::Col<double> elements;
    std::vector<std::vector<std::string>> otus;
};

SquareMatrix loadSquareMatrix(std::string& filename) {
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


SparseMatrix filterMatrix(SquareMatrix table, arma::Col<arma::uword> filtered_element_indices){
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


void writeOutSparseMatrix(struct SparseMatrix& matrix, std::string out_filename) {
    // Get stream handle
    std::ofstream outfile;
    outfile.open(out_filename);
    outfile << "otu_1\totu2\tvalue" << std::endl;
    // Write out values
    for (unsigned int i = 0; i < matrix.elements.n_elem; ++i) {
        outfile << matrix.otus[i][0] << '\t' << matrix.otus[i][1] << '\t' << matrix.elements[i] << std::endl;
    }

}


float getFloatFromChar(const char* optarg) {
    // Check at most the first 8 characters are numerical
    std::string optstring(optarg);
    std::string string_float = optstring.substr(0, 8);
    for (std::string::iterator it = string_float.begin(); it != string_float.end(); ++it) {
        if (!isdigit(*it) && (*it) != '.') {
            std::cerr << "This doesn't look like a float: " << optarg << std::endl;;
            exit(1);
        }
    }
    return std::atof(string_float.c_str());
}


void printHelp() {
    std::cerr << "Program: Filter correlation and p-value table. Output as sparse matrices" << std::endl;
    std::cerr << "Version: 0.1a" << std::endl;
    std::cerr << "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage:" << std::endl;
    std::cerr << "  reduce --correlation_table <cf> --pvalue_table <pf> --correlation <rt> --pvalue <pt>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -c/--correlation_table <cf>       OTU input table" << std::endl;
    std::cerr << "  -p/--pvalue_table <pf>            Number of bootstrap samples to generate" << std::endl;
    std::cerr << "  -r/--correlation <rt>             Absolute (sign is ignored) correlation threshold (default: 0.1)" << std::endl;
    std::cerr << "  -v/--pvalue <pt>                  P-value threshold (default: 0.05)" << std::endl;

}


int main(int argc, char **argv) {
    // Define some variables
    float correlation_threshold = 0.1;
    float pvalue_threshold = 0.05;
    std::string correlation_filename;
    std::string pvalue_filename;

    // Commandline arguments (for getlongtops)
    struct option long_options[] =
        {
            {"correlation_table", required_argument, NULL, 'c'},
            {"pvalue_table", required_argument, NULL, 'p'},
            {"correlation", required_argument, NULL, 'r'},
            {"pvalue", required_argument, NULL, 'v'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, 0, 0}
        };

    // Check if have an attemp at arguments, else print help
    if (argc < 2) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: options -o/--correlation_table and -p/--pvalue_table are required" << std::endl;
        exit(0);
    }

    // Parse commandline arguments
    while (1) {
        int option_index = 0;
        int c;
        c = getopt_long (argc, argv, "hc:p:r:v:", long_options, &option_index);
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
                correlation_threshold = getFloatFromChar(optarg);
                break;
            case 'v':
                pvalue_threshold = getFloatFromChar(optarg);
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
        std::cerr << ": error: argument -o/--correlation_table is required" << std::endl;
        exit(1);
    }
    if (pvalue_filename.empty()) {
        printHelp();
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: argument -p/--pvalue_table is required" << std::endl;
        exit(1);
    }
    // Check that the OTU file exists
    std::ifstream checkfile;
    checkfile.open(correlation_filename);
    if (!checkfile.good()) {
        std::cerr << std::endl << argv[0];
        std::cerr << ": error: OTU table file " << correlation_filename << " does not exist" << std::endl;
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
    SquareMatrix correlation_table = loadSquareMatrix(correlation_filename);
    SquareMatrix pvalue_table = loadSquareMatrix(pvalue_filename);

    // Find index of elements which both pass specified threshold
    arma::Col<arma::uword> filtered_element_indices = arma::find(arma::abs(correlation_table.elements) >= correlation_threshold &&
                                                          pvalue_table.elements <= pvalue_threshold);
    // Will have a generalised function for this
    SparseMatrix filtered_correlation_table = filterMatrix(correlation_table, filtered_element_indices);
    SparseMatrix filtered_pvalue_table = filterMatrix(pvalue_table, filtered_element_indices);

    // Write filtered OTU correlates to file
    writeOutSparseMatrix(filtered_correlation_table, "cor_sparse");
    writeOutSparseMatrix(filtered_pvalue_table, "pval_sparse");
}
