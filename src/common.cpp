#include <iomanip>
#include <string>
#include <vector>


#include "armadillo"
#include "common.h"


// Load an OTU table from file
void OtuTable::load_otu_file(std::string filename) {
    // TODO: Catch now filename or non-existent file
    // Used to store strings from file prior to assignment
    std::string line;
    std::string ele;
    std::stringstream line_stream;
    std::vector<double> temp_counts_vector;
    bool id;
    // Open file stream
    std::ifstream otu_file;
    otu_file.open(filename);
    // Process header
    std::getline(otu_file, line);
    line_stream.str(line);
    // Iterate header columns
    while(std::getline(line_stream, ele, '\t')) {
        //TODO: Add assertion here?
        // Skip the #OTU ID column (first column)
        if (ele == "#OTU ID") {
            continue;
        }
        // Store samples
        sample_names.push_back(ele);
        ++sample_number;
    }
    // Process sample counts, need to get OTU IDS first
    while(std::getline(otu_file, line)) {
        // (Re)sets variables for loop
        id = true;
        line_stream.clear();
        // Add current line to line stream and then split by tabs
        line_stream.str(line);
        while (std::getline(line_stream, ele, '\t')) {
            // Grab the #OTU ID
            if (id) {
                otu_ids.push_back(ele);
                id = false;
                continue;
            }
            // Add current element to OTU count after converting to double; some OTUs may be corrected and therefore a double
            temp_counts_vector.push_back(std::stod(ele));
        }
        ++otu_number;
    }
    // Finally construct the OTU observation matrix and _move_ to struct
    arma::Mat<double> temp_otu_matrix(temp_counts_vector);
    temp_otu_matrix.reshape(sample_number, otu_number);
    counts = std::move(temp_otu_matrix);
}


// Load a correlation table from file
arma::Mat<double> load_correlation_file(std::string& filename) {
    // Used to store strings from file prior to matrix construction and other variables
    std::string line;
    std::string ele;
    std::stringstream line_stream;
    std::vector<double> correlations_vector;
    bool id;

    // Open file stream
    std::ifstream cor_file;
    cor_file.open(filename);

    // Skip header, order SHOULD be the same as input OTU table
    std::getline(cor_file, line);
    line_stream.str(line);

    // Count OTUs in header if needed
    int otu_number = std::count(line.begin(), line.end(), '\t');
    // Reserving space memory
    correlations_vector.reserve(otu_number * otu_number);

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
    correlations.reshape(otu_number, otu_number);
    return correlations;
}


void write_out_square_otu_matrix(arma::Mat<double> & matrix, OtuTable & otu_table, std::string filename) {
    // Get stream handle
    std::ofstream outfile;
    outfile.open(filename);
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
            outfile << std::fixed << std::setw(8) << std::setprecision(4) << "\t" <<
            matrix(i, j);
        }
        outfile << std::endl;
    }
}


int get_int_from_char(const char* optarg) {
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


float get_float_from_char(const char* optarg) {
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
