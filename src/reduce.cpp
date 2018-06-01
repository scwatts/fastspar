#include "reduce.h"


// TODO: Have this as a method? Also I think there is a good amount of code reuse in
// common.cpp load correlation mat
SquareMatrix load_square_matrix(std::string filename) {
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
    std::vector<float> element_vector;
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
            // Add current element to correlation mat after converting to float
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
    arma::Col<float> filtered_table = table.elements(filtered_element_indices);

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


void write_sparse_matrix(struct SparseMatrix &matrix, std::string out_filename) {
    // Get stream handle
    std::ofstream outfile;
    outfile.open(out_filename);
    outfile << "otu_1\totu2\tvalue" << std::endl;
    // Write out values
    for (unsigned int i = 0; i < matrix.elements.n_elem; ++i) {
        outfile << matrix.otus[i][0] << '\t' << matrix.otus[i][1] << '\t' << matrix.elements[i] << std::endl;
    }

}


#if defined(FASTSPAR_CPACKAGE) && !defined(LIBRARY)
int main(int argc, char **argv) {
    // Get commandline arguments
    ReduceOptions reduce_options = get_commandline_arguments(argc, argv);

    // Load correlation and p-value tables
    SquareMatrix correlation_table = load_square_matrix(reduce_options.correlation_filename);
    SquareMatrix pvalue_table = load_square_matrix(reduce_options.pvalue_filename);

    // Find index of elements which both pass specified threshold
    arma::Col<arma::uword> all_filtered_element_indices = arma::find(arma::abs(correlation_table.elements) >= reduce_options.correlation_threshold &&
                                                          pvalue_table.elements <= reduce_options.pvalue_threshold);

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
    write_sparse_matrix(filtered_correlation_table, reduce_options.output_prefix + "_filtered_correlation.tsv");
    write_sparse_matrix(filtered_pvalue_table, reduce_options.output_prefix + "_filtered_pvalue.tsv");
}
#endif
