#include <fstream>
#include <iostream>
#include <sstream>
#include <time.h>
#include <vector>

#include <armadillo>
#include <gsl/gsl_randist.h>

struct OtuTable {
    std::vector<std::string> sample_names;
    std::vector<std::string> otu_ids;
    std::vector<std::vector<int>> otu_observations;
    int otu_number;
    int sample_number;
};

void LoadOtuFile(struct OtuTable * potu_table, std::string filename) {
    // Used to store strings from file prior to assignment
    std::string line;
    std::string ele;
    std::stringstream line_stream;
    int otu_number = 0;
    int sample_number = 0;
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
        potu_table->sample_names.push_back(ele);
        ++sample_number;
    }
    // Process sample counts, need to get OTU IDS first
    while(std::getline(otu_file, line)) {
        // TODO: Is there an alternate design pattern to loop variables as below
        // Loop variables
        std::vector<int> otu_counts;
        bool id;
        // (Re)sets variables for loop
        id = true;
        otu_counts.clear();
        line_stream.clear();
        // Add current line to line stream and then split by tabs
        line_stream.str(line);
        while (std::getline(line_stream, ele, '\t')) {
            // Grab the OTU_id
            if (id) {
                potu_table->otu_ids.push_back(ele);
                id = false;
                continue;
            }
            // Add current element to OTU count after converting to int
            otu_counts.push_back(std::stoi(ele));
        }
        // Add current OTU counts to OtuTable struct instance
        potu_table->otu_observations.push_back(otu_counts);
        ++otu_number;
    // TODO: Check if growing std::vector is sustainable for large tables
    }
    potu_table->otu_number = otu_number;
    potu_table->sample_number = sample_number;
}

int main() {
    // Load the OTU file from
    std::string otu_filename;
    otu_filename = "fake_data.txt";
    struct OtuTable otu_table;
    // TODO: Check if this is the most appropriate design (init struct and then passing to function)
    LoadOtuFile(&otu_table, otu_filename);
    // TODO: It appears from the armadillo docs you can initialise an arma:mat from a 1-d vector.
    //       This means that when loading the OTU table, I could record the dimensions and have a
    //       1-d std::vector<int>. For now I'll leave as is given that datastructure may be useful later
    // Flatten std::vector of otu counts std::vector<double>
    // TODO: Find the most suitable place for typedefs
    typedef std::vector<std::vector<int>> count_container;
    count_container * otu_observations = &otu_table.otu_observations;
    std::vector<int> counts_vector;
    for(count_container::iterator row = otu_observations->begin(); row != otu_observations->end(); ++row) {
        for(std::vector<int>::iterator ele = (*row).begin(); ele != row->end(); ++ele){
            // TODO: use a vector of int pointers rather than copying
            counts_vector.push_back(*ele);
        }
    }
    // Construct armadillo integer matrix, column major construction
    arma::Mat<int> counts(counts_vector);
    counts.reshape(otu_table.sample_number, otu_table.otu_number);

    // STEP 1: Estimate component fractions
    // TODO: Check if fractions matrix is used elsewhere (I don't think it is)
    double fraction_fill = 1 / (double)otu_table.otu_number;
    arma::Mat<double> fractions(otu_table.sample_number, otu_table.otu_number);
    fractions.fill(fraction_fill);
    // Get arma row and add pseudo count (then convert to double vector for rng function)
    int i = 0;
    arma::Row<int> row_pseudocount = counts.row(i) + 1;
    std::vector<double> row_pseudocount_vector = arma::conv_to<std::vector<double>>::from(row_pseudocount);

    // Set up rng environment and seed
    const gsl_rng_type * rng_type;
    gsl_rng_env_setup();
    // gsl_rng_default is a global
    rng_type = gsl_rng_default;
    gsl_rng * p_rng = gsl_rng_alloc(rng_type);
    gsl_rng_set(p_rng, time(NULL));

    // Draw from dirichlet dist, storing results in theta double array
    size_t row_size = row_pseudocount_vector.size();
    double theta[row_size];
    // The function takes double arrays and it seems that you must pass the address of the first element to function
    gsl_ran_dirichlet(p_rng, row_size, &row_pseudocount_vector[0], &theta[0]);
    // Free the rng
    gsl_rng_free(p_rng);

    // TODO: create list of fraction vectors and init arma Mat as above
    // Create arma::Row from double[] and update fractions row
    // Update fractions row
    arma::Mat<double> estimated_fractions_row(theta, 1, 50);
    fractions.row(i) = estimated_fractions_row;

    /*double sum = 0;
    for(int i = 0; i < row_size; ++i){
        sum += theta[i];
        std::cout << theta[i] << " " << row_pseudocount_vector[i] << std::endl;
    }
    std::cout << sum << std::endl;
    std::cout << row_size << std::endl;*/


    //std::vector<int> row;
    //row = counts.col(i)
    //row.for_each() += 1;
    //fractions.col(i) = dirichlet(row);
}