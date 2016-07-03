#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "armadillo"

struct OtuTable {
    std::vector<std::string> samples;
    std::vector<std::string> otus;
    std::vector<std::vector<int>> counts;
};

void LoadOtuFile(struct OtuTable * potu_table, std::string filename) {
    // Used to store strings from file prior to assignment
    std::string line;
    std::string ele;
    std::stringstream line_stream;
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
        potu_table->samples.push_back(ele);
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
                potu_table->otus.push_back(ele);
                id = false;
                continue;
            }
            // Add current element to OTU count after converting to int
            otu_counts.push_back(std::stoi(ele));
        }
        // Add current OTU counts to OtuTable struct instance
        potu_table->counts.push_back(otu_counts);
    // TODO: Check if growing std::vector is sustainable for large tables
    }
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
    std::vector<int> counts;
    for(count_container::iterator row = otu_table.counts.begin(); row != otu_table.counts.end(); ++row) {
        for(std::vector<int>::iterator ele = (*row).begin(); ele != row->end(); ++ele){
            counts.push_back(*ele);
        }
    }
    // Construct armadillo integer matrix, column major construction
    arma::Mat<int> B(counts);
    B.reshape(200, 50);
    arma::Mat<int> Bt = B.t();
}