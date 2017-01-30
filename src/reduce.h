#ifndef __REDUCE_H__
#define __REDUCE_H__


#include <iostream>
#include <sstream>


#include "reduce_opts.h"
#include "common.h"


// Structure to hold information about a square matrix
struct SquareMatrix {
    arma::Mat<double> elements;
    std::vector<std::string> otus;

    // Using initaliser lists to enable calling on Mat/vector constructors at the time of struct construction
    SquareMatrix(std::vector<double> element_vector, std::vector<std::string> otu_vector) : elements(element_vector), otus(otu_vector) {}
};


// Structure to hold information about a sparse matrix
struct SparseMatrix {
    arma::Col<double> elements;
    std::vector<std::vector<std::string>> otus;
};


// Load a square matrix (e.g. correlation matrix) and return a SquareMatrux object
SquareMatrix load_square_matrix(std::string filename);

// Given two criteria (p-value and correlation threshold), filter a matrix and return a SparseMatrix
SparseMatrix filter_matrix(SquareMatrix table, arma::Col<arma::uword> filtered_element_indices);

// Write a SparseMatrix to file
void write_sparse_matrix(struct SparseMatrix &matrix, std::string out_filename);


#endif
