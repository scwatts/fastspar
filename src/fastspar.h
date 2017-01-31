#ifndef __FASTSPAR_H__
#define __FASTSPAR_H__


#include <cmath>
#include <time.h>


#include "common.h"
#if defined(FASTSPAR_CPACKAGE)
    #include "fastspar_opts.h"
#endif


struct FastSpar {
    // Some algorithm parameters
    unsigned int iterations;
    unsigned int exclusion_iterations;
    unsigned int exclusion_threshold;

    // RNG
    gsl_rng *p_rng;

    // Threads to use during iteration
    unsigned int threads;

    // Filenames
    std::string otu_filename;
    std::string correlation_output_filename;
    std::string covariance_output_filename;

    // OTU count matrix
    const OtuTable *otu_table;

    // List of each correlation and covariance matrix calculated during iterations
    std::vector<arma::Mat<double>> correlation_vector;
    std::vector<arma::Col<double>> covariance_vector;

    // The median correlation covariance of all iterations
    arma::Mat<double> median_correlation;
    arma::Mat<double> median_covariance;


    // Construct FastSpar with a given otu_table and other parameters
    FastSpar(const OtuTable *_otu_table, unsigned int _iterations, unsigned int _exclusion_iterations, unsigned int _exclusion_threshold, unsigned int _threads);

    // Infer correlation and covariance for n iterations
    void infer_correlation_and_covariance();

    // Calculate the final correlation and covariance for OTUs
    void calculate_median_correlation_and_covariance();
};


struct FastSparIteration {
    // Some variables
    const OtuTable *otu_table;
    unsigned int exclusion_iterations;
    unsigned int exclusion_threshold;

    // Estimated fractions of OTUs
    arma::Mat<double> fractions;
    // Variance of estimated OTU fractions
    arma::Mat<double> variance;

    // List of highly OTU pairs excluded in this iteration and number of components remaining
    std::vector<arma::uword> excluded;
    unsigned int components_remaining;

    // Basis variance vector
    arma::Col<double> basis_variance;
    // Modifier matrix
    arma::Mat<double> mod;

    // Correlation and covariance for this iteration
    arma::Mat<double> basis_correlation;
    arma::Mat<double> basis_covariance;


    // Construct FastSparIterations with a given otu_table and other parameters
    FastSparIteration(const OtuTable *_otu_table, unsigned int _exclusion_iterations, unsigned int _exclusion_threshold);

    // Estimates component fraction of all OTUs across samples using a dirichlet distribution
    void estimate_component_fractions(gsl_rng *p_rng);

    // Calculate the log ratio variance of the estimated fractions
    void calculate_log_ratio_variance();

    // Calcaulte the component variance
    void calculate_component_variance();

    // Calculate the basis correlation and covariance
    void calculate_correlation_and_covariance();

    // Find and exclude highly correlation OTU pairs
    void find_and_exclude_pairs(float threshold);
};


#endif
