#ifndef __EXACT_PVALUE_H__
#define __EXACT_PVALUE_H__


#include <unordered_map>
#include <string>
#include <vector>


#include <armadillo>


// Collect bootstrap correlation paths by globbing a prefix
std::vector<std::string> get_bootstrap_correlation_paths(std::string & glob_path);

// Count bootstrap correlations more extreme than the observed correlation
void count_values_more_extreme(arma::Mat<double>& abs_observed_correlation, arma::Mat<double>& abs_bootstrap_correlation, arma::Mat<int>& extreme_value_counts);

// Calculate the permutations of count data for an OTU
double calculate_possbile_otu_permutations(std::unordered_map<double, int>& count_frequency, int sample_number);

// Calculate the exact p-value (when total nperm <= 10000)
double calculate_exact_pvalue(double otu_pair_possible_permutations, int& extreme_value_count, unsigned int& permutations);

// Calculate the p-value with an integral estimate (when total nperm > 10000)
double calculate_pvalue_with_integral_estimate(double& otu_pair_possible_permutations, int& extreme_value_count, int& permutations);


#endif
