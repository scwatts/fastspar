#include "catch.hpp"


#include "../src/pvalue.h"


// Framework to provide main implementation
#define CATCH_CONFIG_MAIN


// Integration test
// TODO: break into unittests
TEST_CASE("p-value integration test") {
    // Input data
    std::string otu_fp("data/fake_data.tsv");
    std::string observed_correlation_fp("data/fake_data_cor.tsv");
    std::vector<std::string> bs_correlation_fps = {"data/fake_data_cor_bs_1.tsv", "data/fake_data_cor_bs_2.tsv", "data/fake_data_cor_bs_3.tsv"};

    std::string pseudo_fp("data/fake_data_pseudo_pvalues.tsv");
    std::string exact_fp("data/fake_data_exact_pvalues.tsv");

    // Load data used in calculation
    OtuTable otu_table;
    otu_table.load_otu_file(otu_fp);
    arma::Mat<float> observed_correlation = load_correlation_file(observed_correlation_fp);

    // Calculate pvalues
    arma::Mat<float> test_pseudo_pvalues = calculate_pvalues(otu_table, observed_correlation, bs_correlation_fps, 3, false, 1);
    arma::Mat<float> test_exact_pvalues = calculate_pvalues(otu_table, observed_correlation, bs_correlation_fps, 3, true, 1);

    // Load data for comparison (pseudo p-values from SparCC, exact previously validated)
    arma::Mat<float> pseudo_pvalues = load_correlation_file(pseudo_fp);
    arma::Mat<float> exact_pvalues = load_correlation_file(exact_fp);

    // Tolerating difference of 0.0001 as output types are rounded
    REQUIRE(arma::approx_equal(test_pseudo_pvalues, pseudo_pvalues, "absdiff", 0.001));
    REQUIRE(arma::approx_equal(test_exact_pvalues, exact_pvalues, "absdiff", 0.001));
}


TEST_CASE("factorial division") {
    // Permutations with repetition. Set of n:
    //      element a repeat na times, element b repeated nb times, ...
    //      n = rep(a, na), rep(b, nb), ...
    // permutations:
    //      factorial(n)
    //  -------------------
    //  factorial(na) * factorial(nb) * ...
    //
    //  required is cancelation of factorials during division to avoid overflow
    // Data
    int sample_number = 32;
    std::unordered_map<float, int> count_frequency = {{10, 1}, {0, 20}, {2, 5}, {3, 6}};
    float result = calculate_possbile_otu_permutations(count_frequency, sample_number);
    REQUIRE(1251795504960 == result);
}
