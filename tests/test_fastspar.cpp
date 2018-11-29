#include "catch.hpp"


#include "../src/fastspar.h"


// Integration test
// TODO: break into unittests
TEST_CASE("Correlation, covariance statisitc integration test") {
    // Input files
    std::string data_dir = std::string(DATADIR);
    std::string otu_fp = data_dir + "data/fake_data.tsv";
    std::string correlation_fp = data_dir + "data/fake_data_cor.tsv";
    std::string covariance_fp = data_dir + "data/fake_data_cov.tsv";

    // Run correlations
    OtuTable otu_table;
    otu_table.load_otu_file(otu_fp);
    FastSpar fastspar(&otu_table, 50, 10, 0.1, 1, 1);
    fastspar.infer_correlation_and_covariance();
    fastspar.calculate_median_correlation_and_covariance();

    // Compare to previously validated outputs
    arma::Mat<float> correlation = load_correlation_file(correlation_fp);
    arma::Mat<float> covariance = load_correlation_file(covariance_fp);

    // Tolerating difference of 0.0001 as output types are rounded
    REQUIRE(arma::approx_equal(fastspar.median_correlation, correlation, "absdiff", 0.001));
    REQUIRE(arma::approx_equal(fastspar.median_covariance, covariance, "absdiff", 0.001));
}
