#include "bootstrap.h"


// Get bootstrap samples
arma::Mat<double> get_bootstrap(OtuTable &otu_table, gsl_rng *p_rng) {
    arma::Mat<double> bootstrap(otu_table.sample_number, otu_table.otu_number);
    for (int j = 0; j < otu_table.otu_number; ++j) {
        // Get a random list of sample indices equal to the number of samples
        std::vector<arma::uword> indices(otu_table.sample_number);
        for (std::vector<arma::uword>::iterator it = indices.begin(); it != indices.end(); ++it){
            // With the random row element, determine position in counts matrix
            int r = gsl_rng_uniform_int(p_rng, otu_table.sample_number);
            *it = arma::sub2ind(arma::size(otu_table.counts), r, j);
        }
        // Need to convert to uvec (Col<uword>) for use with arma subviews
        arma::Col<arma::uword> arma_indices(indices);
        // Gather random elements and set row in bootstrap matrix
        bootstrap.col(j) = otu_table.counts(arma_indices);
    }
    return bootstrap;
}


#if defined(FASTSPAR_CPACKAGE)
int main(int argc, char **argv) {
    // Get commandline arguments
    BootstrapOptions bootstrap_options = get_commandline_arguments(argc, argv);

    // OpenMP function from omp.h. This sets the number of threads in a more reliable way but also ignores OMP_NUM_THREADS
    omp_set_num_threads(bootstrap_options.threads);

    // Set up rng environment and seed
    const gsl_rng_type * rng_type;
    gsl_rng_env_setup();
    // gsl_rng_default is a global
    rng_type = gsl_rng_default;
    gsl_rng * p_rng = gsl_rng_alloc(rng_type);
    gsl_rng_set(p_rng, time(NULL));

    // Load the OTU table from file
    printf("Loading OTU count table\n");
    OtuTable otu_table;
    otu_table.load_otu_file(bootstrap_options.otu_filename);

    // Get specified number of bootstrap samples
    printf("Generating bootstrapped samples\n");
#pragma omp parallel for schedule(static, 1)
    for (unsigned int i = 0; i < bootstrap_options.bootstrap_number; ++i) {
        // Get the bootstrap
        arma::Mat<double> bootstrap = get_bootstrap(otu_table, p_rng);

        // Transpose matrix in place before writing out
        printf("\tWriting out bootstrapped %i\n", i);
        arma::inplace_trans(bootstrap);
        std::string bootstrap_filename = bootstrap_options.bootstrap_prefix + "_"  + std::to_string(i) + ".tsv";
        write_out_square_otu_matrix(bootstrap, otu_table, bootstrap_filename);
    }

    // Free rng memory
    gsl_rng_free(p_rng);
}
#endif
