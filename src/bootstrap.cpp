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


void get_and_write_bootstraps(OtuTable &otu_table, unsigned int bootstrap_number, std::string prefix, unsigned int threads, unsigned int seed) {
    // OpenMP function from omp.h. This sets the number of threads in a more reliable way but also ignores OMP_NUM_THREADS
    omp_set_num_threads(threads);

    // Set up rng environment and seed
    gsl_rng *p_rng = get_default_rng_handle(seed);

    // Get specified number of bootstrap samples
    printf("Generating bootstrapped samples\n");
#pragma omp parallel for schedule(static, 1)
    for (unsigned int i = 0; i < bootstrap_number; ++i) {
        // Get the bootstrap
        arma::Mat<double> bootstrap = get_bootstrap(otu_table, p_rng);

        // Transpose matrix in place before writing out
        printf("\tWriting out bootstrapped %i\n", i);
        arma::inplace_trans(bootstrap);
        std::string bootstrap_filename = prefix + "_"  + std::to_string(i) + ".tsv";
        write_out_bootstrap_table(bootstrap, otu_table.otu_ids, bootstrap_filename);
    }

    // Free rng memory
    gsl_rng_free(p_rng);
}


void write_out_bootstrap_table(arma::Mat<double> &bootstrap, std::vector<std::string> otu_ids, std::string filepath) {
    // Get file handle
     FILE *filehandle = fopen(filepath.c_str(), "w");

    // Write out header
    fprintf(filehandle, "#OTU ID");
    for (unsigned int i = 0; i < bootstrap.n_cols; ++i) {
        fprintf(filehandle, "\tbs_%u", i);
    }
    fprintf(filehandle, "\n");

    // Write out values
    for (unsigned int i = 0; i < bootstrap.n_rows; ++i) {
        for (unsigned int j = 0; j < bootstrap.n_cols; ++j) {
            // Write the OTU id as first field in row
            if (j == 0) {
                fprintf(filehandle, "%s", otu_ids[i].c_str());
            }
            fprintf(filehandle, "\t%0.4f", bootstrap(i, j));
        }
        fprintf(filehandle, "\n");
    }

    fclose(filehandle);
}


#if defined(FASTSPAR_CPACKAGE) && !defined(LIBRARY)
int main(int argc, char **argv) {
    // Get commandline arguments
    BootstrapOptions bootstrap_options = get_commandline_arguments(argc, argv);

    // Load the OTU table from file
    printf("Loading OTU count table\n");
    OtuTable otu_table;
    otu_table.load_otu_file(bootstrap_options.otu_filename);

    // Generate bootstraps
    get_and_write_bootstraps(otu_table, bootstrap_options.bootstrap_number, bootstrap_options.bootstrap_prefix, bootstrap_options.threads, bootstrap_options.seed);

}
#endif
