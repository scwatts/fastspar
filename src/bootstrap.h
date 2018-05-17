#ifndef __BOOTSTRAP_H__
#define __BOOTSTRAP_H__


#include <time.h>


#include "common.h"
#if defined(FASTSPAR_CPACKAGE)
    #include "bootstrap_opts.h"
#endif


// Get single bootstrap for OTU table
arma::Mat<double> get_bootstrap(OtuTable &otu_table, gsl_rng *p_rng);

// Get n bootstraps for OTU table
void get_and_write_bootstraps(OtuTable &otu_table, unsigned int bootstrap_number, std::string prefix, unsigned int threads, unsigned int seed);

// Write out a bootstrap count table
void write_out_bootstrap_table(arma::Mat<double> &bootstrap, std::vector<std::string> otu_ids, std::string filepath);


#endif
