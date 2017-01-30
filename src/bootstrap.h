#ifndef __BOOTSTRAP_H__
#define __BOOTSTRAP_H__


#include <time.h>


#include <gsl/gsl_randist.h>


#include "common.h"
#if defined(FASTSPAR_CPACKAGE)
    #include "bootstrap_opts.h"
#endif


// Get bootstrap samples
arma::Mat<double> get_bootstrap(OtuTable &otu_table, gsl_rng *p_rng);


#endif
