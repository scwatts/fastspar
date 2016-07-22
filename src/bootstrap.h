#ifndef __BOOTSTRAP_H__
#define __BOOTSTRAP_H__


#include <gsl/gsl_randist.h>
#include <armadillo>


#include "common.h"


// Get bootstrap samples
arma::Mat<double> get_bootstrap(OtuTable & otu_table, gsl_rng * p_rng);


#endif
