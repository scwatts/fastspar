#ifndef __FASTSPAR_OPTS__
#define __FASTSPAR_OPTS__


#include <getopt.h>


#include "common.h"


struct FastsparOptions {
    // Set some default parameters
    unsigned int iterations = 50;
    unsigned int exclude_iterations = 10;
    float threshold = 0.1;
    unsigned int threads = 1;
    unsigned int seed = 1;
    bool assume_yes = false;

    // Declare some important variables
    std::string otu_filename;
    std::string correlation_filename;
    std::string covariance_filename;
};


void print_help();
void print_version();


FastsparOptions get_commandline_arguments(int argc, char **argv);


#endif
