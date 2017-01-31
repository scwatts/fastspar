#ifndef __EXACT_PVALUE_OPTS_H__
#define __EXACT_PVALUE_OPTS_H__


#include <getopt.h>


#include "common.h"


struct ExactpvalOptions {
    std::string otu_filename;
    std::string correlation_filename;
    std::string bootstrap_prefix;
    std::string out_filename;

    unsigned int permutations = 0;
    unsigned int threads = 1;
};


void print_help();
void print_version();


ExactpvalOptions get_commandline_arguments(int argc, char **argv);


#endif
