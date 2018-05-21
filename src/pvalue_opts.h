#ifndef __PVALUE_OPTS_H__
#define __PVALUE_OPTS_H__


#include <getopt.h>


#include "common.h"


struct PvalOptions {
    std::string otu_filename;
    std::string correlation_filename;
    std::string bootstrap_prefix;
    std::string out_filename;
    bool exact = true;

    unsigned int permutations = 0;
    unsigned int threads = 1;
};


void print_help();
void print_version();


PvalOptions get_commandline_arguments(int argc, char **argv);


#endif
