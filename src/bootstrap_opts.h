#ifndef __BOOTSTRAP_OPTS_H__
#define __BOOTSTRAP_OPTS_H__


#include <getopt.h>


#include "common.h"


struct BootstrapOptions {
    std::string otu_filename;
    std::string bootstrap_prefix;

    unsigned int bootstrap_number = 0;
    unsigned int threads = 1;
    unsigned int seed = 1;
};


void print_help();
void print_version();


BootstrapOptions get_commandline_arguments(int argc, char **argv);


#endif
