#ifndef __REDUCE_OPTS_H__
#define __REDUCE_OPTS_H__


#include <fstream>
#include <thread>


#include <getopt.h>


#include "common.h"


struct ReduceOptions {
    float correlation_threshold = 0.1;
    float pvalue_threshold = 0.05;

    std::string correlation_filename;
    std::string pvalue_filename;
    std::string output_prefix;
};


void print_help();
void print_version();


ReduceOptions get_commandline_arguments(int argc, char **argv);


#endif
