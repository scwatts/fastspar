# Compiler flags
AM_CXXFLAGS=-std=c++11 -Wall -Wextra -fopenmp
AM_FFLAGS=-Wall -Wextra
# Declaring files
bin_PROGRAMS = fastspar fastspar_bootstrap fastspar_exactpvalues fastspar_reduce
fastspar_SOURCES = fastspar.cpp fastspar_opts.cpp common.cpp
fastspar_bootstrap_SOURCES = bootstrap.cpp bootstrap_opts.cpp common.cpp
fastspar_exactpvalues_SOURCES = exact_pvalue.cpp exact_pvalue_opts.cpp gaussq.f common.cpp
fastspar_reduce_SOURCES = reduce.cpp reduce_opts.cpp common.cpp
