# Compiler flags
AM_CXXFLAGS=-std=c++11 -Wall -Wextra -fopenmp
AM_FFLAGS=-Wall -Wextra


# If we're not using the Armadillo wrapper, set macro as described in documentation
if ARMA_DONT_USE_WRAPPER
AM_CXXFLAGS += -DARMA_DONT_USE_WRAPPER
endif

# Compile static binary
if STATIC_COMPILE
AM_CXXFLAGS += -static
endif


# Declaring files
bin_PROGRAMS = fastspar fastspar_bootstrap fastspar_exactpvalues fastspar_reduce
fastspar_SOURCES = fastspar.cpp fastspar_opts.cpp common.cpp
fastspar_bootstrap_SOURCES = bootstrap.cpp bootstrap_opts.cpp common.cpp
fastspar_exactpvalues_SOURCES = exact_pvalue.cpp exact_pvalue_opts.cpp gaussq.f common.cpp
fastspar_reduce_SOURCES = reduce.cpp reduce_opts.cpp common.cpp
