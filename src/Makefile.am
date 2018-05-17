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

# Sources
fastspar_fps = fastspar.cpp fastspar_opts.cpp common.cpp
fastspar_bootstrap_fps = bootstrap.cpp bootstrap_opts.cpp common.cpp
fastspar_exactpvalues_fps = exact_pvalue.cpp exact_pvalue_opts.cpp gaussq.f common.cpp
fastspar_reduce_fps = reduce.cpp reduce_opts.cpp common.cpp

# Libraries for testing
check_LIBRARIES = libfastspar.a libfastspar_bootstrap.a libfastspar_exactpvalues.a libfastspar_reduce.a
libfastspar_a_SOURCES = $(fastspar_fps)
libfastspar_bootstrap_a_SOURCES = $(fastspar_bootstrap_fps)
libfastspar_exactpvalues_a_SOURCES = $(fastspar_exactpvalues_fps)
libfastspar_reduce_a_SOURCES = $(fastspar_reduce_fps)

libfastspar_a_CXXFLAGS = -DLIBRARY $(AM_CXXFLAGS)
libfastspar_bootstrap_a_CXXFLAGS = -DLIBRARY $(AM_CXXFLAGS)
libfastspar_exactpvalues_a_CXXFLAGS = -DLIBRARY $(AM_CXXFLAGS)
libfastspar_reduce_a_CXXFLAGS = -DLIBRARY $(AM_CXXFLAGS)

# ELFs
bin_PROGRAMS = fastspar fastspar_bootstrap fastspar_exactpvalues fastspar_reduce
fastspar_SOURCES = $(fastspar_fps)
fastspar_bootstrap_SOURCES = $(fastspar_bootstrap_fps)
fastspar_exactpvalues_SOURCES = $(fastspar_exactpvalues_fps)
fastspar_reduce_SOURCES = $(fastspar_reduce_fps)
