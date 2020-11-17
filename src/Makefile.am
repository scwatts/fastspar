# Compiler flags
AM_CXXFLAGS = -std=c++11 -Wall -Wextra -fopenmp
AM_FFLAGS = -Wall -Wextra -fPIE

# If we're not using the Armadillo wrapper, set macro as described in documentation
if ARMA_DONT_USE_WRAPPER
AM_CXXFLAGS += -DARMA_DONT_USE_WRAPPER
endif

# Compile static binary
if STATIC_COMPILE
AM_CXXFLAGS += -static
endif

# Sources
common_fps = common.cpp common.h local_config.h
fastspar_fps = fastspar.cpp fastspar.h fastspar_opts.cpp fastspar_opts.h $(common_fps)
fastspar_bootstrap_fps = bootstrap.cpp bootstrap.h bootstrap_opts.cpp bootstrap_opts.h $(common_fps)
fastspar_pvalues_fps = pvalue.cpp pvalue.h pvalue_opts.cpp pvalue_opts.h gaussq.f $(common_fps)
fastspar_reduce_fps = reduce.cpp reduce.h reduce_opts.cpp reduce_opts.h $(common_fps)

# Libraries for testing
check_LIBRARIES = libfastspar.a libfastspar_bootstrap.a libfastspar_pvalues.a libfastspar_reduce.a
libfastspar_a_SOURCES = $(fastspar_fps)
libfastspar_bootstrap_a_SOURCES = $(fastspar_bootstrap_fps)
libfastspar_pvalues_a_SOURCES = $(fastspar_pvalues_fps)
libfastspar_reduce_a_SOURCES = $(fastspar_reduce_fps)

libfastspar_a_CXXFLAGS = -DLIBRARY $(AM_CXXFLAGS)
libfastspar_bootstrap_a_CXXFLAGS = -DLIBRARY $(AM_CXXFLAGS)
libfastspar_pvalues_a_CXXFLAGS = -DLIBRARY $(AM_CXXFLAGS)
libfastspar_reduce_a_CXXFLAGS = -DLIBRARY $(AM_CXXFLAGS)

# ELFs
bin_PROGRAMS = fastspar fastspar_bootstrap fastspar_pvalues fastspar_reduce
fastspar_SOURCES = $(fastspar_fps)
fastspar_bootstrap_SOURCES = $(fastspar_bootstrap_fps)
fastspar_pvalues_SOURCES = $(fastspar_pvalues_fps)
fastspar_reduce_SOURCES = $(fastspar_reduce_fps)
