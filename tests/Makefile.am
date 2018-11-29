# Build flags
AM_CXXFLAGS = -std=c++11 -Wall -Wextra -fopenmp
AM_CPPFLAGS = -DDATADIR=\"$(srcdir)/\"
AM_FFLAGS = -Wall -Wextra

EXTRA_DIST = data

LDADD = ../src/libfastspar.a ../src/libfastspar_bootstrap.a ../src/libfastspar_pvalues.a ../src/libfastspar_reduce.a

check_LIBRARIES = libtest_main.a
libtest_main_a_SOURCES = test_main.cpp catch.hpp

LDADD += libtest_main.a

check_PROGRAMS = test_fastspar test_pvalues
AM_DEFAULT_SOURCE_EXT = .cpp
TESTS = $(check_PROGRAMS)
