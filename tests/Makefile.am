# Build flags
AM_CXXFLAGS=-std=c++11 -Wall -Wextra -fopenmp
AM_FFLAGS=-Wall -Wextra

check_PROGRAMS = run_tests
TESTS = $(check_PROGRAMS)
run_tests_SOURCES = run_tests.cpp test_fastspar.cpp
run_tests_LDADD = ../src/libfastspar.a ../src/libfastspar_bootstrap.a ../src/libfastspar_exactpvalues.a ../src/libfastspar_reduce.a
