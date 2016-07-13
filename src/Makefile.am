AM_CXXFLAGS=-std=c++11 -Wall -Wextra
AM_FFLAGS=-Wall -Wextra
#LDADD = -larmadillo -lgsl -lgslcblas -lm $(INTLLIBS)

bin_PROGRAMS = sparcpp sparcpp_bootstrap sparcpp_exact_pvalues
sparcpp_SOURCES = sparcpp.cpp
sparcpp_bootstrap_SOURCES = bootstrap.cpp
sparcpp_exact_pvalues_SOURCES = exact_pvalue.cpp gaussq.f
