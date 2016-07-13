# Fortran compiler and flags
FC=gfortran
FFLAGS=-Wall -O2

# C++ compiler and flags
CXX=g++
CFLAGS=-std=c++11 -Wall -O2

# Common shared libraries
LDLIBS=-larmadillo -lgsl -lgslcblas -lm

# In and out files
SOURCES=sparcpp.cpp bootstrap.cpp exact_pvalue.cpp gaussq.f
# TODO: Replace with expression
OBJECTS=sparcpp.o bootstrap.o exact_pvalue.o gaussq.o
EXECUTABLES=sparcpp bootstrap exact_pvalues

all: $(SOURCES) $(EXECUTABLES)

# TODO: Generalise linking of executables
sparcpp: sparcpp.o
	$(CXX) $(LDLIBS) $< -o $@

bootstrap: bootstrap.o
	$(CXX) $(LDLIBS) $< -o $@

exact_pvalues: exact_pvalue.o gaussq.o
	$(FC) $(LDLIBS) -lstdc++ $? -o $@

.cpp.o:
	$(CXX) -c $(CFLAGS) $< -o $@

.f.o:
	$(FC) -c $(FFLAGS) $< -o $@

.PHONY: clean
clean:
	rm $(OBJECTS) $(EXECUTABLES)
