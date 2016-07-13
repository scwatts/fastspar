FC=gfortran
FFLAGS=-Wall -O2

CXX=g++
CFLAGS=-std=c++11 -Wall -O2
LDLIBS=-larmadillo -lgsl -lgslcblas -lm

SOURCES=sparcpp.cpp bootstrap.cpp exact_pvalue.cpp gaussq.f
# TODO: Replace with expression
OBJECTS=sparcpp.o bootstrap.o exact_pvalue.o gaussq.o
EXECUTABLES=sparcpp bootstrap exact_pvalues

all: $(SOURCES) $(EXECUTABLES)

#$(EXECUTABLE): $(OBJECTS)
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


#g++ -Wall -O2 -std=c++11 -c exact_pvalue.cpp
#gfortran -Wall -O2 -c gaussq.f
#gfortran -Wall -O2 -o exactpvalues exact_pvalue.o gaussq.o /usr/lib/gcc/x86_64-linux-gnu/5/libstdc++.a -lgsl -lgslcblas -lm
