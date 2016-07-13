AC_INIT([sparcpp], [0.1], [step.wat@gmail.com])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CXX
AC_PROG_F77([gfortran])
LIBS="-larmadillo -lgsl -lgslcblas -lm"
# On helix, we need to link against mpi libraries else building fails
AC_CHECK_LIB(mpi, MPI_Init, [LIBS="$LIBS -lmpi"])
AC_CHECK_LIB(mpi_cxx, MPI_Init, [LIBS="$LIBS -lmpi_cxx"])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
