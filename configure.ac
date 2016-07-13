AC_INIT([sparcpp], [0.1], [step.wat@gmail.com])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
# Check for c++ and fortran compilers
AC_PROG_CXX
AC_PROG_F77([gfortran])
# Check libraries are installed
AC_CHECK_LIB(armadillo, main, [LIBS="$LIBS -larmadillo"], AC_MSG_ERROR(could not find the library Armadillo))
AC_CHECK_LIB(m, cos, [LIBS="$LIBS -lm"], AC_MSG_ERROR(could not find the library LibM))
AC_CHECK_LIB(gslcblas, cblas_dgemm, [LIBS="$LIBS -lgslcblas"], AC_MSG_ERROR(could not find the library GSL (gslcblas)))
AC_CHECK_LIB(gsl, gsl_ran_dirichlet, [LIBS="$LIBS -lgsl"], AC_MSG_ERROR(could not find the library GSL))
# On helix, we need to link against mpi libraries else building fails
AC_CHECK_LIB(mpi, MPI_Init, [LIBS="$LIBS -lmpi"])
AC_CHECK_LIB(mpi_cxx, MPI_Init, [LIBS="$LIBS -lmpi_cxx"])
# Set config header output
AC_CONFIG_HEADERS([config.h])
# Decalre Makefile outputs
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
