AC_INIT([sparcpp], [0.1], [step.wat@gmail.com])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
# Check for c++ and fortran compilers
AC_PROG_CXX
AC_PROG_F77
# Check libraries are installed
# Manually check for HDF5 wrappers in PATH
AC_PATH_PROGS([HDF5], [h5c++ h5cc h5fc h5pcc], ["no"])
if test "$HDF5" = "no"; then
        AC_MSG_ERROR(could not find the library HDF5)
fi
# Use built-in macros to check for other libraries
AC_CHECK_LIB(gomp, main, [LIBS="$LIBS -lgomp"], AC_MSG_ERROR(could not find the library OpenMP))
AC_CHECK_LIB(armadillo, main, [LIBS="$LIBS -larmadillo"], AC_MSG_ERROR(could not find the library Armadillo))
AC_CHECK_LIB(m, cos, [LIBS="$LIBS -lm"], AC_MSG_ERROR(could not find the library LibM))
AC_CHECK_LIB(gslcblas, cblas_dgemm, [LIBS="-lgslcblas $LIBS"], AC_MSG_ERROR(could not find the library GSL (gslcblas)))
AC_CHECK_LIB(gsl, gsl_ran_dirichlet, [LIBS="-lgsl $LIBS"], AC_MSG_ERROR(could not find the library GSL))
# On helix, we need to link against mpi libraries else building fails
AC_CHECK_LIB(mpi, MPI_Init, [LIBS="$LIBS -lmpi"])
AC_CHECK_LIB(mpi_cxx, MPI_Init, [LIBS="$LIBS -lmpi_cxx"])
# Set config header output
AC_CONFIG_HEADERS([config.h])
# Decalre Makefile outputs
AC_CONFIG_FILES([Makefile src/Makefile])
# Disables rebuild rules (fix issue when cloning github repo and Makefile.am, configure.ac, and aclocal.m4 have borked datastamps)
AM_MAINTAINER_MODE([disable])
AC_OUTPUT
