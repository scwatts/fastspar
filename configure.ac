# Initialisation
AC_INIT([fastspar], [0.0.10], [s.watts2@student.unimelb.edu.au])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_LANG([C++])

# Add option for static compilation and option to not use Arma wrapper
AC_ARG_ENABLE([static-compile], AS_HELP_STRING([--enable-static-compile], [Compile static binary; forces --disable-arma-wrapper]),
	[case "${enableval}" in
	  no)  static_compile=false ;;
	  *)   static_compile=true ;;
	esac], [static_compile=false])

AC_ARG_ENABLE([arma-wrapper], AS_HELP_STRING([--disable-arma-wrapper], [Disables use of the Armadillo wrapper]),
	[case "${enableval}" in
	  *)   arma_wrapper=false;;
	esac], [arma_wrapper=true])

# Set conditionals to use in the Makefile.am
AM_CONDITIONAL([ARMA_DONT_USE_WRAPPER],
               [test "x$static_compile" = "xtrue" || test "x$arma_wrapper" = "xfalse"])
AM_CONDITIONAL([STATIC_COMPILE],
               [test "x$static_compile" = "xtrue"])

# Check required compilers and utilities
AC_PROG_CXX
AC_PROG_F77
AC_PROG_RANLIB
AM_PROG_AR

# Check dependencies
if $CXX --version 2>&1 | grep Open64 > /dev/null ; then
  # Open64
  AC_CHECK_LIB([openmp], [omp_get_num_threads],
               [OMP_LIBS="-lopenmp"], [])
elif $CXX --version 2>&1 | grep clang > /dev/null ; then
  # Clang
  AC_CHECK_LIB([omp], [GOMP_parallel_start],
               [OMP_LIBS="-lomp"], [])
else
  # GCC
  AC_CHECK_LIB([gomp], [omp_get_num_threads],
               [OMP_LIBS="-lgomp"], [])
fi
AS_IF([test "x$OMP_LIBS" == "x"], [AC_MSG_ERROR([could not find an openmp library])])

AC_CHECK_LIB([m], [cos],
             [], [AC_MSG_ERROR([could not find the libm library])])

AC_CHECK_LIB([gslcblas], [cblas_dgemm],
             [], [AC_MSG_ERROR([could not find the gslcblas library])])

AC_CHECK_LIB([gsl], [gsl_ran_dirichlet],
             [], [AC_MSG_ERROR([could not find the gsl library])])

AC_CHECK_LIB([quadmath], [quadmath_snprintf],
             [], [AC_MSG_ERROR([could not find the quadmath library])])

AC_CHECK_LIB([gfortran], [main],
             [], [AC_MSG_ERROR([could not find the gfortran library])])

AC_CHECK_LIB([arpack], [main],
             [], [AC_MSG_ERROR([could not find the arpack library])])

AC_F77_FUNC(sgemm)
AC_CHECK_LIB([openblas], [$sgemm],
             [], [AC_MSG_ERROR([could not find the openblas library])])

# Check that the Armadillo runtime library is installed if we need it
AS_IF([test "x$static_compile" = "xfalse" && test "x$arma_wrapper" = "xtrue"],
      [AC_CHECK_LIB([armadillo], [main],
                    [], AC_MSG_ERROR([could not find the library armadillo]))])

# Add openblas and openmp libraries
LIBS="-lopenblas $LIBS $OMP_LIBS"

# Set config headers and makefiles
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
