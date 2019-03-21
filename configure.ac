# Initialisation
AC_INIT([fastspar], [0.0.10], [s.watts2@student.unimelb.edu.au])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

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

# Check that the Armadillo runtime library is installed if we need it
AS_IF([test "x$static_compile" = "xfalse" && test "x$arma_wrapper" = "xtrue"],
      [AC_CHECK_LIB([armadillo], [main],
                    [], AC_MSG_ERROR([could not find the library Armadillo]))])

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
AC_CHECK_LIB([gomp], [main],
             [], [AC_MSG_ERROR([could not find the library OpenMP])])

AC_CHECK_LIB([m], [cos],
             [], [AC_MSG_ERROR([could not find the library LibM])])

AC_CHECK_LIB([gslcblas], [cblas_dgemm],
             [], [AC_MSG_ERROR([could not find the library GSL \(gslcblas\)])])

AC_CHECK_LIB([gsl], [gsl_ran_dirichlet],
             [], [AC_MSG_ERROR([could not find the library GSL])])

AC_CHECK_LIB([quadmath], [quadmath_snprintf],
             [], [AC_MSG_ERROR([could not find the library quadmath])])

AC_CHECK_LIB([gfortran], [main],
             [], [AC_MSG_ERROR([could not find the library gfortran])])

AX_BLAS([], [AC_MSG_ERROR([could not find the library BLAS])])

# Set libraries as specified by AX_BLAS documentation
LIBS="$BLAS_LIBS $LIBS $FLIBS"

# Set config headers and makefiles
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
