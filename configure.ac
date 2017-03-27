# Initialisation
AC_INIT([fastspar], [0.0.3], [step.wat@gmail.com])
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
AS_IF([test "x$static_compile" = "xfalse" && test "x$arma_wrapper" = "xtrue"], [AC_CHECK_LIB(armadillo, main, [LIBS="$LIBS -larmadillo"], AC_MSG_ERROR(could not find the library Armadillo))])

# Set conditionals to use in the Makefile.am
AM_CONDITIONAL([ARMA_DONT_USE_WRAPPER], [test "x$static_compile" = "xtrue" || test "x$arma_wrapper" = "xfalse"])
AM_CONDITIONAL([STATIC_COMPILE], [test "x$static_compile" = "xtrue"])


# Check for c++ and fortran compilers
AC_PROG_CXX
AC_PROG_F77


# Use built-in macros to check for libraries
AC_CHECK_LIB(gomp, main, [LIBS="$LIBS -lgomp"], AC_MSG_ERROR(could not find the library OpenMP))
AC_CHECK_LIB(m, cos, [LIBS="$LIBS -lm"], AC_MSG_ERROR(could not find the library LibM))
AC_CHECK_LIB(gslcblas, cblas_dgemm, [LIBS="-lgslcblas $LIBS"], AC_MSG_ERROR(could not find the library GSL (gslcblas)))
AC_CHECK_LIB(gsl, gsl_ran_dirichlet, [LIBS="-lgsl $LIBS"], AC_MSG_ERROR(could not find the library GSL))
AC_CHECK_LIB(quadmath, quadmath_snprintf, [LIBS="-lquadmath $LIBS"], AC_MSG_ERROR(could not find the library quadmath))
AC_CHECK_LIB(gfortran, main, [LIBS="-lgfortran $LIBS"], AC_MSG_ERROR(could not find the library gfortran))
AC_CHECK_LIB(openblas, main, [LIBS="-lopenblas $LIBS"], AC_MSG_ERROR(could not find the library OpenBLAS))

# Set config header output
AC_CONFIG_HEADERS([config.h])

# Declare Makefile outputs
AC_CONFIG_FILES([Makefile src/Makefile])

# Disables rebuild rules (fix issue when cloning github repo and Makefile.am, configure.ac, and aclocal.m4 have borked datastamps)
AM_MAINTAINER_MODE([disable])
AC_OUTPUT
