#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.59])
AC_INIT([finestructure], [0.0.4], [dan.lawson@bristol.ac.uk])
AC_CONFIG_AUX_DIR(config)
#AC_CONFIG_SRCDIR([fines.h])
#AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

AC_PROG_CXX
# Checks for libraries.
# find gsl, wherever it's hiding on the system
#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm])
#AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AX_PATH_GSL(,
			,
			 AC_MSG_ERROR(could not find required version of GSL))

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([floor pow])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
