#                                               -*- Autoconf -*-
# Goby C++ libraries
# Process this file with autoconf to produce a configure script.
# In general, the safest way to proceed is to run ./autogen.sh

AC_PREREQ([2.59])
AC_INIT([Goby],[1.7-dev],[mas2062@med.cornell.edu])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_HEADER_STDC
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([memset])

AC_CONFIG_FILES([Makefile src/Makefile goby.pc])
AC_OUTPUT
