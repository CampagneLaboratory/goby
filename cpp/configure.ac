#                                               -*- Autoconf -*-
# Goby C++ libraries
# Process this file with autoconf to produce a configure script.
# In general, the safest way to proceed is to run ./autogen.sh

AC_PREREQ(2.59)
AC_INIT([Goby],[1.7-dev],[mas2062@med.cornell.edu])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_HOST
AC_CANONICAL_TARGET
#AC_DISABLE_STATIC
AM_INIT_AUTOMAKE

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL

AC_SUBST(version_info, [0:0:0])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_HEADER_STDC
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([memset])

ACX_PTHREAD
#AC_CXX_STL_HASH

# use pkg-config to make sure that protobuf is available
PKG_CHECK_MODULES([PROTOBUF], [protobuf >= 2.3.0]) 

AC_CONFIG_FILES([Makefile src/Makefile goby.pc])
AC_OUTPUT
