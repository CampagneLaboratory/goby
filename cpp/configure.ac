slib#                                               -*- Autoconf -*-
# Goby C++ libraries
# Process this file with autoconf to produce a configure script.
# In general, the safest way to proceed is to run ./autogen.sh

AC_PREREQ(2.59)
AC_INIT([Goby],[1.8-dev],[icb@med.cornell.edu])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE

AC_SUBST(version_info, [1:0:0])

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_HEADER_STDC
AC_C_CONST
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T

AC_CHECK_HEADERS([tr1/unordered_map])
AC_CHECK_HEADERS([tr1/unordered_set])
AC_CHECK_HEADERS([ext/hash_map])
AC_CHECK_HEADERS([ext/hash_set])

AC_CHECK_HEADERS([fcntl.h])
AC_CHECK_HEADERS([malloc.h])

# Checks for library functions.
AC_CHECK_FUNCS([getcwd])
AC_CHECK_FUNCS([memset])
AC_CHECK_FUNCS([putenv])
AC_CHECK_FUNCS([realpath])
AC_CHECK_FUNCS([setenv])
AC_CHECK_FUNCS([strerror])

AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT

AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK

# check for the boost libraries
AX_BOOST_BASE
AX_BOOST_PROGRAM_OPTIONS

# use pkg-config to make sure that protobuf is available
PKG_CHECK_MODULES([PROTOBUF], [protobuf >= 2.3.0]) 

AC_CONFIG_FILES([Makefile src/Makefile goby.pc])
AC_OUTPUT
