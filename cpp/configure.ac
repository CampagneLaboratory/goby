#                                               -*- Autoconf -*-
# Goby C++ libraries
# Process this file with autoconf to produce a configure script.
# In general, the safest way to proceed is to run ./autogen.sh

AC_PREREQ(2.59)
AC_INIT([Goby],[2.2.1],[icb@med.cornell.edu])
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
AC_CHECK_SIZEOF(unsigned long)
AC_CHECK_SIZEOF(unsigned long long)

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
#
# Boost is OPTIONAL and will likely be removed in the near future
# For GMAP, you should make sure boost is NOT included
#
#AX_BOOST_BASE
#AX_BOOST_DATE_TIME
#AX_BOOST_FILESYSTEM
#AX_BOOST_PROGRAM_OPTIONS
#AX_BOOST_SYSTEM

# use pkg-config to make sure that protobuf is available
# see http://code.google.com/p/protobuf/
PKG_CHECK_MODULES([PROTOBUF], [protobuf >= 2.4.1]) 
# use pkg-config to make sure that PCRE is available
# see http://www.pcre.org/
PKG_CHECK_MODULES([LIBPCRECPP], [libpcrecpp >= 8.10]) 

AC_CONFIG_FILES([Makefile src/Makefile goby.pc])
AC_OUTPUT
