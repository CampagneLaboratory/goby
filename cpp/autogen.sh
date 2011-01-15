#!/bin/sh -e

test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.
#autoreconf --force -I /opt/local/share/aclocal --install --verbose "$srcdir"
# The following seems to work better on MacOs than autoreconf. Needs testing on other platforms as well.
aclocal -I /opt/local/share/aclocal
autoheader
automake
autoconf
test -n "$NOCONFIGURE" || "$srcdir/configure" "$@"
