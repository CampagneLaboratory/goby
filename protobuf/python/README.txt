This directory contains the Python API for reading binary data files
created using the Goby next-gen data management framework.

Normally, this directory comes as part of the complete Goby package,
available from:

  http://goby.campagnelab.org/

The complete package includes the Java source code.  If you downloaded
this package from PyPI or some other Python-specific source, you may
have received only the Python part of the code.  In this case, you
will need to obtain the Protocol Compiler from some other source
before you can use this package.

Development Warning
===================

The Goby Python libraries are not as mature as the Java
implementation.  It may be more buggy and is not intended to provide
the complete set of features that are found in the Java version.

Installation
============

1) Make sure you have Python 2.5 or newer.  If in doubt, run:

     $ python -V

2) Download and install the prerequsite python packages:

2a) Google Protocol Buffers

    - Avaiable from http://code.google.com/p/protobuf/ or PyPI

2b) pyjavaproperties - Python replacement for java.util.Properties

   - Available from http://pypi.python.org/pypi/pyjavaproperties

3) Install the Goby package:

     $ python setup.py install

   This step may require superuser privileges.

Usage
=====

The complete documentation for Goby is available via the
web at:

  http://goby.campagnelab.org/
