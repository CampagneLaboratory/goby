This directory contains the Python API for reading binary data files
created using the Goby next-gen data management framework.

Normally, this directory comes as part of the complete Goby package,
available from:

  http://goby.campagnelab.org/

The complete package includes the Java source code.  If you downloaded
this package from PyPI or some other Python-specific source, you may
have received only the Python part of the code.

Development Warning
===================

The Goby Python libraries are not as mature as the Java
implementation.  It may be more buggy and is not intended to provide
the complete set of features that are found in the Java version.

Installation
============

1) Make sure you have Python 2.5 or newer.  If in doubt, run:

     $ python -V

2) Download and install the prerequisite python packages:

a) Protocol Buffers

   Available from http://code.google.com/p/protobuf/ or PyPI

b) pyjavaproperties - Python replacement for java.util.Properties

   Available from http://pypi.python.org/pypi/pyjavaproperties

3) Install the Goby package:

     $ python setup.py install

   This step may require superuser privileges.

Usage
=====

Example scripts are provided to demonstrate how to access the content
of Goby files in Python.

- Here is how to scan a Goby alignment file:

  GobyAlignmentStats.py basename

(The files basename.entries and basename.header must exit.)

- The next command will print the content of an alignment file as text:

GobyAlignmentToText.py basename

- The next command will convert a compact reads file to fasta format:

GobyCompactToFasta.py file.compact-reads

- The next command will print statistics about the content of a
  compact reads file: 

GobyReadsStats.py file.compact-reads

Documentation
=============

The complete documentation for Goby is available online at:

  http://goby.campagnelab.org/

