#!/usr/bin/env python

#
# Copyright (C) 2010 Institute for Computational Biomedicine,
#                    Weill Medical College of Cornell University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import getopt
import os
import struct
import sys
import StringIO
from gzip import GzipFile

import Alignments_pb2

def read_utf(fd):
    "Python implementation of Java's DataInputStream.readUTF method."
    length = struct.unpack('>H', fd.read(2))[0]
    return fd.read(length).decode('utf8')

def read_int(fd):
    "Python implementation of Java's DataInputStream.readInt method."
    return struct.unpack('>I', fd.read(4))[0]

def usage():
    print "usage:", sys.argv[0], "[-h|--help] [-v|--verbose] <basename>" 

def main():
    DELIMITER_CONTENT = 0xFF;
    DELIMITER_LENGTH = 8;

    verbose = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hv", ["help", "verbose"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(1)

    # Collect options
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-v", "--verbose"):
            verbose = True
        
    if len(args) != 1:
        usage()
        sys.exit(2)

    basename = args[0]
    if verbose:
        print "Compact Alignment basename =", basename

    # read the entries file
    entries_filename = basename + ".entries"
    if verbose:
        print "Reading entries from", entries_filename
    f = open(entries_filename, "rb")
    try:
        collection = Alignments_pb2.AlignmentCollection()
        num_bytes = None
        while num_bytes != 0:
            f.seek(DELIMITER_LENGTH, os.SEEK_CUR)
            num_bytes = read_int(f)
            if (num_bytes != 0):
                buf = f.read(num_bytes)
                buf = GzipFile("", "rb", 0, StringIO.StringIO(buf)).read()
                collection.ParseFromString(buf)
                print collection
        
    finally:
        f.close()

if __name__ == "__main__":
    main()
        
