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
import sys

from AlignmentReader import AlignmentReader

verbose = False

def usage():
    print "usage:", sys.argv[0], "[-h|--help] [-v|--verbose] <basename>" 

def main():
    global verbose

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
    print "Compact Alignment basename =", basename

    reader = AlignmentReader(basename, verbose)
    header = reader.header

    print "Info from header:"
    print "Number of target sequences =", header.number_of_targets
    print "Number of target length entries =", len(header.target_length)
    print "Min target length =", min(header.target_length)
    print "Max target length =", max(header.target_length)
    print "Mean target length =", sum(header.target_length) / len(header.target_length)

    print "Number of query sequences =", header.number_of_queries
    print "Number of query length entries =", len(header.query_length)
#    print "Min query length =", min(header.query_length)
#    print "Max query length =", max(header.query_length)
#    print "Mean query length =", sum(header.query_length) / len(header.query_length)

    print "Constant query lengths =", header.constantQueryLength
    print "Has query identifiers =", len(header.query_name_mapping.mappings) > 0
    print "Has target identifiers =", len(header.target_name_mapping.mappings) > 0

if __name__ == "__main__":
    main()
        
