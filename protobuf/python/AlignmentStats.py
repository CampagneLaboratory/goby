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
import re
import sys

from goby.Alignments import AlignmentReader, TooManyHitsReader

verbose = False

commify_regex = re.compile(r'^(-?\d+)(\d{3})')

#
# commify(num, separator) ->  string
#
# Return a string representing the number num with separator inserted
# for every power of 1000.   Separator defaults to a comma.
# E.g., commify(1234567) ->  '1,234,567'
#
def commify(num, separator=','):
    num = str(num)  # just in case we were passed a numeric value
    more_to_do = 1
    while more_to_do:
        (num, more_to_do) = commify_regex.subn(r'\1%s\2' % separator,num)
    return num

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

    alignment_reader = AlignmentReader(basename, verbose)
    header = alignment_reader.header

    tmh_reader = TooManyHitsReader(basename, verbose)
    tmh = tmh_reader.tmh

    print "Info from header:"
    print "Number of target sequences = %s" % commify(header.number_of_targets)

    # target length stats
    target_length = len(header.target_length)
    if target_length > 0:
        min_target_length = min(header.target_length)
        max_target_length = max(header.target_length)
        mean_target_length =  sum(header.target_length) / float(target_length)
    else:
        min_target_length = 0
        max_target_length = 0
        mean_target_length = 0

    print "Number of target length entries = %s" % commify(target_length)
    print "Min target length = %s" % commify(min_target_length)
    print "Max target length = %s" % commify(max_target_length)
    print "Mean target length = %s" % commify(mean_target_length)
    print

    print "Number of query sequences = %s" % commify(header.number_of_queries)

    # query length stats
    has_constant_query_length = header.HasField("constantQueryLength")

    # special case if query lengths are constant to reduce storage
    if has_constant_query_length:
        query_length =  header.number_of_queries
        min_query_length = header.constantQueryLength
        max_query_length =  header.constantQueryLength
        mean_query_length =  header.constantQueryLength
    else:
        query_length = len(header.query_length)
        if query_length > 0:
            min_query_length = min(header.query_length)
            max_query_length = max(header.query_length)
            mean_query_length =  sum(header.query_length) / query_length
        else:
            min_query_length = 0
            max_query_length = 0
            mean_query_length = 0

    print "Number of query length entries = %s" % commify(query_length)
    print "Min query length = %s" % commify(min_query_length)
    print "Max query length = %s" % commify(max_query_length)
    print "Mean query length = %s" % commify(mean_query_length)

    print "Constant query lengths = %s" % has_constant_query_length
    print "Has query identifiers = %s" % (len(header.query_name_mapping.mappings) > 0)
    print "Has target identifiers = %s" % (len(header.target_name_mapping.mappings) > 0)
    print

    print "TMH: aligner threshold = %d" % tmh.alignerThreshold
    print "TMH: number of ambiguous matches = %s" % commify(len(tmh_reader.queryindex_to_numhits))
    print "TMH: %%ambiguous matches = %f %%" % (len(tmh_reader.queryindex_to_numhits) * 100.0 / header.number_of_queries)

if __name__ == "__main__":
    main()
        
