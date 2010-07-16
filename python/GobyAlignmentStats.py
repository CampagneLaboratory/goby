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
import stat
import sys
import goby

from goby.Alignments import AlignmentReader, TooManyHitsReader
from goby.utils import commify

def usage():
    print "usage:", sys.argv[0], "[-h|--help] [-v|--verbose] <basename>"

def main():
    verbose = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hv", ["help", "verbose"])
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
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

    basename = goby.Alignments.get_basename(args[0])
    print "Compact Alignment basename =", basename

    alignment_reader = AlignmentReader(basename, verbose)
    header = alignment_reader.header
    tmh_reader = TooManyHitsReader(basename, verbose)
    tmh = tmh_reader.tmh
    entries_filesize = os.stat(basename + ".entries")[stat.ST_SIZE]
    print "Info from header:"
    print "Sorted:", header.sorted
    print "Indexed: ", header.indexed
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
    has_constant_query_length = header.HasField("constant_query_length")

    # special case if query lengths are constant to reduce storage
    if has_constant_query_length:
        query_length =  header.number_of_queries
        min_query_length = header.constant_query_length
        max_query_length = header.constant_query_length
        mean_query_length = float(header.constant_query_length)
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

    print "TMH: aligner threshold = %d" % tmh.aligner_threshold
    print "TMH: number of ambiguous matches = %s" % commify(len(tmh_reader.queryindex_to_numhits))
    print "TMH: %%ambiguous matches = %f %%" % (len(tmh_reader.queryindex_to_numhits) * 100.0 / header.number_of_queries)

    max_query_index = -1
    max_target_index = -1
    number_of_entries = 0
    number_of_logical_alignment_entries = 0
    total = 0
    average_score = 0.0
    number_of_variations = 0

    aligned_query_indices = set(tmh_reader.queryindex_to_numhits.keys())
    
    for entry in alignment_reader:
        number_of_entries += 1      # Across this file
        number_of_logical_alignment_entries += entry.multiplicity
        total += entry.query_aligned_length
        average_score += entry.score
        max_query_index = max(max_query_index, entry.query_index)
        max_target_index = max(max_target_index, entry.target_index)
        number_of_variations += len(entry.sequence_variations)
        aligned_query_indices.add(entry.query_index)

    number_of_query_sequences = max_query_index + 1
    number_of_target_sequences = max_target_index + 1
    average_score = average_score / float(number_of_logical_alignment_entries)

    print "num query indices = %s" % commify(number_of_query_sequences)
    print "num target indices = %s" % commify(number_of_target_sequences)
    print "Number of alignment entries = %s" % commify(number_of_logical_alignment_entries)
    print "Number of query indices that matched = %s" % commify(len(aligned_query_indices))
    print "Percent matched = %s %%" % commify(len(aligned_query_indices) * 100.0 / number_of_query_sequences)
    print "Avg query alignment length = %s" % commify(total / float(number_of_entries))
    print "Avg score alignment = %s" % commify(average_score)
    print "Avg number of variations per query sequence = %s" % commify(number_of_variations / float(number_of_query_sequences))
    print "Average bytes per entry = %s" % commify(entries_filesize / float(number_of_logical_alignment_entries))
    print

if __name__ == "__main__":
    main()
