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

from goby.Reads import ReadsReader
from goby.utils import commify

def usage():
    print "usage:", sys.argv[0], "[-h|--help] [-v|--verbose] <fileame>" 

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

    filename = args[0]
    print "Compact reads filename = %s" % filename
    print

    filesize = os.stat(filename)[stat.ST_SIZE]
    number_of_entries = 0
    number_of_identifiers = 0
    number_of_descriptions = 0
    number_of_sequences = 0
    number_of_sequence_pairs = 0
    number_of_quality_scores = 0
    number_of_quality_score_pairs = 0
    min_read_length = sys.maxint
    max_read_length = -sys.maxint - 1
    total_read_length = 0
    total_read_length_pair = 0

    reads_reader = ReadsReader(filename, verbose)
    for entry in reads_reader:
        number_of_entries += 1
        read_length = entry.read_length
        total_read_length += read_length
        total_read_length_pair += entry.read_length_pair

        if entry.HasField("read_identifier"):
            number_of_identifiers += 1
        if entry.HasField("description"):
            number_of_descriptions += 1
        if entry.HasField("sequence"):
            number_of_sequences += 1
        if entry.HasField("sequence_pair"):
            number_of_sequences_pairs += 1
        if entry.HasField("quality_scores"):
            number_of_quality_scores += 1
        if entry.HasField("quality_scores_pair"):
            number_of_quality_score_pairs += 1

        min_read_length = min(min_read_length, read_length)
        max_read_length = max(max_read_length, read_length)

    print "Average bytes per entry: %s" % commify(filesize / float(number_of_entries))
    print "Average bytes per base:  %s" % commify(filesize / float(total_read_length))
    print "Has identifiers = %s (%s)" % (number_of_identifiers > 0, commify(number_of_identifiers))
    print "Has descriptions = %s (%s)" % (number_of_descriptions > 0, commify(number_of_descriptions))
    print "Has sequences = %s (%s)" % (number_of_sequences > 0, commify(number_of_sequences))
    print "Has sequence pairs = %s (%s)" % (number_of_sequence_pairs > 0, commify(number_of_sequence_pairs))
    print "Has quality scores = %s (%s)" % (number_of_quality_scores > 0, commify(number_of_quality_scores))
    print "Has quality score pairs = %s (%s)" % (number_of_quality_score_pairs > 0, commify(number_of_quality_score_pairs))
    print "Number of entries = %s" % commify(number_of_entries)
    print "Min read length = %s" % commify(min_read_length)
    print "Max read length = %s" % commify(max_read_length)
    print "Avg read length = %s" % commify(total_read_length / float(number_of_entries))
    print "Avg read pair length = %s" % commify(total_read_length_pair / float(number_of_entries))

if __name__ == "__main__":
    main()
