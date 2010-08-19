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

import csv
import goby
import getopt
import sys

from goby.Alignments import AlignmentReader, TooManyHitsReader

def usage():
    print "usage:", sys.argv[0], "[-h|--help] [-v|--verbose] [-f|--format <plain|sam>] [-o|--output <file>] <basename>" 

def main():
    format = "plain"
    output = sys.stdout
    verbose = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:ho:v", ["format=", "help", "output=", "verbose"])
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
        elif opt in ("-f", "--format"):
            format = arg
        elif opt in ("-o", "--output"):
            output = open(arg, "w")

    if format != "plain" and format != "sam":
        print >> sys.stderr, "Format", format, "is not supported"
        usage()
        sys.exit(3)

    if len(args) != 1:
        usage()
        sys.exit(2)

    basename = goby.Alignments.get_basename(args[0])
    if verbose:
        print "Compact Alignment basename =", basename

    reader = AlignmentReader(basename, verbose)
    header = reader.header
    writer = csv.writer(output, delimiter='\t')

    # write the header lines
    if format == "plain":
        writer.writerow(["queryId", "referenceId", "referenceLength", "numberOfIndels", "numberOfMismatches", "score", "startPosition", "alignmentLength", "matchesReverseStrand"])

    writer.writerow(["@HD","VN:1.0"])
    writer.writerow(["@PG", "Goby", "VN:python"])

    # write the target identifiers and lengths using the target mapping information
    for mapping in header.target_name_mapping.mappings:
        writer.writerow(["@SQ", "SN:" + mapping.name, "LN:" + str(header.target_length[mapping.index])])

    # write the actual alignment information
    for entry in reader:
        # use the query name if we have it
# TODO - write the query name
#        if header.query_name_mapping:
#            query_id = header.query_name_mapping[entry.query_index].name
#        else:
        query_id = entry.query_index

        target_id = header.target_name_mapping.mappings[entry.target_index].name
        target_length = header.target_length[entry.target_index]
        writer.writerow([query_id, target_id, target_length, entry.number_of_indels, entry.number_of_mismatches, entry.score, entry.position, entry.query_aligned_length, entry.matching_reverse_strand])

if __name__ == "__main__":
    main()
