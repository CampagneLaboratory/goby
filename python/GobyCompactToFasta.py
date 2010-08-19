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
import struct
import sys
import textwrap

from goby.Reads import ReadsReader

def usage():
    print "usage:", sys.argv[0], "[-h|--help] [-v|--verbose] [-f|--format <fasta|fastq>] [-o|--output <output-filename>] <filename>" 

def main():
    verbose = False
    format = "fasta"
    output = sys.stdout

    fake_quality_score = 40

    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:o:hv", ["format=", "output=", "help", "verbose"])
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

    if format != "fasta" and format != "fastq":
        print >> sys.stderr, "Format", format, "is not supported"
        usage()
        sys.exit(3)

    if len(args) != 1:
        usage()
        sys.exit(2)

    filename = args[0]
    if verbose:
        print "Processing file =", filename

    if format == "fasta":
        new_entry_character = ">"
    else:
        new_entry_character = "@"

    reads_reader = ReadsReader(filename, verbose)
    for entry in reads_reader:
        # if the entry has a description, use it otherwise use the index
        if (entry.HasField("description")):
            description = entry.description
        else:
            description = entry.read_index
        print >> output, "%c%s" % (new_entry_character, description)
        for line in textwrap.wrap(entry.sequence, 60):
            print >> output, line

        if format == "fastq":
            print >> output, "+"
            quality_string = ""
            if entry.HasField("qualityScores"):
                for quality_score in entry.quality_scores:
                    quality_string += chr(struct.unpack("B", quality_score)[0] + 64)  # Assume Illumina encoding
                
            else:
                # fill the string with a constant quality code
                quality_string = quality_string.ljust(len(entry.sequence), chr(fake_quality_score + 64))

            for line in textwrap.wrap(quality_string, 60):
                print >> output, line
                
if __name__ == "__main__":
    main()
