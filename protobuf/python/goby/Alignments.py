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

import gzip
import Alignments_pb2
from MessageChunks import MessageChunksReader

# Java properties - http://pypi.python.org/pypi/pyjavaproperties/
from pyjavaproperties import Properties

""" Contains classes that can parse binary alignment data
stored in the Goby "compact" format.
"""

def get_basename(filename):
    """ Return the basename corresponding to the input alignment filename.
    Note that if the filename does have the extension known to be a
    compact alignemt the returned value is the original filename.
    """
    for ext in [".entries", ".header", ".tmh", ".stats", ".counts"]:
        if filename.endswith(ext):
            return filename[:-len(ext)]
    return filename

class AlignmentReader(object):
    """ Reads alignments written in the Goby "compact" format.
    The AlignmentReader is actually an iterator over indiviual
    alignment entries stored in the file.
    """

    # basename for this alignment
    basename = None

    # statistics from this alignment
    statistics = Properties()

    # alignment header
    header = Alignments_pb2.AlignmentHeader()

    # reader for the alignment entries (interally stored in chunks)
    entries_reader = None

    # Current chunk of alignment entries
    entries = []

    # current entry index
    current_entry_index = 0
    
    def __init__(self, basename, verbose = False):
        """ Initialize the AlignmentReader using the
        basename of the alignment files.  Goby alignments
        are made up of files ending with:

        ".entries", ".header", ".tmh", ".stats",

        The basename of this alignment is the full path to
        one of these files including everything but the
        extension.
        """

        # store the basename
        self.basename = basename

        # read the "stats" file
        stats_filename = basename + ".stats"
        if verbose:
            print "Reading properties from", stats_filename
        try:
            self.statistics.load(open(stats_filename))
        except IOError, err:
            print str(err)
            pass

        # read the header
        header_filename = basename + ".header"
        if verbose:
            print "Reading header from", header_filename
        f = gzip.open(header_filename, "rb")
        self.header.ParseFromString(f.read())
        f.close()

        # open the entries
        self.entries_reader = AlignmentCollectionReader(basename)

    def next(self):
        """ Return next alignment entry from the entries file
        """
        # is it time to get the next chunk from the entries file?
        if not self.entries or self.current_entry_index >= len(self.entries):
            self.entries = self.entries_reader.next().alignmentEntries
            self.current_entry_index = 0

            # If the entries came back empty, we're done
            if not self.entries:
                raise StopIteration

        entry = self.entries[self.current_entry_index]
        self.current_entry_index += 1
        return entry

    def __iter__(self):
        return self

    def __str__(self):
        return self.basename

class TooManyHitsReader(object):
    """ Reads "Too Many Hits information from alignments
    written in the Goby "compact" format
    """

    # basename for this alignment
    basename = None

    # too many hits
    tmh = Alignments_pb2.AlignmentTooManyHits()

    # query index to number of hits
    queryindex_to_numhits = dict()

    # query index to depth/length of match.
    queryindex_to_depth = dict()

    def __init__(self, basename, verbose = False):
        """ Initialize the TooManyHitsReader using the basename
        of the alignment files.
        """

        # store the basename
        self.basename = basename

        # read the "too many hits" info
        tmh_filename = basename + ".tmh"
        if verbose:
            print "Reading too many hits info from", tmh_filename

        try:
            f = open(tmh_filename, "rb")
            self.tmh.ParseFromString(f.read())
            f.close()

            for hit in self.tmh.hits:
                self.queryindex_to_numhits[hit.query_index] = hit.at_least_number_of_hits
                if hit.HasField("length_of_match"):
                    self.queryindex_to_depth[hit.query_index] = hit.length_of_match
        except IOError, err:
            print str(err)
            pass

    def __str__(self):
        return self.basename

class AlignmentCollectionReader(MessageChunksReader):
    """ Iterator for alignment collections within the alignment entries."""

    def __init__(self, basename, verbose = False):
        MessageChunksReader.__init__(self, basename + ".entries", verbose)

    def next(self):
        """ Return next alignment collection from the entries file."""
        buf = MessageChunksReader.next(self)
        collection = Alignments_pb2.AlignmentCollection()
        collection.ParseFromString(buf)
        return collection

    def __iter__(self):
        return self

    def __str__(self):
        return self.basename
