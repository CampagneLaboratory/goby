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

""" Contains classes that can parse binary alignment data
stored in the Goby "compact" format.
"""

import gzip
import sys

import Alignments_pb2
from MessageChunks import MessageChunksReader

# Java properties - http://pypi.python.org/pypi/pyjavaproperties/
from pyjavaproperties import Properties

def get_basename(filename):
    """ Return the basename corresponding to the input alignment filename.
    Note that if the filename does have the extension known to be a
    compact alignment the returned value is the original filename.
    """
    for ext in [".entries", ".header", ".tmh", ".stats", ".counts", ".index"]:
        if filename.endswith(ext):
            return filename[:-len(ext)]
    return filename

class AlignmentReader(object):
    """ Reads alignments written in the Goby "compact" format.
    The AlignmentReader is actually an iterator over indiviual
    alignment entries stored in the file.
    """

    basename = None
    """ basename for this alignment """

    statistics = Properties()
    """ statistics from this alignment """

    header = Alignments_pb2.AlignmentHeader()
    """ alignment header """

    entries_reader = None
    """ reader for the alignment entries (interally stored in chunks) """

    entries = []
    """ Current chunk of alignment entries """

    current_entry_index = 0
    """ current entry index """

    def __init__(self, basename, verbose = False):
        """ Initialize the AlignmentReader using the
        basename of the alignment files.  Goby alignments
        are made up of files ending with:

        ".entries", ".header", ".tmh", ".stats", ".counts", ".index"

        The basename of this alignment is the full path to
        one of these files including everything but the
        extension.  If a basename is passed in with one of
        these extensions, it will be stripped.
        """

        # store the basename
        self.basename = get_basename(basename);

        # read the "stats" file
        stats_filename = self.basename + ".stats"
        if verbose:
            print "Reading properties from", stats_filename
        try:
            self.statistics.load(open(stats_filename))
        except IOError, err:
            # the "stats" file is optional
            print >>sys.stderr, "Could not read alignment statistics from", stats_filename
            print >>sys.stderr, str(err)
            pass

        # read the header
        header_filename = self.basename + ".header"
        if verbose:
            print "Reading header from", header_filename
        f = gzip.open(header_filename, "rb")
        self.header.ParseFromString(f.read())
        f.close()

        # open the entries
        self.entries_reader = AlignmentCollectionReader(self.basename, verbose)

    def next(self):
        """ Return next alignment entry from the entries file
        """
        # is it time to get the next chunk from the entries file?
        if not self.entries or self.current_entry_index >= len(self.entries):
            self.entries = self.entries_reader.next().alignment_entries
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

    basename = None
    """ basename for this alignment """

    tmh = Alignments_pb2.AlignmentTooManyHits()
    """ too many hits """

    queryindex_to_numhits = dict()
    """ query index to number of hits """

    queryindex_to_depth = dict()
    """ query index to depth/length of match. """

    def __init__(self, basename, verbose = False):
        """ Initialize the TooManyHitsReader using the basename
        of the alignment files.
        """

        # store the basename
        self.basename = get_basename(basename)

        # read the "too many hits" info
        tmh_filename = self.basename + ".tmh"
        if verbose:
            print "Reading too many hits info from", tmh_filename

        try:
            f = gzip.open(tmh_filename, "rb")
            self.tmh.ParseFromString(f.read())
            f.close()
            for hit in self.tmh.hits:
                self.queryindex_to_numhits[hit.query_index] = hit.at_least_number_of_hits
                if hit.HasField("length_of_match"):
                    self.queryindex_to_depth[hit.query_index] = hit.length_of_match
        except IOError, err:
            # the "Too Many Hits" file is optional
            print >>sys.stderr, "Could not read \"Too Many Hits\" information from", tmh_filename
            print >>sys.stderr, "Assuming no queries have too many hits."
            print >>sys.stderr, str(err)
            pass

    def __str__(self):
        return self.basename

class AlignmentCollectionReader(MessageChunksReader):
    """ Iterator for alignment collections within the alignment entries."""

    def __init__(self, basename, verbose = False):
        MessageChunksReader.__init__(self, get_basename(basename) + ".entries", verbose)

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
