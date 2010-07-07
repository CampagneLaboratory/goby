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

""" Contains classes that can parse binary sequence data
stored in the Goby "compact" format.
"""

import gzip
import Reads_pb2
from MessageChunks import MessageChunksReader

class ReadsReader(object):
    """ Reads sequence "reads" written in the Goby "compact" format.
    The ReadsReader is actually an iterator over indiviual
    read entries stored in the file.
    """

    filename = None
    """ name of the compact reads file """

    entries_reader = None
    """ reader for the read entries (interally stored in chunks) """

    entries = []
    """ Current chunk of read entries """

    current_entry_index = 0
    """ current entry index """
    
    def __init__(self, filename, verbose = False):
        """ Initialize the ReadsReader using the name
        of the compact file.  Goby compact reads have
        the extension ".compact-reads"
        """

        # store the filename
        self.filename = filename

        # open the entries
        self.entries_reader = ReadsCollectionReader(filename, verbose)

    def next(self):
        """ Return next read entry from the file.
        """

        # is it time to get the next chunk from the file?
        if not self.entries or self.current_entry_index >= len(self.entries):
            self.entries = self.entries_reader.next().reads
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
        return self.filename

class ReadsCollectionReader(MessageChunksReader):
    """ Iterator for read collections within a compact
    reads file
    """

    def __init__(self, filename, verbose = False):
        MessageChunksReader.__init__(self, filename, verbose)

    def next(self):
        """ Return next read collection from the
        compact reads file
        """
        buf = MessageChunksReader.next(self)
        collection = Reads_pb2.ReadCollection()
        collection.ParseFromString(buf)
        return collection

    def __iter__(self):
        return self

    def __str__(self):
        return self.filename
