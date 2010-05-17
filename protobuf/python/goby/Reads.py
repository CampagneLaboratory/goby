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
import Reads_pb2
from MessageChunks import MessageChunksReader

#
# Reads sequence "reads" written in the Goby "compact" format
#
class ReadsReader():
    # name of this file
    filename = None

    # reader for the alignment entries (interally stored in chunks)
    entries_reader = None

    # Current chunk of alignment entries
    entries = None

    # current entry index
    current_entry_index = 0
    
    def __init__(self, filename, verbose = False):
        # store the filename
        self.filename = filename

        # open the entries
        self.entries_reader = ReadsCollectionReader(filename)

    #
    # Return next read entry from the file
    #
    def next(self):
        # is it time to get the next chunk from the file?
        if self.entries is None or self.current_entry_index >= len(self.entries):
            self.entries = self.entries_reader.next().reads
            self.current_entry_index = 0

        entry = self.entries[self.current_entry_index]
        self.current_entry_index += 1
        return entry

    def __iter__(self):
        return self

    def __str__(self):
        return self.filename

#
# Iterator for read collections
#
class ReadsCollectionReader(MessageChunksReader):
    def __init__(self, filename, verbose = False):
        MessageChunksReader.__init__(self, filename, verbose)

    #
    # Return next alignment collection from the entries file
    #
    def next(self):
        buf = MessageChunksReader.next(self)
        collection = Reads_pb2.ReadCollection()
        collection.ParseFromString(buf)
        return collection

    def __iter__(self):
        return self

    def __str__(self):
        return self.filename
