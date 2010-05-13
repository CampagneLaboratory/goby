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

import os
import stat
import struct
import StringIO
from gzip import GzipFile

import Alignments_pb2

# Python implementation of Java's DataInputStream.readUTF method.
def read_utf(fd):
    length = struct.unpack('>H', fd.read(2))[0]
    return fd.read(length).decode('utf8')

# Python implementation of Java's DataInputStream.readInt method.
def read_int(fd):
    return struct.unpack('>I', fd.read(4))[0]

#
# Class to parse goby "chunked" data
#
class MessageChunksReader():
    # verbose messages
    verbose = False

    # basename for this alignment
    basename = None

    # the name of the entries file for this alignment
    entries_filename = None

    # the size of the entries file (in bytes)
    entries_filesize = None

    # the current index into the entries file
    entries_fileindex = None;

    # length of the delimiter tag (in bytes)
    DELIMITER_LENGTH = 8;

    def __init__(self, basename, verbose = False):
        self.verbose = verbose

        # store the basename
        self.basename = basename

        # read the entries file
        self.entries_filename = basename + ".entries"
        if self.verbose:
            print "Reading entries from", self.entries_filename

        self.entries_filesize = os.stat(self.entries_filename)[stat.ST_SIZE]
        self.entries_fileindex = 0

    #
    # Return next chunk of bytes from the file
    #
    def next(self):
        # stop iteration if the file index has reached the end of file
        if self.entries_fileindex >= self.entries_filesize:
            raise StopIteration

        # open the entries file
        f = open(self.entries_filename, "rb")
        try:
            #  position to point just after the next delimiter
            self.entries_fileindex += self.DELIMITER_LENGTH
            if self.verbose:
                print "seeking to position", self.entries_fileindex
            f.seek(self.entries_fileindex, os.SEEK_SET)

            # get the number of bytes expected in the next chunk
            num_bytes = read_int(f)
            if self.verbose:
                print "expecting", num_bytes, "in next chunk"

            # if there are no more bytes, we're done
            if (num_bytes == 0):
                raise StopIteration
            else:
                buf = f.read(num_bytes)
                buf = GzipFile("", "rb", 0, StringIO.StringIO(buf)).read()
                return buf
        finally:
            self.entries_fileindex = f.tell()
            f.close()

    def __iter__(self):
        return self

    def __str__(self):
        return self.basename
