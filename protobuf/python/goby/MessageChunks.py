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

""" Support for reading Goby "compact" files.
"""

def read_utf(fd):
    """ Python implementation of Java's DataInputStream.readUTF method.
    """
    length = struct.unpack('>H', fd.read(2))[0]
    return fd.read(length).decode('utf8')

def read_int(fd):
    """ Python implementation of Java's DataInputStream.readInt method.
    """
    return struct.unpack('>I', fd.read(4))[0]

class MessageChunksReader(object):
    """ Class to parse a file that contains goby "chunked" data.
    The MessageChunksReader is actually an iterator over indiviual
    entries stored in the compressed file.
    """

    # verbose messages
    verbose = False

    # the name of the chunked file
    filename = None

    # the size of the chunked file (in bytes)
    filesize = None

    # the current index into the chunked file
    fileindex = None;

    # length of the delimiter tag (in bytes)
    DELIMITER_LENGTH = 8;

    def __init__(self, filename, verbose = False):
        self.verbose = verbose

        # read the entries file
        self.filename = filename
        if self.verbose:
            print "Reading data from", self.filename

        self.filesize = os.stat(self.filename)[stat.ST_SIZE]
        self.fileindex = 0

    def next(self):
        """ Return next chunk of bytes from the file. """

        # stop iteration if the file index has reached the end of file
        if self.fileindex >= self.filesize:
            raise StopIteration

        # open the entries file
        f = open(self.filename, "rb")
        try:
            #  position to point just after the next delimiter
            self.fileindex += self.DELIMITER_LENGTH
            if self.verbose:
                print "seeking to position", self.fileindex
            f.seek(self.fileindex, os.SEEK_SET)

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
            self.fileindex = f.tell()
            f.close()

    def __iter__(self):
        return self

    def __str__(self):
        return self.filename
