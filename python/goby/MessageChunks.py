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

""" Support for reading Goby chunked "compact" files. """

import os
import struct
import StringIO
from gzip import GzipFile

def read_int(fd):
    """ Python implementation of Java's DataInputStream.readInt method.
    """
    buf = fd.read(4)
    if len(buf) == 4:
        length = struct.unpack('>I', buf)[0]
    else:
        length = 0
    return length

class MessageChunksReader(object):
    """ Class to parse a file that contains goby "chunked" data.
    The MessageChunksReader is actually an iterator over indiviual
    entries stored in the compressed file.
    """
 
    verbose = False
    """ print verbose messages """

    filename = None
    """ the name of the chunked file """

    fileobject = None
    """ handle to the actual chunked file """

    fileindex = None
    """ the current index into the chunked file """

    DELIMITER_LENGTH = 8
    """ length of the delimiter tag (in bytes) """

    def __init__(self, filename, verbose = False):
        """ Initialize the MessageChunksReader using the name
        of the of the compact file (reads, alignment entries, etc.)
        """

        self.verbose = verbose

        if self.verbose:
            print "Reading data from", filename

        # store the file information
        self.filename = filename

        # open the file
        self.fileobject = open(self.filename, "rb")

        # and set the current file pointer to the beginning of the file
        self.fileindex = 0

    def __del__(self):
        """ Close the underlying file when the object is destroyed """
        if self.verbose:
            print "closing file:", self.filename
        self.fileobject.close()

    def next(self):
        """ Return next chunk of bytes from the file. """
        try:
            #  position to point just after the next delimiter
            self.fileindex += self.DELIMITER_LENGTH
            if self.verbose:
                print "seeking to position", self.fileindex
            self.fileobject.seek(self.DELIMITER_LENGTH, os.SEEK_CUR)

            # get the number of bytes expected in the next chunk
            num_bytes = read_int(self.fileobject)
            if self.verbose:
                print "expecting", num_bytes, "in next chunk"

            # if there are no more bytes, we're done
            if (num_bytes == 0):
                raise StopIteration
            else:
                # each chunk is compressed
                buf = self.fileobject.read(num_bytes)
                buf = GzipFile("", "rb", 0, StringIO.StringIO(buf)).read()
                return buf
        finally:
            self.fileindex = self.fileobject.tell()

    def __iter__(self):
        """ The message chunks reader is actually an iterator over chunks
        in the compact file.  Each iteration will return a chunk of bytes raw
        bytes which can then be processed """
        return self

    def __str__(self):
        return self.filename
