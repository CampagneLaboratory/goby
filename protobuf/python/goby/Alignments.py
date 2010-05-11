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

import gzip
import Alignments_pb2

# Java properties - http://pypi.python.org/pypi/pyjavaproperties/
from pyjavaproperties import Properties

#
# Reads alignments written in the Goby "compact" format
#
class AlignmentReader():
    # basename for this alignment
    basename = None

    # statistics from this alignment
    statistics = None

    # alignment header
    header = Alignments_pb2.AlignmentHeader()

    # too many hits
    tmh = Alignments_pb2.AlignmentTooManyHits()

    def __init__(self, basename, verbose = False):
        # store the basename
        self.basename = basename

        # read the "stats" file
        stats_filename = basename + ".stats"
        if verbose:
            print "Reading properties from", stats_filename

        self.statistics = Properties()
        self.statistics.load(open(stats_filename))

        # read the header (TODO: support old format?)
        f = gzip.open(basename + ".header", "rb")
        self.header.ParseFromString(f.read())
        f.close()

        # read the "too many hits" info
        f = open(basename + ".tmh", "rb")
        self.tmh.ParseFromString(f.read())
        f.close()

    def __str__(self):
        return self.basename
