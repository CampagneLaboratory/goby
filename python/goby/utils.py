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

""" General purpose utiliies used by Goby. """

import re

commify_regex = re.compile(r'^(-?\d+)(\d{3})')

def commify(num, separator=','):
    """
    commify(num, separator) ->  string
    
    Return a string representing the number num with separator inserted
    for every power of 1000.   Separator defaults to a comma.
    E.g., commify(1234567) ->  '1,234,567'
    """

    num = str(num)  # just in case we were passed a numeric value
    more_to_do = 1
    while more_to_do:
        (num, more_to_do) = commify_regex.subn(r'\1%s\2' % separator,num)
    return num
