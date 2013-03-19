/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.goby.alignments;

/**
 * Describes a location on a reference sequence.
 *
 * @author Fabien Campagne
 *         Date: Dec 14, 2010
 *         Time: 6:00:39 PM
 */
public class ReferenceLocation implements Comparable {

    public final int targetIndex;
    public int position;
    /**
     * The amount of compressed data since the previous location.
     */

    public long compressedByteAmountSincePreviousLocation;

    public ReferenceLocation(int referenceIndex, int position) {
        this.targetIndex = referenceIndex;
        this.position = position;
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof ReferenceLocation)) {
            return false;
        }
        final ReferenceLocation other = (ReferenceLocation) o;
        return targetIndex == other.targetIndex && position == other.position;
    }

    @Override
    public int hashCode() {
        return targetIndex ^ position;
    }

    public int compareTo(Object o) {
        if (!(o instanceof ReferenceLocation)) {
            return -1;
        }

        ReferenceLocation other = (ReferenceLocation) o;
        if (other.targetIndex == targetIndex) {

            return position - other.position;

        } else {
            return targetIndex - other.targetIndex;
        }
    }
}
