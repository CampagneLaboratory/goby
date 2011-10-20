/*
 * Copyright (C) 2009-2011 Institute for Computational Biomedicine,
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

package edu.cornell.med.icb.goby.algorithmic.data.ranges;

/**
 * Represents a range (start-end).
 * @author Fabien Campagne
 *         Date: 10/19/11
 *         Time: 11:19 AM
 */
public class Range  implements Comparable<Range>{
    public int min=Integer.MIN_VALUE;
    public int max=Integer.MAX_VALUE;

    @Override
    public int hashCode() {
        return min ^ max;
    }

    @Override
    public boolean equals(Object o) {
        Range otherRange=(Range)o;
        return min==otherRange.min && max==otherRange.max;

    }

    public Range(int min, int max) {
        this.min=min;
        this.max=max;
    }

    public Range() {

    }

    @Override
    public String toString() {
        return String.format("min=%d max=%d ",min,max);
    }

    @Override
    public int compareTo(final Range range) {
        return min - range.min ;
    }
}
