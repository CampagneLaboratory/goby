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

package edu.cornell.med.icb.goby.alignments.processors;

import java.text.Normalizer;

/**
 * @author Fabien Campagne
 *         Date: May 14, 2011
 *         Time: 10:50:34 AM
 */
public class ObservedIndel {
    /**
     * Start position is zero-based.
     */
    final int startPosition;
    /**
     * End position is zero-based.
     */
    final int endPosition;
    final String from;
    final String to;

    /**
     * Return the length of the indel, in bases (e.g., --- has a length of 3).
     *
     * @return
     */
    public int length() {
        return endPosition - startPosition;
    }

    /**
     * Construct an indel observation.
     *
     * @param startPosition The position where the indel starts, zero-based, position of the base at the left of the first gap.
     * @param endPosition   The position where the indel ends, zero-based, position of the base at the right of the first gap.
     * @param from          Bases in the reference
     * @param to            Bases in the read
     */
    public ObservedIndel(final int startPosition, final int endPosition, final String from, final String to) {
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.from = from;
        this.to = to;
    }

    public int getStart() {
        return startPosition;
    }

    public int getEnd() {
        return endPosition;
    }

    public String from() {
        return from;
    }

    public String to() {
        return to;
    }
    public final boolean isReadInsertion() {
        return to.contains("-");
    }

    @Override
    public int hashCode() {
        return startPosition ^ endPosition ^ from.hashCode() ^ to.hashCode();
    }

    @Override
    public boolean equals(final Object o) {
        if (!(o instanceof ObservedIndel)) {
            return false;
        }
        final ObservedIndel other = (ObservedIndel) o;
        return other.startPosition == startPosition &&
                other.endPosition == endPosition &&
                other.from().equals(from) && other.to.equals(to);
    }

    @Override
    public String toString() {
        return String.format("%s/%s %d-%d",from,to ,startPosition,endPosition);
    }
}
