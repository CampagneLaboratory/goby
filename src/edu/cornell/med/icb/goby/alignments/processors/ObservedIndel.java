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

/**
 * @author Fabien Campagne
 *         Date: May 14, 2011
 *         Time: 10:50:34 AM
 */
public class ObservedIndel {
    /**
     * Start position is zero-based.
     */
    int startPosition;
    /**
     * End position is zero-based.
     */
    int endPosition;
    String from;
    String to;

    /**
     * Return the length of the indel, in bases (e.g., --- has a length of 3).
     * @return
     */
    public int length() {return endPosition-startPosition;}

    /**
     * Construct an indel observation.
     * @param startPosition The position where the indel starts, zero-based, position of the base at the left of the first gap.
     * @param endPosition   The position where the indel ends, zero-based, position of the base at the right of the first gap.
     * @param from          Bases in the reference
     * @param to            Bases in the read
     */
    public ObservedIndel(int startPosition, int endPosition, String from, String to) {
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.from=from;
        this.to=to;
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

    public boolean isReadInsertion() {
        return to.contains("-");
    }
   
}
