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

package edu.cornell.med.icb.goby.algorithmic.data;

import it.unimi.dsi.lang.MutableString;

/**
 * Stores information about the span of equivalent indels.
 *
 * @author Fabien Campagne
 *         Date: 6/7/11
 *         Time: 6:15 PM
 */
public class EquivalentIndelRegion {
    /**
     * Index of the reference sequence on which the indel is observed.
     */
    public int referenceIndex;
    /**
     * Start position (zero-based) of the position after which the indel occurs.
     * ACT^{GG} bases   in this case, the start position of the GG insertion is 2 and its
     * 0123456  pos     end position is 3
     */
    public int startPosition;
    public int endPosition;
    public String from;
    public String to;
    public String flankLeft;
    public String flankRight;

    /**
     * Return the from bases, surrounded by flankLeft and flankRight bases.
     *
     * @return from bases in context of the flanking sequence.
     */
    public String fromInContext() {
        MutableString fromC = new MutableString();
        fromC.append(flankLeft);
        fromC.append(from);
        fromC.append(flankRight);
        return fromC.toString();
    }

    /**
     * Return the to bases, surrounded by flankLeft and flankRight bases.
     *
        * @return to bases in context of the flanking sequence.
     */
    public String toInContext() {
        MutableString toC = new MutableString();
        toC.append(flankLeft);
        toC.append(to);
        toC.append(flankRight);
        return toC.toString();
    }
}
