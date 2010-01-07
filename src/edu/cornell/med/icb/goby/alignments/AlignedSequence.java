/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
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

import it.unimi.dsi.lang.MutableString;

/**
 * Encodes information about one aligned sequence (support for MAF alignment format).
 *
 * @author Fabien Campagne
 *         Date: May 4, 2009
 *         Time: 3:54:17 PM
 */
public class AlignedSequence {
    public AlignedSequence() {
        sequenceIdentifier = new MutableString();
        alignment = new MutableString();
    }

    public MutableString sequenceIdentifier;
    public int alignedStart;
    public int alignedLength;
    public char strand;
    public int sequenceLength;
    public MutableString alignment;

}
